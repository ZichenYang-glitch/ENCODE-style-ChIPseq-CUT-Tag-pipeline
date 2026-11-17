// @vitest-environment node
import { mkdtemp, mkdir, rm, writeFile } from 'node:fs/promises';
import { tmpdir } from 'node:os';
import path from 'node:path';
import { createRequire } from 'node:module';
import { build, type Plugin } from 'vite';
import { expect, it } from 'vitest';
import config from '../../vite.config';

it('keeps CommonJS bytes and require semantics stable across dependency load order', async () => {
  const root = await mkdtemp(path.join(tmpdir(), 'helixweave-commonjs-build-'));
  const dependency = path.join(root, 'node_modules', 'conditional-require');
  try {
    await mkdir(dependency, { recursive: true });
    await writeFile(path.join(dependency, 'value.cjs'), 'module.exports = { value: 42 };');
    await writeFile(
      path.join(dependency, 'conditional.cjs'),
      'exports.read = () => require("./value.cjs");',
    );
    await writeFile(
      path.join(dependency, 'unconditional.cjs'),
      'module.exports = require("./value.cjs");',
    );
    await writeFile(path.join(root, 'entry.js'), [
      'import conditional from "./node_modules/conditional-require/conditional.cjs";',
      'import unconditional from "./node_modules/conditional-require/unconditional.cjs";',
      'export const value = conditional.read().value;',
      'export const same = conditional.read() === unconditional;',
    ].join('\n'));

    const outputs: string[] = [];
    for (const first of ['conditional', 'unconditional']) {
      const preload: Plugin = {
        name: 'ordered-commonjs-fixture-load',
        async buildStart() {
          // Complete one dependency first; no sleeps or probabilistic race.
          await this.load({ id: path.join(dependency, `${first}.cjs`) });
        },
      };
      const result = await build({
        configFile: false,
        root,
        logLevel: 'silent',
        plugins: [preload],
        build: {
          ...config.build,
          rollupOptions: {
            ...config.build?.rollupOptions,
            input: undefined,
          },
          write: false,
          minify: false,
          lib: { entry: path.join(root, 'entry.js'), formats: ['cjs'], fileName: 'fixture' },
        },
      });
      const builds = Array.isArray(result) ? result : [result];
      expect(builds).toHaveLength(1);
      const output = builds[0];
      if (!('output' in output)) throw new Error('Expected a completed build');
      const chunk = output.output[0];
      if (chunk.type !== 'chunk') throw new Error('Expected JavaScript');
      outputs.push(chunk.code);
      const filename = path.join(root, `${first}.cjs`);
      await writeFile(filename, chunk.code);
      expect(createRequire(import.meta.url)(filename)).toMatchObject({ value: 42, same: true });
    }
    expect(outputs[1]).toBe(outputs[0]);
  } finally {
    await rm(root, { recursive: true, force: true });
  }
});

import { builtinEnvironments } from 'vitest/environments';

const jsdom = builtinEnvironments.jsdom;

interface NodeAbort {
  AbortController: typeof AbortController;
  AbortSignal: typeof AbortSignal;
}

// jsdom replaces AbortSignal after this setup starts. Node's Request only
// accepts the signal class it was built with, so keep that constructor.
export default {
  ...jsdom,
  name: 'jsdom',
  async setup(global: typeof globalThis, options: { jsdom?: object }) {
    const nodeAbort: NodeAbort = {
      AbortController: global.AbortController,
      AbortSignal: global.AbortSignal,
    };
    const result = await jsdom.setup(global, options);
    Object.defineProperty(global, '__nodeAbort', {
      value: nodeAbort,
      configurable: true,
      enumerable: false,
      writable: true,
    });
    return {
      teardown(next: typeof globalThis & { __nodeAbort?: NodeAbort }) {
        delete next.__nodeAbort;
        return result.teardown(next);
      },
    };
  },
};

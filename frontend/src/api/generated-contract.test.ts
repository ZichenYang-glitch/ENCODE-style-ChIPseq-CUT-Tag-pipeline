import { describe, it, expect } from 'vitest';
import { readFileSync } from 'node:fs';
import { fileURLToPath } from 'node:url';
import { resolve, dirname } from 'node:path';

const openapi = JSON.parse(
  readFileSync(
    resolve(dirname(fileURLToPath(import.meta.url)), '../../openapi.json'),
    'utf-8',
  ),
);

describe('generated OpenAPI client coverage', () => {
  const paths = Object.keys(openapi.paths);

  it('exports all public API paths in the OpenAPI contract', () => {
    expect(paths).toContain('/api/v1/workflows/');
    expect(paths).toContain('/api/v1/workflows/{workflow_id}/schema');
    expect(paths).toContain('/api/v1/workflows/{workflow_id}/validate');
    expect(paths).toContain('/api/v1/workflows/{workflow_id}/agent/chat');
    expect(paths).toContain('/api/v1/workflows/{workflow_id}/runs');
    expect(paths).toContain('/api/v1/runs/{run_id}');
    expect(paths).toContain('/api/v1/runs/{run_id}/cancel');
    expect(paths).toContain('/api/v1/runs/{run_id}/events');
    expect(paths).toContain('/api/v1/runs/{run_id}/logs');
    expect(paths).toContain('/api/v1/runs/{run_id}/preflight');
  });
});

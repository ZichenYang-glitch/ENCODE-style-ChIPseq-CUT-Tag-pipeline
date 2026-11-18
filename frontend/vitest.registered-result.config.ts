import { defineConfig } from 'vitest/config';

export default defineConfig({
  test: {
    environment: 'node',
    include: ['src/api/registeredResultDownload.acceptance.test.ts'],
  },
});

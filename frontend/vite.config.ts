/// <reference types="vitest" />
import { defineConfig } from 'vite';
import react from '@vitejs/plugin-react';
import path from 'node:path';

export default defineConfig({
  plugins: [react()],
  build: {
    // Auto-detecting conditional requires races with module loading (notably Ajv).
    // Preserve CommonJS initialization semantics independently of load order.
    commonjsOptions: { strictRequires: true },
    rollupOptions: {
      input: {
        main: path.resolve(__dirname, 'index.html'),
        viewer: path.resolve(__dirname, 'viewer.html'),
      },
    },
  },
  server: {
    proxy: {
      '/api': process.env.VITE_API_PROXY_TARGET ?? 'http://127.0.0.1:8000',
    },
  },
  resolve: {
    alias: {
      '@': path.resolve(__dirname, './src'),
    },
  },
  test: {
    globals: true,
    environment: './src/test/jsdom-environment.ts',
    setupFiles: ['./src/test/setup.ts'],
    exclude: ['e2e/**', 'e2e-synthetic/**', 'node_modules/**', 'dist/**'],
  },
});

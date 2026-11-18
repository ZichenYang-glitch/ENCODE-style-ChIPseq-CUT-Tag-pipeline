/** @type {import('orval').Options} */
module.exports = {
  workflow: {
    input: {
      target: './openapi.json',
    },
    output: {
      client: 'fetch',
      target: './src/api/generated/',
      schemas: './src/api/generated/models',
      mode: 'tags-split',
      clean: true,
      // Catalog writes declare a required Idempotency-Key header. Orval
      // only emits header parameters when this flag is on.
      headers: true,
      override: {
        mutator: {
          path: './src/api/fetcher.ts',
          name: 'fetcher',
        },
        fetch: {
          includeHttpResponseReturnType: false,
        },
        operations: {
          downloadRunArtifact: {
            mutator: {
              path: './src/api/fetcher.ts',
              name: 'blobFetcher',
            },
          },
          readRunBrowserTrackBytes: {
            mutator: {
              path: './src/api/fetcher.ts',
              name: 'blobFetcher',
            },
          },
          headRunBrowserTrackBytes: {
            mutator: {
              path: './src/api/fetcher.ts',
              name: 'blobFetcher',
            },
          },
        },
      },
    },
  },
};

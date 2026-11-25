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
          // OpenAPI form arrays repeat the query key. Orval otherwise joins
          // the values with a comma, which catalog filters reject.
          explode: true,
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
          catalogReadArchivedTrackContent: {
            mutator: {
              path: './src/api/fetcher.ts',
              name: 'blobFetcher',
            },
          },
          catalogHeadArchivedTrackContent: {
            mutator: {
              path: './src/api/fetcher.ts',
              name: 'blobFetcher',
            },
          },
          catalogReadArchivedTrackIndex: {
            mutator: {
              path: './src/api/fetcher.ts',
              name: 'blobFetcher',
            },
          },
          catalogHeadArchivedTrackIndex: {
            mutator: {
              path: './src/api/fetcher.ts',
              name: 'blobFetcher',
            },
          },
          catalogDownloadRegisteredResultRevision: {
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

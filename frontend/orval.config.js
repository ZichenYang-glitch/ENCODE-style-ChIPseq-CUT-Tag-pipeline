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
      override: {
        mutator: {
          path: './src/api/fetcher.ts',
          name: 'fetcher',
        },
        fetch: {
          includeHttpResponseReturnType: false,
        },
      },
    },
  },
};

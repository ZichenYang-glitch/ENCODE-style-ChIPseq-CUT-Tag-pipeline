export default (context) => {
  const file = typeof context.file === 'string' ? context.file : '';
  if (file.includes('/node_modules/')) {
    return { plugins: {} };
  }
  return {
    plugins: {
      tailwindcss: {},
      autoprefixer: {},
    },
  };
};

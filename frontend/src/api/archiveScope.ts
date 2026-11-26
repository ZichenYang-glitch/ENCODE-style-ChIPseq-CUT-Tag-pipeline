export const PROJECT_ID = /^prj_[0-9a-f]{32}$/;
export const EXPERIMENT_ID = /^exp_[0-9a-f]{32}$/;
export const ARCHIVE_VERSION_ID = /^aver_[0-9a-f]{32}$/;

export function sameOriginApiUrl(path: string): string | null {
  if (!path.startsWith('/api/')) return null;
  const base = import.meta.env.VITE_API_BASE_URL?.trim().replace(/\/$/, '') ?? '';
  try {
    const url = new URL(base + path, window.location.origin);
    return url.origin === window.location.origin ? url.href : null;
  } catch {
    return null;
  }
}

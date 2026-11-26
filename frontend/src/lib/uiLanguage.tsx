import { useSyncExternalStore } from 'react';

export type UiLanguage = 'en' | 'zh-CN';

const preferenceKey = 'helixweave.ui-language';
const changeEvent = 'helixweave-ui-language-change';
let inMemoryChoice: UiLanguage | null = null;

export function browserLanguage(preferences: readonly string[]): UiLanguage {
  for (const tag of preferences) {
    const primary = tag.toLowerCase().split('-')[0];
    if (primary === 'en') return 'en';
    if (primary === 'zh') return 'zh-CN';
  }
  return 'en';
}

function currentLanguage(): UiLanguage {
  if (inMemoryChoice !== null) return inMemoryChoice;
  try {
    const saved = window.localStorage.getItem(preferenceKey);
    if (saved === 'zh-CN' || saved === 'en') return saved;
  } catch {
    // The browser can disable site storage while the page is still usable.
  }
  const preferences = navigator.languages;
  return browserLanguage(preferences?.length ? preferences : [navigator.language]);
}

export function setUiLanguage(choice: UiLanguage) {
  try {
    window.localStorage.setItem(preferenceKey, choice);
    inMemoryChoice = null;
  } catch {
    inMemoryChoice = choice;
  }
  window.dispatchEvent(new Event(changeEvent));
}

function subscribe(notify: () => void) {
  function onStorage(event: StorageEvent) {
    if (event.key === preferenceKey || event.key === null) {
      inMemoryChoice = null;
      notify();
    }
  }
  window.addEventListener('storage', onStorage);
  window.addEventListener('languagechange', notify);
  window.addEventListener(changeEvent, notify);
  return () => {
    window.removeEventListener('storage', onStorage);
    window.removeEventListener('languagechange', notify);
    window.removeEventListener(changeEvent, notify);
  };
}

export function useUiLanguage(): UiLanguage {
  return useSyncExternalStore(subscribe, currentLanguage, () => 'en');
}

export function UiLanguageSelect() {
  const language = useUiLanguage();
  return (
    <label className="inline-flex items-center gap-2 text-sm" lang={language}>
      {language === 'en' ? 'Language' : '语言'}
      <select
        className="rounded border border-[var(--color-border)] bg-[var(--color-surface)] px-2 py-1"
        value={language}
        onChange={(event) => setUiLanguage(event.target.value === 'zh-CN' ? 'zh-CN' : 'en')}
      >
        <option value="en" lang="en">English</option>
        <option value="zh-CN" lang="zh-CN">简体中文</option>
      </select>
    </label>
  );
}

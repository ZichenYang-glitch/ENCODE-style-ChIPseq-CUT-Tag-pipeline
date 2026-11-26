import { act, render, screen } from '@testing-library/react';
import userEvent from '@testing-library/user-event';
import { afterEach, expect, it, vi } from 'vitest';
import { browserLanguage, UiLanguageSelect } from './uiLanguage';

afterEach(() => {
  localStorage.removeItem('helixweave.ui-language');
  vi.restoreAllMocks();
});

it('picks the first supported browser preference and falls back to English', () => {
  expect(browserLanguage(['fr-CA', 'en-GB', 'zh-CN'])).toBe('en');
  expect(browserLanguage(['fr', 'zh-Hans-CN', 'en-US'])).toBe('zh-CN');
  expect(browserLanguage(['de-DE'])).toBe('en');
  expect(browserLanguage([])).toBe('en');
});

it('remembers an explicit choice and responds to changes from another tab', async () => {
  vi.spyOn(navigator, 'languages', 'get').mockReturnValue(['zh-CN']);
  const first = render(<UiLanguageSelect />);
  await userEvent.setup().selectOptions(screen.getByRole('combobox', { name: '语言' }), 'en');
  first.unmount();
  render(<UiLanguageSelect />);
  expect(screen.getByRole('combobox', { name: 'Language' })).toHaveValue('en');
  act(() => {
    localStorage.removeItem('helixweave.ui-language');
    window.dispatchEvent(new StorageEvent('storage', { key: 'helixweave.ui-language' }));
  });
  expect(screen.getByRole('combobox', { name: '语言' })).toHaveValue('zh-CN');
});

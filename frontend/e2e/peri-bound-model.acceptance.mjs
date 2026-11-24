import { readFileSync } from 'node:fs';
import { execFileSync } from 'node:child_process';
import { chromium, expect } from '@playwright/test';

const fixture = JSON.parse(readFileSync(0, 'utf8'));
const browser = await chromium.launch();
try {
  const ctx = await browser.newContext({baseURL: fixture.url});
  await ctx.addCookies([fixture.session,fixture.csrf].map(([name,value])=>({name,value,url:fixture.url})));
  const page = await ctx.newPage();
  await page.goto(`/workflows/${fixture.workflow}/configuration-drafts`);
  await page.getByRole('combobox',{name:'Project'}).selectOption(fixture.projectId);
  await page.getByRole('button',{name:/Peri run action fixture/}).click();
  await page.getByRole('button',{name:'Connect assistant'}).click();
  await expect(page.getByText('Assistant available')).toBeVisible();
  // Run as the local test operator, never in the browser or in a web route.
  execFileSync(fixture.grant[0],fixture.grant.slice(1));
  await page.getByRole('radio',{name:'Suggest values'}).check();
  await page.getByRole('checkbox',{name:'Threads'}).check();
  await page.getByLabel('Preferred value').fill('3');
  const wait = page.waitForResponse(r=>r.request().method()==='POST' && r.url().endsWith('/turns'));
  await page.getByRole('button',{name:'Request suggestion'}).click();
  const r=await wait, b=await r.json();
  expect(r.status()).toBe(200);
  expect(b.result_kind).toBe('proposal');
  expect(b.adoptable).toBe(true);
  await expect(page.getByRole('table',{name:'Suggested changes'})).toContainText('7');
  const next = page.waitForResponse(response =>
    response.request().method() === 'POST' && response.url() === r.url());
  await page.getByRole('button', {name: 'Request suggestion'}).click();
  const denied = await next;
  const second = await denied.json();
  expect(denied.status()).toBe(200);
  expect(second.status).toBe('failed');
  expect(second.failure_reason).toBe('SCOPE_DENIED');
  expect(second.adoptable).toBe(false);
  const previous = await ctx.request.get(r.url() + '/' + b.turn_id);
  expect(previous.status()).toBe(200);
  const pending = await previous.json();
  console.log(JSON.stringify({second: second.failure_reason, first: pending.proposal_state, adoptable: pending.adoptable}));
  expect(pending.proposal_state).toBe('pending');
  expect(pending.adoptable).toBe(true);
  await expect(page.getByRole('alert').filter({hasText: 'The suggestion request did not succeed'})).toBeVisible();
  await expect(page.getByRole('table', {name: 'Suggested changes'})).toContainText('7');
  const receipt = page.waitForResponse(response => response.url().endsWith('/adopt'));
  await page.getByRole('button', {name: 'Adopt suggestion'}).click();
  expect((await receipt).status()).toBe(200);
  console.log(JSON.stringify({api:new URL(r.url()).pathname,http:r.status(),requested:3,suggested:7,adopted:true,second:second.status}));
} finally { await browser.close(); }

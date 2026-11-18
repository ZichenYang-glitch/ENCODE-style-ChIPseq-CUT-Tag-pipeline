import { describe, expect, it } from 'vitest';
import { washuAcceptFailed } from './washu-accept-verdict.mjs';

const positive = {
  code: 0,
  output: 'WASHU_PROOF ok {"window":"chr1:500-2700"}\n',
};

const negative = {
  code: 0,
  output: 'WASHU_PROOF rejected protected track read\nproof work /tmp/negative\n',
};

describe('WashU accept verdict', () => {
  it('passes a drawn positive proof and a confirmed protected-track rejection', () => {
    expect(washuAcceptFailed(positive, negative)).toBe(false);
  });

  it('fails when the negative proof exits non-zero, even if it mentions rejection', () => {
    expect(
      washuAcceptFailed(positive, {
        code: 1,
        output: 'WASHU_PROOF rejected protected track read\nunexpected browser crash\n',
      }),
    ).toBe(true);
  });

  it('fails when the negative proof exits 0 without the rejection result', () => {
    expect(washuAcceptFailed(positive, { code: 0, output: 'browser timed out\n' })).toBe(true);
    expect(washuAcceptFailed(positive, { code: 0, tail: 'proof work only\n' })).toBe(true);
  });

  it('fails when a negative result also claims the positive proof passed', () => {
    expect(
      washuAcceptFailed(positive, {
        code: 0,
        output: 'WASHU_PROOF ok {"window":"chr1"}\nWASHU_PROOF rejected protected track read\n',
      }),
    ).toBe(true);
  });

  it('fails when the positive proof does not both exit 0 and report ok', () => {
    expect(washuAcceptFailed({ code: 1, output: 'WASHU_PROOF ok {}\n' }, negative)).toBe(true);
    expect(washuAcceptFailed({ code: 0, output: 'viewer never drew\n' }, negative)).toBe(true);
    expect(washuAcceptFailed({ code: 0, output: '' }, negative)).toBe(true);
  });
});

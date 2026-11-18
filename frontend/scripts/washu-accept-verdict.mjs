/**
 * Decide whether the WashU accept summary failed.
 *
 * A negative proof passes only when the protected-track rejection was the
 * expected result: exit 0, text contains `WASHU_PROOF rejected`, and text
 * does not contain `WASHU_PROOF ok`. Any other negative outcome is an
 * unexpected proof error and fails the summary. The positive proof passes
 * only with exit 0 and `WASHU_PROOF ok`.
 *
 * @param {{ code?: number | null, output?: string, tail?: string }} positive
 * @param {{ code?: number | null, output?: string, tail?: string }} negative
 */
export function washuAcceptFailed(positive, negative) {
  const positiveText = proofText(positive);
  const negativeText = proofText(negative);
  const positivePassed = positive?.code === 0 && positiveText.includes('WASHU_PROOF ok');
  const negativePassed =
    negative?.code === 0 &&
    negativeText.includes('WASHU_PROOF rejected') &&
    !negativeText.includes('WASHU_PROOF ok');
  return !positivePassed || !negativePassed;
}

function proofText(result) {
  if (typeof result?.output === 'string') return result.output;
  if (typeof result?.tail === 'string') return result.tail;
  return '';
}

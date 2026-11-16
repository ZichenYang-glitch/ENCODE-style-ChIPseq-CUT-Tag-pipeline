import { Panel } from '../../components/Panel';
import type { ExplanationView } from './model';

interface ExplanationNoticeProps {
  explanation: ExplanationView | null;
}

export function ExplanationNotice({ explanation }: ExplanationNoticeProps) {
  return (
    <Panel title="Explanation">
      {!explanation || !explanation.contentAvailable ? (
        <p role="status">
          The explanation was not included in the response. An explanation
          cannot be adopted.
        </p>
      ) : explanation.labels.length === 0 ? (
        <p role="status">
          The explanation did not include any approved reasons. An explanation
          cannot be adopted.
        </p>
      ) : (
        <div className="space-y-2">
          <p>
            This explains a validation issue. It does not change the draft.
          </p>
          <ul aria-label="Explanation reasons" className="list-disc space-y-1 pl-5 text-sm">
            {explanation.labels.map((label) => (
              <li key={label}>{label}</li>
            ))}
          </ul>
        </div>
      )}
    </Panel>
  );
}

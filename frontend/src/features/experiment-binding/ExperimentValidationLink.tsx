import { Link } from 'react-router-dom';
import { Button } from '../../components/Button';
import { Panel } from '../../components/Panel';
import {
  EXPERIMENT_VALIDATION_HANDOFF_KEY,
  type DraftCarryResult,
} from './handoff';

interface ExperimentValidationLinkProps {
  workflowId: string;
  result: DraftCarryResult;
  disabled: boolean;
  disabledReason: string | null;
}

export function ExperimentValidationLink({
  workflowId,
  result,
  disabled,
  disabledReason,
}: ExperimentValidationLinkProps) {
  const handoff = result.handoff;
  const blocked = disabled || handoff === null;
  return (
    <Panel title="Experiment validation">
      <div className="space-y-3">
        <p>
          This opens the ordinary validation page with the saved draft&apos;s
          configuration, project, and sample revisions. It does not validate
          the draft, does not create a run, and does not bind an experiment.
        </p>
        <ul className="list-disc space-y-1 pl-5 text-sm">
          {result.notes.map((note) => (
            <li key={note}>{note}</li>
          ))}
        </ul>
        {blocked ? (
          <p className="text-sm text-[var(--color-text-muted)]">
            {disabledReason ??
              'This saved draft cannot open experiment validation yet.'}
          </p>
        ) : (
          <Button asChild variant="secondary">
            <Link
              to={`/workflows/${encodeURIComponent(workflowId)}/new-run?step=review`}
              state={{ [EXPERIMENT_VALIDATION_HANDOFF_KEY]: handoff }}
            >
              Open experiment validation
            </Link>
          </Button>
        )}
      </div>
    </Panel>
  );
}

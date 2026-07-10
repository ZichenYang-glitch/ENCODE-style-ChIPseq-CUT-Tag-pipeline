import { Link } from 'react-router-dom';
import { Button } from '@/components/ui/button';

export function NotFoundPage() {
  return (
    <div className="flex flex-1 flex-col items-center justify-center gap-4 p-6 text-center">
      <h2 className="text-2xl font-semibold text-[var(--color-text)]">404</h2>
      <p className="text-sm text-[var(--color-text-muted)]">
        This page does not exist.
      </p>
      <Link to="/workflows">
        <Button variant="primary">Back to workflows</Button>
      </Link>
    </div>
  );
}

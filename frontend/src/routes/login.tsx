import { useState, type FormEvent } from 'react';
import { Navigate, useLocation, useNavigate } from 'react-router-dom';
import { Dna } from 'lucide-react';
import { useAuth } from '../app/auth';
import { ApiError } from '../api/fetcher';
import { Button } from '../components/Button';

interface LocationState {
  from?: string;
}

export function LoginPage() {
  const { loading, setupRequired, principal, login } = useAuth();
  const navigate = useNavigate();
  const location = useLocation();
  const from = (location.state as LocationState | null)?.from ?? '/workflows';
  const [username, setUsername] = useState('');
  const [password, setPassword] = useState('');
  const [error, setError] = useState<string | null>(null);
  const [submitting, setSubmitting] = useState(false);

  if (!loading && principal !== null) {
    return <Navigate to={from} replace />;
  }

  async function handleSubmit(event: FormEvent<HTMLFormElement>) {
    event.preventDefault();
    setSubmitting(true);
    setError(null);
    try {
      await login(username, password);
      navigate(from, { replace: true });
    } catch (caught) {
      if (caught instanceof ApiError) {
        setError(caught.message);
      } else {
        setError('Login failed.');
      }
    } finally {
      setSubmitting(false);
    }
  }

  return (
    <div className="flex min-h-screen items-center justify-center bg-[var(--color-bg)] px-4">
      <section
        aria-labelledby="login-heading"
        className="w-full max-w-sm rounded-lg border border-[var(--color-border)] bg-[var(--color-surface)] p-6 shadow-sm"
      >
        <h1
          id="login-heading"
          className="flex items-center gap-2 text-lg font-semibold text-[var(--color-accent)]"
        >
          <Dna aria-hidden="true" size={20} />
          HelixWeave
        </h1>
        {setupRequired ? (
          <p className="mt-4 text-sm text-[var(--color-text-muted)]" role="status">
            No administrator exists yet. Create the first one from the operator
            console with{' '}
            <code className="rounded bg-slate-100 px-1 py-0.5 text-xs">
              helixweave admin account bootstrap
            </code>
            , then sign in here.
          </p>
        ) : (
          <form className="mt-6 space-y-4" onSubmit={handleSubmit}>
            <div>
              <label
                htmlFor="login-username"
                className="block text-sm font-medium text-[var(--color-text)]"
              >
                Username
              </label>
              <input
                id="login-username"
                name="username"
                type="text"
                autoComplete="username"
                required
                minLength={3}
                maxLength={64}
                value={username}
                onChange={(event) => setUsername(event.target.value)}
                className="mt-1 block w-full rounded-md border border-[var(--color-border)] bg-white px-3 py-2 text-sm focus:border-[var(--color-accent)] focus:outline-none focus:ring-1 focus:ring-[var(--color-accent)]"
              />
            </div>
            <div>
              <label
                htmlFor="login-password"
                className="block text-sm font-medium text-[var(--color-text)]"
              >
                Password
              </label>
              <input
                id="login-password"
                name="password"
                type="password"
                autoComplete="current-password"
                required
                value={password}
                onChange={(event) => setPassword(event.target.value)}
                className="mt-1 block w-full rounded-md border border-[var(--color-border)] bg-white px-3 py-2 text-sm focus:border-[var(--color-accent)] focus:outline-none focus:ring-1 focus:ring-[var(--color-accent)]"
              />
            </div>
            {error !== null && (
              <p role="alert" className="text-sm text-red-600">
                {error}
              </p>
            )}
            <Button type="submit" disabled={submitting} className="w-full">
              {submitting ? 'Signing in…' : 'Sign in'}
            </Button>
          </form>
        )}
      </section>
    </div>
  );
}

import { lazy, Suspense } from 'react';
import {
  createBrowserRouter,
  Navigate,
  useParams,
  type RouteObject,
} from 'react-router-dom';
import { AppShell } from './AppShell';
import { NotFoundPage } from '../routes/not-found';
import { WorkflowsLayout } from '../routes/workflows/layout';
import { WorkflowsIndexPage } from '../routes/workflows/index';
import { WorkflowDetailPage } from '../routes/workflows/detail';
import { RunDetailPage } from '../routes/runs/detail';

const NewRunWorkbenchPage = lazy(() =>
  import('../routes/workflows/new-run').then((module) => ({
    default: module.NewRunWorkbenchPage,
  })),
);

const RunHistoryPage = lazy(() =>
  import('../routes/runs').then((module) => ({
    default: module.RunHistoryPage,
  })),
);

function WorkflowDetailRoute() {
  const { workflowId } = useParams<{ workflowId: string }>();
  if (!workflowId) {
    return <NotFoundPage />;
  }
  return <WorkflowDetailPage key={workflowId} workflowId={workflowId} />;
}

function NewRunWorkbenchRoute() {
  const { workflowId } = useParams<{ workflowId: string }>();
  if (!workflowId) {
    return <NotFoundPage />;
  }
  return (
    <Suspense fallback={<WorkbenchRouteLoading />}>
      <NewRunWorkbenchPage key={workflowId} workflowId={workflowId} />
    </Suspense>
  );
}

function WorkbenchRouteLoading() {
  return (
    <section
      data-testid="input-workbench-loading"
      className="min-h-[38rem] min-w-0 flex-1 animate-pulse rounded-lg border border-[var(--color-border)] bg-[var(--color-surface)] p-4"
      aria-label="Loading input workbench"
    >
      <div className="h-5 w-48 rounded bg-slate-200" />
      <div className="mt-4 h-10 rounded bg-slate-100" />
      <div className="mt-4 h-96 rounded bg-slate-100" />
    </section>
  );
}

export const appRoutes: RouteObject[] = [
  {
    path: '/',
    element: <AppShell />,
    children: [
      { index: true, element: <Navigate to="/workflows" replace /> },
      {
        path: 'workflows',
        element: <WorkflowsLayout />,
        children: [
          { index: true, element: <WorkflowsIndexPage /> },
          { path: ':workflowId/new-run', element: <NewRunWorkbenchRoute /> },
          { path: ':workflowId', element: <WorkflowDetailRoute /> },
        ],
      },
      {
        path: 'runs',
        element: (
          <Suspense fallback={<RunHistoryRouteLoading />}>
            <RunHistoryPage />
          </Suspense>
        ),
      },
      { path: 'runs/:runId', element: <RunDetailPage /> },
      { path: '*', element: <NotFoundPage /> },
    ],
  },
];

export function createAppRouter() {
  return createBrowserRouter(appRoutes);
}

function RunHistoryRouteLoading() {
  return (
    <section
      className="min-h-[30rem] min-w-0 flex-1 animate-pulse rounded-lg border border-[var(--color-border)] bg-[var(--color-surface)] p-4"
      aria-label="Loading run history page"
    >
      <div className="h-5 w-40 rounded bg-slate-200" />
      <div className="mt-4 h-10 rounded bg-slate-100" />
      <div className="mt-4 h-64 rounded bg-slate-100" />
    </section>
  );
}

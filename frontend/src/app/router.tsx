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

function WorkflowDetailRoute() {
  const { workflowId } = useParams<{ workflowId: string }>();
  if (!workflowId) {
    return <NotFoundPage />;
  }
  return <WorkflowDetailPage key={workflowId} workflowId={workflowId} />;
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
          { path: ':workflowId', element: <WorkflowDetailRoute /> },
        ],
      },
      { path: 'runs/:runId', element: <RunDetailPage /> },
      { path: '*', element: <NotFoundPage /> },
    ],
  },
];

export function createAppRouter() {
  return createBrowserRouter(appRoutes);
}

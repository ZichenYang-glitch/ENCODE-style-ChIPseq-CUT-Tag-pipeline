import {
  createContext,
  useCallback,
  useContext,
  useMemo,
  type ReactNode,
} from 'react';
import { useQuery, useQueryClient } from '@tanstack/react-query';
import {
  fetchSessionState,
  login as loginRequest,
  logout as logoutRequest,
  type AuthPrincipal,
} from '../api/authClient';

export interface AuthContextValue {
  loading: boolean;
  setupRequired: boolean;
  principal: AuthPrincipal | null;
  login: (username: string, password: string) => Promise<void>;
  logout: () => Promise<void>;
  refresh: () => Promise<void>;
}

const AuthContext = createContext<AuthContextValue | null>(null);

export function AuthProvider({ children }: { children: ReactNode }) {
  const queryClient = useQueryClient();
  const session = useQuery({
    queryKey: ['auth', 'session'],
    queryFn: fetchSessionState,
    retry: false,
    staleTime: 30_000,
  });

  const refresh = useCallback(async () => {
    await queryClient.invalidateQueries({ queryKey: ['auth', 'session'] });
  }, [queryClient]);

  const login = useCallback(
    async (username: string, password: string) => {
      await loginRequest(username, password);
      await queryClient.invalidateQueries({ queryKey: ['auth', 'session'] });
    },
    [queryClient],
  );

  const logout = useCallback(async () => {
    await logoutRequest();
    queryClient.clear();
    await queryClient.invalidateQueries({ queryKey: ['auth', 'session'] });
  }, [queryClient]);

  const value = useMemo<AuthContextValue>(
    () => ({
      loading: session.isPending,
      setupRequired: session.data?.setup_required ?? false,
      principal: session.data?.authenticated ? session.data.principal : null,
      login,
      logout,
      refresh,
    }),
    [session.isPending, session.data, login, logout, refresh],
  );

  return <AuthContext.Provider value={value}>{children}</AuthContext.Provider>;
}

export function useAuth(): AuthContextValue {
  const value = useContext(AuthContext);
  if (value === null) {
    throw new Error('useAuth must be used inside AuthProvider');
  }
  return value;
}

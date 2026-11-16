declare module 'wuepgg' {
  import type { ComponentType } from 'react';

  export interface WashuTrackConfig {
    type: string;
    name?: string;
    url?: string;
    options?: Record<string, unknown>;
  }

  export const GenomeHub: ComponentType<{
    storeConfig: { storeId: string; enablePersistence?: boolean };
    genomeName: string;
    tracks: WashuTrackConfig[];
    viewRegion?: string;
    showGenomeNavigator?: boolean;
    showNavBar?: boolean;
    showToolBar?: boolean;
    showDisclosure?: boolean;
  }>;
}

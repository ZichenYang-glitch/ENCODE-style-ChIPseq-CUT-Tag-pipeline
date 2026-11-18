/**
 * Version 1 messages between the platform host and the same-origin viewer page.
 *
 * The instance id identifies one iframe document. It is not a login session,
 * cookie, or credential. Track revisions and byte URLs stay on the authorized
 * browser API; this payload only names the run, optional artifact ids, the
 * initial region, and display options.
 *
 * The framed page posts `listening` only after its message handler is
 * installed, and repeats that announcement until configure or close arrives.
 * One document load uses one epoch, so a repeat is not a new document. The
 * host sends configure in response. `accepted` confirms that this document
 * applied the configure. A retired instance id is ignored, and a destroyed
 * frame is closed instead of left waiting on the old exchange.
 */

export const VIEWER_PROTOCOL_VERSION = 1;
export const VIEWER_MESSAGE_LIMIT = 4096;
export const VIEWER_TRACK_LIMIT = 32;

const INSTANCE_ID =
  /^[0-9a-f]{8}-[0-9a-f]{4}-4[0-9a-f]{3}-[89ab][0-9a-f]{3}-[0-9a-f]{12}$/i;
const RUN_ID = /^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$/;
const ARTIFACT_ID = /^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$/;
const REGION = /^[A-Za-z0-9_.:-]{1,128}:[0-9]{1,10}-[0-9]{1,10}$/;
const STATUS = /^(?:bound|unbound|unavailable)$/;
const HOST_TYPES = new Set([
  'helixweave.viewer.configure',
  'helixweave.viewer.close',
]);
const VIEWER_TYPES = new Set([
  'helixweave.viewer.ready',
  'helixweave.viewer.failed',
  'helixweave.viewer.listening',
  'helixweave.viewer.accepted',
]);

export type ViewerReferenceStatus = 'bound' | 'unbound' | 'unavailable';

export interface ViewerConfigurePayload {
  runId: string;
  region: string | null;
  artifactIds: string[] | null;
  display: { showToolbar: boolean };
}

export interface ViewerReadyPayload {
  runId: string;
  referenceStatus: ViewerReferenceStatus;
  artifactIds: string[];
}

export interface ViewerFailedPayload {
  runId: string;
  message: string;
}

export type ViewerHostCommand =
  | {
      version: 1;
      type: 'helixweave.viewer.configure';
      instanceId: string;
      payload: ViewerConfigurePayload;
    }
  | {
      version: 1;
      type: 'helixweave.viewer.close';
      instanceId: string;
      payload: Record<string, never>;
    };

export type ViewerPageNotice =
  | {
      version: 1;
      type: 'helixweave.viewer.ready';
      instanceId: string;
      payload: ViewerReadyPayload;
    }
  | {
      version: 1;
      type: 'helixweave.viewer.failed';
      instanceId: string;
      payload: ViewerFailedPayload;
    }
  | {
      version: 1;
      type: 'helixweave.viewer.listening';
      instanceId: string;
      payload: { epoch: string };
    }
  | {
      version: 1;
      type: 'helixweave.viewer.accepted';
      instanceId: string;
      payload: { runId: string };
    };

export function newViewerInstanceId(): string {
  return crypto.randomUUID();
}

export function configureViewerMessage(
  instanceId: string,
  payload: ViewerConfigurePayload,
): ViewerHostCommand {
  return {
    version: VIEWER_PROTOCOL_VERSION,
    type: 'helixweave.viewer.configure',
    instanceId,
    payload,
  };
}

export function closeViewerMessage(instanceId: string): ViewerHostCommand {
  return {
    version: VIEWER_PROTOCOL_VERSION,
    type: 'helixweave.viewer.close',
    instanceId,
    payload: {},
  };
}

export function listeningViewerMessage(instanceId: string, epoch: string): ViewerPageNotice {
  return {
    version: VIEWER_PROTOCOL_VERSION,
    type: 'helixweave.viewer.listening',
    instanceId,
    payload: { epoch },
  };
}

export function acceptedViewerMessage(instanceId: string, runId: string): ViewerPageNotice {
  return {
    version: VIEWER_PROTOCOL_VERSION,
    type: 'helixweave.viewer.accepted',
    instanceId,
    payload: { runId },
  };
}

export function postViewerMessage(
  target: Window,
  message: ViewerHostCommand | ViewerPageNotice,
  origin: string,
): void {
  if (origin === '*' || origin === '' || origin === 'null') return;
  target.postMessage(message, origin);
}

export function readHostCommand(
  data: unknown,
  instanceId: string,
  event: { origin: string; source: unknown },
  expected: { origin: string; source: object | null },
): ViewerHostCommand | null {
  if (!sameChannel(event, expected) || !sameInstance(data, instanceId)) return null;
  return parseHostCommand(data);
}

export function readViewerNotice(
  data: unknown,
  instanceId: string,
  event: { origin: string; source: unknown },
  expected: { origin: string; source: object | null },
): ViewerPageNotice | null {
  if (!sameChannel(event, expected) || !sameInstance(data, instanceId)) return null;
  return parseViewerNotice(data);
}

export function standaloneRunId(search: string): string | null {
  const params = new URLSearchParams(search);
  const runId = params.get('run');
  if (runId === null || !RUN_ID.test(runId)) return null;
  return runId;
}

export function framedInstanceId(hash: string): string | null {
  const params = new URLSearchParams(hash.startsWith('#') ? hash.slice(1) : hash);
  const instanceId = params.get('instance');
  if (instanceId === null || !INSTANCE_ID.test(instanceId)) return null;
  return instanceId;
}

function sameChannel(
  event: { origin: string; source: unknown },
  expected: { origin: string; source: object | null },
): boolean {
  return (
    expected.source !== null &&
    event.origin === expected.origin &&
    expected.origin !== '*' &&
    expected.origin !== '' &&
    event.source === expected.source
  );
}

function sameInstance(data: unknown, instanceId: string): boolean {
  if (!isRecord(data)) return false;
  return data.instanceId === instanceId && INSTANCE_ID.test(instanceId);
}

function parseHostCommand(data: unknown): ViewerHostCommand | null {
  if (!envelope(data, HOST_TYPES)) return null;
  const record = data as Record<string, unknown>;
  if (record.type === 'helixweave.viewer.close') {
    if (!emptyRecord(record.payload)) return null;
    return {
      version: 1,
      type: 'helixweave.viewer.close',
      instanceId: record.instanceId as string,
      payload: {},
    };
  }
  const payload = configurePayload(record.payload);
  if (payload === null) return null;
  return {
    version: 1,
    type: 'helixweave.viewer.configure',
    instanceId: record.instanceId as string,
    payload,
  };
}

function parseViewerNotice(data: unknown): ViewerPageNotice | null {
  if (!envelope(data, VIEWER_TYPES)) return null;
  const record = data as Record<string, unknown>;
  if (record.type === 'helixweave.viewer.listening') {
    if (!isRecord(record.payload) || Object.keys(record.payload).length !== 1) return null;
    const epoch = text(record.payload.epoch);
    if (epoch === null || !INSTANCE_ID.test(epoch)) return null;
    return {
      version: 1,
      type: 'helixweave.viewer.listening',
      instanceId: record.instanceId as string,
      payload: { epoch },
    };
  }
  const payload = isRecord(record.payload) ? record.payload : null;
  if (payload === null || !RUN_ID.test(text(payload.runId) ?? '')) return null;
  if (record.type === 'helixweave.viewer.accepted') {
    if (Object.keys(payload).length !== 1) return null;
    return {
      version: 1,
      type: 'helixweave.viewer.accepted',
      instanceId: record.instanceId as string,
      payload: { runId: payload.runId as string },
    };
  }
  if (record.type === 'helixweave.viewer.failed') {
    const message = text(payload.message);
    if (message === null || message.length < 1 || message.length > 300) return null;
    if (Object.keys(payload).length !== 2) return null;
    return {
      version: 1,
      type: 'helixweave.viewer.failed',
      instanceId: record.instanceId as string,
      payload: { runId: payload.runId as string, message },
    };
  }
  const status = text(payload.referenceStatus);
  const artifactIds = idList(payload.artifactIds);
  if (status === null || !STATUS.test(status) || artifactIds === null) return null;
  if (Object.keys(payload).length !== 3) return null;
  return {
    version: 1,
    type: 'helixweave.viewer.ready',
    instanceId: record.instanceId as string,
    payload: {
      runId: payload.runId as string,
      referenceStatus: status as ViewerReferenceStatus,
      artifactIds,
    },
  };
}

function envelope(data: unknown, types: Set<string>): boolean {
  if (!isRecord(data) || Object.keys(data).length !== 4) return false;
  if (data.version !== VIEWER_PROTOCOL_VERSION) return false;
  if (typeof data.type !== 'string' || !types.has(data.type)) return false;
  if (typeof data.instanceId !== 'string' || !INSTANCE_ID.test(data.instanceId)) return false;
  const encoded = encodedLength(data);
  return encoded !== null && encoded <= VIEWER_MESSAGE_LIMIT;
}

function configurePayload(value: unknown): ViewerConfigurePayload | null {
  if (!isRecord(value) || Object.keys(value).length !== 4) return null;
  const runId = text(value.runId);
  if (runId === null || !RUN_ID.test(runId)) return null;
  if (value.region !== null && (text(value.region) === null || !REGION.test(String(value.region)))) {
    return null;
  }
  let artifactIds: string[] | null;
  if (value.artifactIds === null) {
    artifactIds = null;
  } else {
    const parsed = idList(value.artifactIds);
    if (parsed === null) return null;
    artifactIds = parsed;
  }
  if (!isRecord(value.display) || Object.keys(value.display).length !== 1) return null;
  if (typeof value.display.showToolbar !== 'boolean') return null;
  return {
    runId,
    region: value.region === null ? null : String(value.region),
    artifactIds: artifactIds ?? null,
    display: { showToolbar: value.display.showToolbar },
  };
}

function idList(value: unknown): string[] | null {
  if (!Array.isArray(value) || value.length > VIEWER_TRACK_LIMIT) return null;
  const ids: string[] = [];
  for (const item of value) {
    if (typeof item !== 'string' || !ARTIFACT_ID.test(item) || ids.includes(item)) return null;
    ids.push(item);
  }
  return ids;
}

function encodedLength(value: unknown): number | null {
  try {
    const encoded = JSON.stringify(value);
    if (typeof encoded !== 'string') return null;
    return encoded.length;
  } catch {
    return null;
  }
}

function emptyRecord(value: unknown): boolean {
  return isRecord(value) && Object.keys(value).length === 0;
}

function text(value: unknown): string | null {
  if (typeof value !== 'string' || value.trim() !== value) return null;
  return value;
}

function isRecord(value: unknown): value is Record<string, unknown> {
  return value !== null && typeof value === 'object' && !Array.isArray(value);
}

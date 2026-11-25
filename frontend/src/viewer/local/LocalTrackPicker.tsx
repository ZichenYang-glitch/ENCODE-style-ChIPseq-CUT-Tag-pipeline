/**
 * Interface seam for local tracks. IRENE-64 owns the browser reader that
 * replaces this file. Keep that implementation if the two commits conflict.
 * The panel passes a WashU scaffold such as hg38 and receives the full list
 * whenever the user adds or removes a file.
 */
export type LocalWashuTrack = {
  id: string;
  label: string;
  assembly: string;
  model: { type: string; name: string; [key: string]: unknown };
};

export default function LocalTrackPicker({
  assembly,
  onTracksChange,
}: {
  assembly: string;
  onTracksChange: (tracks: LocalWashuTrack[]) => void;
}) {
  void onTracksChange;
  return (
    <div data-testid="local-track-picker" data-assembly={assembly}>
      <p role="status">
        Local tracks stay in this browser. Choose files here; nothing is uploaded.
      </p>
    </div>
  );
}

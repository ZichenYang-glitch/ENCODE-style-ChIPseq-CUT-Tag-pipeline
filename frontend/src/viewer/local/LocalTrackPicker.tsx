import { useEffect, useRef, useState } from 'react'
import {
  PLAIN_PEAK_MAX_BYTES,
  assemblySignature,
  buildLocalTrack,
  type LocalWashuTrack,
} from './contract'
import {
  classifyLocalFiles,
  type LocalSelectionChoice,
  type LocalSelectionIssue,
} from './classifyLocalFiles'

export type { LocalWashuTrack }
export { PLAIN_PEAK_MAX_BYTES, assemblySignature }

interface OpenChoice extends LocalSelectionChoice {
  assembly: string
}

export default function LocalTrackPicker({
  assembly,
  onTracksChange,
}: {
  assembly: string
  onTracksChange: (tracks: LocalWashuTrack[]) => void
}) {
  const [tracks, setTracks] = useState<LocalWashuTrack[]>([])
  const [issues, setIssues] = useState<LocalSelectionIssue[]>([])
  const [choices, setChoices] = useState<OpenChoice[]>([])
  const [busy, setBusy] = useState(false)
  const request = useRef(0)
  const onTracksChangeRef = useRef(onTracksChange)
  onTracksChangeRef.current = onTracksChange

  useEffect(() => {
    onTracksChangeRef.current(tracks)
  }, [tracks])

  async function acceptFiles(list: FileList | null) {
    if (!list || list.length === 0) return
    const ticket = ++request.current
    const selectedAssembly = assembly
    const files = Array.from(list)
    setBusy(true)
    try {
      const result = await classifyLocalFiles(files)
      if (ticket !== request.current) return
      setIssues(result.issues)
      if (result.choices.length > 0) {
        setChoices((current) => [
          ...current,
          ...result.choices.map((choice) => ({
            ...choice,
            assembly: selectedAssembly,
          })),
        ])
      }
      if (result.ready.length > 0) {
        setTracks((current) => [
          ...current,
          ...result.ready.map((prepared) => buildLocalTrack(prepared, selectedAssembly)),
        ])
      }
    } finally {
      if (ticket === request.current) setBusy(false)
    }
  }

  function confirmChoice(choice: OpenChoice, optionId: string) {
    const option = choice.options.find((item) => item.id === optionId)
    setChoices((current) => current.filter((item) => item.id !== choice.id))
    if (!option) return
    setTracks((current) => [
      ...current,
      buildLocalTrack(option.prepared, choice.assembly),
    ])
  }

  function dismissChoice(choice: OpenChoice) {
    setChoices((current) => current.filter((item) => item.id !== choice.id))
  }

  const signature = assemblySignature(assembly)

  return (
    <section className="min-w-0 space-y-3" data-testid="local-track-picker" aria-label="Local tracks">
      <p data-testid="local-track-assembly">{signature}</p>
      <p className="text-sm">
        文件本身不能证明组装。已添加的轨道保持添加时的参考，不会随当前选择改标。
      </p>
      <ul className="text-sm" data-testid="local-track-formats" data-plain-max-bytes={PLAIN_PEAK_MAX_BYTES}>
        <li data-format="longrange">预建 longrange：.longrange.gz + .tbi，按窗口读取，不转换整份。</li>
        <li data-format="bigwig">bigWig：.bw / .bigWig，按索引窗口读取。</li>
        <li data-format="bed">
          峰：BED3/BED4、narrowPeak、broadPeak，不超过 8 MiB；.bed.gz + .tbi 按窗口读取。
        </li>
      </ul>
      <label className="block text-sm">
        添加本地轨道
        <input
          className="mt-1 block w-full max-w-md text-sm"
          data-testid="local-track-file-input"
          aria-label="Local track files"
          type="file"
          multiple
          onChange={(event) => {
            const list = event.target.files
            event.target.value = ''
            void acceptFiles(list)
          }}
        />
      </label>
      {busy ? (
        <p role="status" data-testid="local-track-busy">
          正在检查文件头…
        </p>
      ) : null}
      {issues.map((item) => (
        <div key={item.id} role="alert" data-reason={item.reason}>
          <p>{item.message}</p>
          <button
            type="button"
            aria-label="Dismiss local track message"
            onClick={() => setIssues((current) => current.filter((issue) => issue.id !== item.id))}
          >
            关闭
          </button>
        </div>
      ))}
      {choices.map((choice) => (
        <div key={choice.id} data-testid="local-track-choice" data-reason={choice.reason}>
          <p>{choice.message}</p>
          <p>确认后记为 {assemblySignature(choice.assembly)}，不会跟着之后的参考改标。</p>
          <div className="flex flex-wrap gap-2">
            {choice.options.map((option) => (
              <button
                key={option.id}
                type="button"
                data-choice={option.prepared.washUType}
                onClick={() => confirmChoice(choice, option.id)}
              >
                {option.label}
              </button>
            ))}
            <button type="button" data-choice="skip" onClick={() => dismissChoice(choice)}>
              不添加
            </button>
          </div>
        </div>
      ))}
      <ul className="space-y-2">
        {tracks.map((track) => {
          const readMethod = track.model.metadata.readMethod
          const mismatched = track.assembly.trim().toLowerCase() !== assembly.trim().toLowerCase()
          return (
            <li
              key={track.id}
              data-testid="local-track"
              data-track-id={track.id}
              data-track-type={track.model.type}
              data-read-method={readMethod}
              data-format={track.model.metadata.format}
              data-assembly={track.assembly}
            >
              <p>{track.label}</p>
              <p>{assemblySignature(track.assembly)}</p>
              <p>{readMethodLabel(readMethod)}</p>
              {mismatched ? (
                <p role="status" data-reason="assembly-not-relabeled">
                  {track.label} 按 {assemblySignature(track.assembly)} 添加。当前参考是 {signature}
                  ，未改标。
                </p>
              ) : null}
              <button
                type="button"
                aria-label={`Remove local track ${track.label}`}
                onClick={() => setTracks((current) => current.filter((item) => item.id !== track.id))}
              >
                移除
              </button>
            </li>
          )
        })}
      </ul>
    </section>
  )
}

function readMethodLabel(readMethod: string): string {
  if (readMethod === 'wuepgg3-track-text-file') return '整文件读取，受 8 MiB 上限约束'
  return '按窗口读取，不上传'
}

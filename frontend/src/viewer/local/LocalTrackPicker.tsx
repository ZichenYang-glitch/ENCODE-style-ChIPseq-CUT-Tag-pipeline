import { useEffect, useRef, useState } from 'react'
import {
  PLAIN_PEAK_MAX_BYTES,
  assemblyKey,
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
  checking?: boolean
  problem?: LocalSelectionIssue
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

  async function acceptFiles(files: File[]) {
    if (!files.length) return
    const ticket = ++request.current
    const selectedAssembly = assembly
    setBusy(true)
    try {
      const result = await classifyLocalFiles(files, selectedAssembly)
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

  async function confirmChoice(choice: OpenChoice, optionId: string) {
    const option = choice.options.find((item) => item.id === optionId)
    if (!option || choice.checking) return
    setChoices(current => current.map(item => item.id === choice.id ? { ...item, checking: true } : item))
    const prepared = await option.prepare()
    if ('issue' in prepared) {
      setChoices(current => current.map(item => item.id === choice.id ? { ...item, checking: false, problem: prepared.issue } : item))
      return
    }
    setChoices((current) => current.filter((item) => item.id !== choice.id))
    setTracks((current) => [
      ...current,
      buildLocalTrack(prepared, choice.assembly),
    ])
  }

  function dismissChoice(choice: OpenChoice) {
    setChoices((current) => current.filter((item) => item.id !== choice.id))
  }

  const signature = assemblySignature(assembly)

  return (
    <section className="min-w-0 space-y-3" data-testid="local-track-picker" aria-label="Local tracks">
      <h2 className="font-medium">打开本地轨道</h2>
      <p data-testid="local-track-assembly">{signature}</p>
      <p className="text-sm">
        请选择文件实际使用的组装。缺少组装声明时以您的选择登记；文件留在此浏览器。
      </p>
      <ul className="text-sm" data-testid="local-track-formats" data-plain-max-bytes={PLAIN_PEAK_MAX_BYTES}>
        <li data-format="longrange">预建 longrange 使用 .longrange.gz 和配套的 BGZF 压缩 .tbi，按窗口读取。</li>
        <li data-format="bigwig">bigWig（.bw 或 .bigWig）按索引窗口读取。</li>
        <li data-format="bed">
          BED3/BED4、narrowPeak、broadPeak 文本上限为 8 MiB。带 .tbi 的 BGZF 峰文件按窗口读取。
        </li>
        <li>本地 .tbi 解压后的上限为 32 MiB。</li>
      </ul>
      <label className="block text-sm">
        添加本地轨道
        <input
          className="mt-1 block w-full max-w-md text-sm"
          data-testid="local-track-file-input"
          aria-label="Local track files"
          type="file"
          multiple
          disabled={busy}
          onChange={(event) => {
            const files = Array.from(event.currentTarget.files ?? [])
            event.currentTarget.value = ''
            void acceptFiles(files)
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
          {choice.problem && <p role="alert" data-reason={choice.problem.reason}>{choice.problem.message}</p>}
          <div className="flex flex-wrap gap-2">
            {choice.options.map((option) => (
              <button
                key={option.id}
                type="button"
                data-choice={option.washUType}
                disabled={choice.checking}
                onClick={() => void confirmChoice(choice, option.id)}
              >
                {option.label}
              </button>
            ))}
            <button type="button" data-choice="skip" disabled={choice.checking} onClick={() => dismissChoice(choice)}>
              不添加
            </button>
          </div>
        </div>
      ))}
      <ul className="space-y-2">
        {tracks.map((track) => {
          const readMethod = track.model.metadata.readMethod
          const mismatched = assemblyKey(track.assembly) !== assemblyKey(assembly)
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

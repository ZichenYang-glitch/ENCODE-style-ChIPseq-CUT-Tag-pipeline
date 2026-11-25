"""Synthetic local tracks built with the locked pysam and pyBigWig.

    python test/browser/local_overlay_fixtures.py --dest DIR

The chr1 filler uses a fixed seed and is large enough to pass both online
converter caps. chr7 anchors, the peak and the 4.5 signal match the fixture
that already drew in GenomeHub. chr17 adds one TP53 pair, peak and signal.
A second indexed file supplies the foreign TBI used as the negative case.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
from random import Random
import time

ONLINE_SOURCE_BYTES = 32 * 1024 * 1024
ONLINE_SOURCE_RECORDS = 2_000_000
DEFAULT_FILLER_PAIRS = 1_010_000

_ROOT = Path(__file__).resolve().parents[2]
_CONSENSUS = _ROOT / "test" / "fixtures" / "browser" / "local-overlay-consensus.json"


def load_consensus() -> dict:
    return json.loads(_CONSENSUS.read_text(encoding="utf-8"))


def pair_lines(chrom: str, pair: dict) -> list[str]:
    start, end = pair["a"]
    mate_start, mate_end = pair["b"]
    return [
        f"{chrom}\t{start}\t{end}\t{chrom}:{mate_start}-{mate_end},1",
        f"{chrom}\t{mate_start}\t{mate_end}\t{chrom}:{start}-{end},1",
    ]


def fetch_window(
    data: Path,
    index: Path,
    chrom: str,
    start: int,
    end: int,
) -> list[str]:
    import pysam

    with pysam.TabixFile(str(data), index=str(index)) as table:
        return list(table.fetch(chrom, start, end))


def _write_longrange(path: Path, spec: dict, filler_pairs: int) -> int:
    import pysam

    random = Random(spec["seed"])
    pending: list[str] = []
    count = 0
    with pysam.BGZFile(str(path), "wb") as out:
        for index in range(filler_pairs):
            start = index * 100 + 100
            mate = start + 45
            tag = f"{random.getrandbits(128):032x}"
            pending.append(
                f"chr1\t{start}\t{start + 10}\tchr1:{mate}-{mate + 10},1\t{tag}\n"
            )
            pending.append(
                f"chr1\t{mate}\t{mate + 10}\tchr1:{start}-{start + 10},1\t{tag}\n"
            )
            count += 2
            if len(pending) >= 8192:
                out.write("".join(pending).encode())
                pending.clear()
        for key in ("tp53", "egfr"):
            chrom = spec["genes"]["TP53" if key == "tp53" else "EGFR"]["chrom"]
            for line in pair_lines(chrom, spec["longrange_pairs"][key][0]):
                pending.append(line + "\n")
                count += 1
        if pending:
            out.write("".join(pending).encode())
    pysam.tabix_index(str(path), preset="bed", force=True)
    return count


def _write_other(path: Path, line: str) -> None:
    import pysam

    with pysam.BGZFile(str(path), "wb") as out:
        out.write((line + "\n").encode())
    pysam.tabix_index(str(path), preset="bed", force=True)


def _write_signal(path: Path, spec: dict) -> list[list]:
    import pyBigWig

    sizes = spec["chrom_sizes_hg38"]
    with pyBigWig.open(str(path), "w") as writer:
        writer.addHeader(
            [(item["chrom"], sizes[item["chrom"]]) for item in spec["signals"]]
        )
        for item in spec["signals"]:
            writer.addEntries(
                [item["chrom"]],
                [item["start"]],
                ends=[item["end"]],
                values=[float(item["value"])],
            )
    confirmed: list[list] = []
    with pyBigWig.open(str(path)) as reader:
        for item in spec["signals"]:
            rows = reader.intervals(item["chrom"], item["start"], item["end"]) or []
            confirmed.extend(
                [item["chrom"], start, end, value] for start, end, value in rows
            )
    return confirmed


def _write_peaks(path: Path, spec: dict) -> str:
    text = "".join(
        f"{item['chrom']}\t{item['start']}\t{item['end']}\t{item['name']}\n"
        for item in spec["peaks"]
    )
    path.write_text(text, encoding="ascii")
    return text


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def build_local_overlay_fixtures(
    dest: Path,
    *,
    filler_pairs: int = DEFAULT_FILLER_PAIRS,
) -> dict:
    spec = load_consensus()
    dest.mkdir(parents=True, exist_ok=True)
    names = spec["filenames"]
    data = dest / names["longrange"]
    started = time.perf_counter()
    record_count = _write_longrange(data, spec, filler_pairs)
    index = dest / names["longrange_index"]
    window = spec["windows"]["egfr"]
    query_started = time.perf_counter()
    hits = fetch_window(data, index, window["chrom"], window["start"], window["end"])
    query_ms = (time.perf_counter() - query_started) * 1000
    tp53 = spec["windows"]["tp53"]
    tp53_hits = fetch_window(data, index, tp53["chrom"], tp53["start"], tp53["end"])
    other = dest / names["other_longrange"]
    _write_other(other, spec["other_longrange_line"])
    other_index = dest / names["other_index"]
    try:
        foreign_hits = fetch_window(
            data, other_index, window["chrom"], window["start"], window["end"]
        )
        wrong_index = {"raised": False, "hits": foreign_hits}
    except ValueError as error:
        wrong_index = {"raised": True, "error": str(error), "hits": []}
    signal = _write_signal(dest / names["signal"], spec)
    peaks = _write_peaks(dest / names["peaks"], spec)
    data_bytes = data.stat().st_size
    index_bytes = index.stat().st_size
    observation = {
        "filenames": names,
        "seed": spec["seed"],
        "filler_pairs": filler_pairs,
        "record_count": record_count,
        "feature_records": 4,
        "data_file_bytes": data_bytes,
        "index_file_bytes": index_bytes,
        "signal_file_bytes": (dest / names["signal"]).stat().st_size,
        "peaks_file_bytes": (dest / names["peaks"]).stat().st_size,
        "online_source_byte_limit": ONLINE_SOURCE_BYTES,
        "online_source_record_limit": ONLINE_SOURCE_RECORDS,
        "exceeds_source_bytes": data_bytes > ONLINE_SOURCE_BYTES,
        "exceeds_source_records": record_count > ONLINE_SOURCE_RECORDS,
        "data_sha256": _sha256(data),
        "index_sha256": _sha256(index),
        "build_seconds": round(time.perf_counter() - started, 3),
        "window": {
            "region": f"{window['chrom']}:{window['start']}-{window['end']}",
            "hits": hits,
            "elapsed_ms": round(query_ms, 3),
        },
        "tp53_window": {
            "region": f"{tp53['chrom']}:{tp53['start']}-{tp53['end']}",
            "hits": tp53_hits,
        },
        "signal_intervals": signal,
        "peaks_text": peaks,
        "wrong_index": wrong_index,
    }
    (dest / "observation.json").write_text(
        json.dumps(observation, indent=2) + "\n",
        encoding="utf-8",
    )
    return observation


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dest", required=True, type=Path)
    parser.add_argument("--filler-pairs", type=int, default=DEFAULT_FILLER_PAIRS)
    args = parser.parse_args()
    spec = load_consensus()
    observation = build_local_overlay_fixtures(
        args.dest, filler_pairs=args.filler_pairs
    )
    expected = pair_lines("chr7", spec["longrange_pairs"]["egfr"][0])
    if observation["window"]["hits"] != expected:
        raise SystemExit("index spot-check missed the EGFR anchors")
    if observation["wrong_index"]["raised"] is not True:
        raise SystemExit("foreign TBI was accepted for the EGFR window")
    if args.filler_pairs >= DEFAULT_FILLER_PAIRS and not (
        observation["exceeds_source_bytes"] and observation["exceeds_source_records"]
    ):
        raise SystemExit("synthetic longrange does not exceed the online converter cap")
    summary = {
        "records": observation["record_count"],
        "data_file_bytes": observation["data_file_bytes"],
        "index_file_bytes": observation["index_file_bytes"],
        "window_hits": len(observation["window"]["hits"]),
        "tp53_hits": len(observation["tp53_window"]["hits"]),
        "window_elapsed_ms": observation["window"]["elapsed_ms"],
        "wrong_index_raised": observation["wrong_index"]["raised"],
    }
    print("LOCAL_OVERLAY_FIXTURE " + json.dumps(summary))


if __name__ == "__main__":
    main()

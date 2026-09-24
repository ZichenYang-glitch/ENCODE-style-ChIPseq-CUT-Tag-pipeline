#!/usr/bin/env python3
"""确定性 Hi-TrAC 微型输入；源端点是设计值，必须以真实比对核对。

只生成合成序列，不运行工具，不读用户数据。输出目录必须不存在。
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import random
from pathlib import Path

SEED = 20260924
LINKER = "CTGTCTCTTATACACATCT"


def rc(sequence: str) -> str:
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def linker_position(sequence: str) -> int:
    """固定原脚本的设计预期，独立运行时须再与原脚本比较。"""
    prefixes = (LINKER[:9], rc(LINKER)[:9])
    return next(
        (i for i in range(len(sequence) - 9) if sequence[i : i + 9] in prefixes), -1
    )


def reference() -> tuple[dict[str, str], int]:
    rng = random.Random(SEED)
    motifs = {LINKER[:9], rc(LINKER)[:9], rc(LINKER[:9]), rc(rc(LINKER)[:9])}
    for attempt in range(1, 1001):
        sequences = {
            name: "".join(rng.choices("ACGT", k=100_000)) for name in ("chrA", "chrB")
        }
        # 真正同序列的双位点参考，用于验证低 MAPQ；不预设 Bowtie2 实际 MAPQ。
        sequences["chrB"] = (
            sequences["chrB"][:70_000]
            + sequences["chrA"][70_000:72_000]
            + sequences["chrB"][72_000:]
        )
        if all(all(motif not in seq for motif in motifs) for seq in sequences.values()):
            return sequences, attempt
    raise RuntimeError("未找到无意外 linker seed 的参考")


def end(chrom: str, start: int, length: int = 120, strand: str = "+") -> dict:
    return {"chrom": chrom, "start": start, "end": start + length, "strand": strand}


def make_pair(
    ref: dict[str, str],
    name: str,
    first: dict | None,
    second: dict | None,
    suffix1: str = "",
    suffix2: str = "",
    note: str = "",
    lane: str | None = None,
) -> dict:
    seqs = []
    for coordinate, suffix in ((first, suffix1), (second, suffix2)):
        if coordinate is None:
            seq = "N" * 120
        else:
            seq = ref[coordinate["chrom"]][coordinate["start"] : coordinate["end"]]
            if coordinate["strand"] == "-":
                seq = rc(seq)
        seqs.append(seq + suffix)
    positions = [linker_position(seq) for seq in seqs]
    trimmed = [seq if pos < 0 else seq[:pos] for seq, pos in zip(seqs, positions)]
    result = {
        "id": name,
        "sequences": seqs,
        "source_endpoints": [first, second],
        "expected_linker_positions": positions,
        "expected_trim_lengths": list(map(len, trimmed)),
        "expected_trim_retained": all(len(seq) >= 10 for seq in trimmed),
        "design_note": note,
        "lane_before_user_premerge": lane,
    }
    if first and second and first["chrom"] == second["chrom"]:
        mids = [(item["start"] + item["end"]) / 2 for item in (first, second)]
        result["designed_float_distance"] = abs(mids[1] - mids[0])
        result["designed_qc_integer_distance"] = abs(int(mids[1]) - int(mids[0]))
    return result


def baseline(ref: dict[str, str]) -> list[dict]:
    def p(name, first, second, linker=False, note="", lane="lane1"):
        return make_pair(
            ref,
            name,
            end("chrA", first),
            end("chrA", second, strand="-"),
            LINKER if linker else "",
            note=note,
            lane=lane,
        )

    rows = [
        p("a1", 1000, 1300, True, "cis300，有 linker；与 a2 同端点", "lane1"),
        p("a2", 1000, 1300, True, "跨设计 lane 的同样本重复", "lane2"),
        p("b", 3000, 3300, False, "cis300，无 linker；背景应删除"),
        p("c", 5000, 6000, False, "cis1000；背景保留且 QC close"),
        p("d", 10000, 15000, False, "cis5000；middle", "lane2"),
        p("e", 20000, 35000, False, "cis15000；distal", "lane2"),
        make_pair(
            ref,
            "f",
            end("chrA", 40000),
            end("chrB", 1000, strand="-"),
            LINKER,
            note="trans 有 linker；保留",
            lane="lane1",
        ),
        make_pair(
            ref,
            "g",
            end("chrA", 45000),
            end("chrB", 5000, strand="-"),
            note="trans 无 linker；背景删除",
            lane="lane2",
        ),
    ]
    assert [row["expected_linker_positions"] for row in rows] == [
        [120, -1],
        [120, -1],
        [-1, -1],
        [-1, -1],
        [-1, -1],
        [-1, -1],
        [120, -1],
        [-1, -1],
    ]
    return rows


def write_pair(directory: Path, sample: str, rows: list[dict]) -> dict:
    directory.mkdir(parents=True, exist_ok=True)
    files = []
    for mate in (1, 2):
        path = directory / f"{sample}_R{mate}.fastq.gz"
        # 固定 mtime 且无 filename 字段，使生成器输出可重复；验收不要求 gzip 字节一致。
        with (
            path.open("xb") as raw,
            gzip.GzipFile(fileobj=raw, mode="wb", mtime=0, filename="") as handle,
        ):
            for row in rows:
                sequence = row["sequences"][mate - 1]
                record = f"@{sample}_{row['id']}/{mate}\n{sequence}\n+\n{'I' * len(sequence)}\n"
                handle.write(record.encode("ascii"))
        files.append(
            {
                "path": str(path),
                "bytes": path.stat().st_size,
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    return {
        "sample": sample,
        "raw_pairs": len(rows),
        "files": files,
        "pairs": [
            {key: value for key, value in row.items() if key != "sequences"}
            for row in rows
        ],
    }


def generate(root: Path) -> dict:
    root = root.absolute()
    root.mkdir(parents=True, exist_ok=False)
    ref, reference_draws = reference()
    fasta = root / "reference.fa"
    with fasta.open("x", encoding="ascii") as handle:
        for name, sequence in ref.items():
            handle.write(f">{name}\n")
            handle.write(
                "\n".join(sequence[i : i + 80] for i in range(0, len(sequence), 80))
                + "\n"
            )
    base = baseline(ref)
    examples: dict[str, dict] = {}

    def scenario(
        name: str,
        samples: list[list[dict]],
        description: str,
        policy: str = "needs_real_mapping",
    ):
        examples[name] = {
            "description": description,
            "policy_expectation": policy,
            "samples": [
                write_pair(root / name / "fastq", f"s{i:06d}", rows)
                for i, rows in enumerate(samples, 1)
            ],
        }

    scenario(
        "positive",
        [base, base],
        "两个独立样本各八对；每样本已由用户预合 lane，一对 gzip 输入",
        "success_if_actual_alignment_matches_design",
    )
    examples["positive"]["conditional_expected_summary"] = [
        8,
        8,
        1,
        8,
        0.125,
        5 / 7,
        3 / 5,
        1 / 5,
        1 / 5,
        5,
        5 / 8,
        4 / 5,
        1 / 2,
        1 / 4,
        1 / 4,
    ]
    examples["positive"]["conditional_expected_counts"] = {
        "all": 8,
        "all_qc_unique": 7,
        "all_qc_cis": 5,
        "noBg": 5,
        "noBg_qc_cis": 4,
    }
    links = [base[4]]
    specs = [
        ("single", 51000, 120, LINKER, "", [120, -1]),
        ("both", 52000, 120, LINKER, LINKER, [120, 120]),
        ("reverse_complement", 53000, 120, rc(LINKER), "", [120, -1]),
        ("multiple", 54000, 120, LINKER + "GATTACAGAT" + LINKER, "", [120, -1]),
        ("last_nine_seed", 55000, 120, LINKER[:9], "", [-1, -1]),
        ("last_nine_rc_seed", 56000, 120, rc(LINKER)[:9], "", [-1, -1]),
        ("trim_to_nine", 57000, 9, LINKER, "", [9, -1]),
        ("trim_to_ten", 58000, 10, LINKER, "", [10, -1]),
    ]
    for name, start, length, suffix1, suffix2, expected in specs:
        pair = make_pair(
            ref,
            name,
            end("chrA", start, length),
            end("chrA", start + 300, strand="-"),
            suffix1,
            suffix2,
            note="只预期 trim；短 reads 和末尾 seed 的实际比对须实测",
        )
        assert pair["expected_linker_positions"] == expected, (
            name,
            pair["expected_linker_positions"],
        )
        links.append(pair)
    distances = [base[4]]
    for distance, start in (
        (999, 60000),
        (1000, 62000),
        (1001, 64000),
        (10000, 1000),
        (10001, 1500),
    ):
        distances.append(
            make_pair(
                ref,
                f"distance_{distance}",
                end("chrA", start),
                end("chrA", start + distance, strand="-"),
                note="无 linker；区分背景与 QC 边界",
            )
        )
    for name, first, second in (
        ("odd_999_5", 66000, 67000),
        ("odd_1000_5", 68000, 69001),
    ):
        distances.append(
            make_pair(
                ref,
                name,
                end("chrA", first, 121),
                end("chrA", second, strand="-"),
                note="不同读长导致半整数距离；QC 使用整数中点",
            )
        )
    low_mapq = make_pair(
        ref,
        "duplicated_reference",
        end("chrA", 70200),
        end("chrA", 70500, strand="-"),
        LINKER,
        note="chrB:70000-72000 与 chrA 同序列；必须核实际 MAPQ，不能凭设计断言",
    )
    unmapped = make_pair(
        ref, "unmapped", None, None, note="两端各120个N；真实工具核对未比对"
    )
    orphan1 = make_pair(
        ref,
        "only_mate1_mapped",
        end("chrA", 74000),
        None,
        LINKER,
        note="真实工具核对 orphan BAM/bamToBed 行为，不预设过滤后一定零 PET",
    )
    orphan2 = make_pair(
        ref,
        "only_mate2_mapped",
        None,
        end("chrA", 75000, strand="-"),
        suffix2=LINKER,
        note="mate1未比对，mate2可比对；核 bedtools 处理",
    )
    scenario(
        "boundaries",
        [links, distances, [base[4], low_mapq, unmapped, orphan1, orphan2]],
        "三样本分别核 linker/长度、距离、低MAPQ/未比对/orphan；各保留一个有效 cis anchor",
    )
    scenario("empty", [[]], "合法 gzip 内零 FASTQ 记录", "reject_whole_attempt")
    scenario("all_trimmed", [[links[-2]]], "唯一pair修剪至9bp", "reject_whole_attempt")
    scenario("all_unmapped", [[unmapped]], "全部未比对", "reject_whole_attempt")
    scenario(
        "noBg_empty",
        [[base[2]]],
        "有all cis，但无linker的短cis被背景过滤",
        "reject_whole_attempt",
    )
    scenario(
        "all_trans",
        [[base[6]]],
        "有linker的trans；all和noBg无cis",
        "reject_whole_attempt",
    )
    scenario(
        "all_trans_no_linker",
        [[base[7]]],
        "trans无linker；all无cis且noBg空",
        "reject_whole_attempt",
    )
    scenario(
        "only_low_mapq",
        [[low_mapq]],
        "重复参考导致低MAPQ候选；须核实际工具结果",
        "needs_real_mapping",
    )
    scenario(
        "singlemate_only",
        [[orphan1, orphan2]],
        "完整FASTQ配对但仅单mate能比对；核实际BAM和BEDPE",
        "needs_real_mapping",
    )
    scenario(
        "multisample_one_empty",
        [base, [base[2]]],
        "第二样本noBg空，验证整attempt拒绝而非部分成功",
        "reject_whole_attempt",
    )
    metadata = {
        "schema": 1,
        "seed": SEED,
        "reference_draws": reference_draws,
        "upstream_script_sha256": "c3c4ef4e6287fa4ade97a6f5980d20c81ea345b8e8e13c88ae3b7c4bd3a67aec",
        "reference": {
            "path": str(fasta),
            "sha256": hashlib.sha256(fasta.read_bytes()).hexdigest(),
            "contigs": {key: len(value) for key, value in ref.items()},
            "duplicated_region": "chrA:70000-72000 equals chrB:70000-72000",
        },
        "evidence_boundary": "仅生成合成输入与手推预期，未执行bowtie2-build/tracPre2；真实MAPQ/端点/QC是后续证据",
        "scenarios": examples,
    }
    (root / "design.json").write_text(
        json.dumps(metadata, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    return metadata


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    result = generate(args.output)
    print(
        json.dumps(
            {
                "reference_sha256": result["reference"]["sha256"],
                "reference_draws": result["reference_draws"],
                "scenarios": list(result["scenarios"]),
            },
            ensure_ascii=False,
        )
    )

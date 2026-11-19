# 同事试用：从已核验样本找到 PET 和 QC

状态：试用单设计。写给一位没有参与本预处理实现、本文件、以及该批次计算或核验的同事。
完成下表只说明这位同事从受权入口找到了目标 PET、说清了 QC 分母，并打开了自己勾选的
WashU 轨道。它不是用户签字，不是整批发布，也不是生物学确认。

设计对照的 main 是 `852202096db519a5e6ffdb39418fdb14e371887e`（2026-11-19）。
当时该人源 Hi-TrAC 预处理批次的正式结果包尚未发布，执行入口是本机资格脚本，
没有平台 run，也没有外部结果修订。数据主人尚未授权同事打开该批次。因此本文件
只规定动作；下面没有填写过的通过记录。

## 复用的专题

分母、产物类型和轨道含义沿用已有文档。本单不另定义一套。

| 同事要做的事 | 沿用 |
| --- | --- |
| QC 十五列的分子和分母 | [预处理契约第 3 节](hitrac-preprocess-contract.md#3-pet输出和-qc-单位) |
| 公开产物、政策 A、内部样本 token | [操作说明第 4 节](hitrac-preprocess-operations.md#4-结果输出与通知契约) |
| 上游 `tracPre2.py` 身份 | [操作说明第 1 节](hitrac-preprocess-operations.md#1-固定身份)：cLoops2 `de6cc732fa00b408551b9f4272933640c08447f1`，脚本 SHA256 `c3c4ef4e6287fa4ade97a6f5980d20c81ea345b8e8e13c88ae3b7c4bd3a67aec` |
| PET 在 WashU 里是什么 | [WashU 只读契约「科学语义」](washu-browser-read-contract.md#科学语义) |
| 从项目找到实验 | [实验目录 M0](../architecture/experiment-catalog-m0-contract.md) 的实验追溯信封。目录页实现在复审中的 [PR #235](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/235)（`37e04b9582a390712db5f116eaa1498a16c90c87`） |
| checkout 产物的身份 | [PR #228](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/228)（`ab9bbd1a7957defbc98929730c9e0af93767ddac`）的 `rreg_` / `rrev_`。该契约不在本提交的 main 上 |
| 在 WashU 里勾选、移除轨道 | 设计 [PR #227](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/227)（`fada3ed129dfcd014a3c7825c77c315667faf978`）；实现复审 [PR #237](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/237)（`9c162a7bdda3ce182a44baa7617f69d041e87683`）。本分支的 WashU 文没有自选一节 |
| 软件行为与科学结论分开 | [PR #231](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/231)（`15650d77a3ecc16c87f636d384f7c18320857520`） |

平台 QC 行的 `metric_key` 以 `src/encode_pipeline/adapters/hitrac_preprocess/results.py`
里的 `METRIC_KEYS` 为准：`raw_pairs`、`trimmed_pairs`、`mapping_ratio`、`all_pets`、
`all_redundancy`、`all_cis_ratio`、`all_close_ratio`、`all_middle_ratio`、
`all_distal_ratio`、`no_bg_pets`、`yield`、`no_bg_cis_ratio`、`no_bg_close_ratio`、
`no_bg_middle_ratio`、`no_bg_distal_ratio`。来源产物类型是 `hitrac_summary`。

## 人和数据前置

三条都成立之后才从第 1 步做起。缺任何一条，就把记录停在这一条，试用未完成。

1. **人。** 执行人没有提交过 Hi-TrAC 预处理 adapter、本试用单，也没有担任该批次的计算或核验。开发者自己填表不算这次试用。
2. **授权。** 数据主人点名允许这位同事打开这一份已核验样本。允许看汇总 QC，不等于允许打开 PET 或 WashU。
3. **受权入口。** 样本能从项目里打开，并且入口自带核验证据。两条路只取实际存在的那条。

平台 run 的入口是：成员登录后，从项目打开实验；实验追溯信封的 `trace_status` 为
`experiment_snapshot`；该运行的 workflow id 是 `hitrac-preprocess`；同一 generation
里能看到 `hitrac_bedpe_all` 或 `hitrac_bedpe_no_bg`，以及 `hitrac_summary`。
发布记录写明这份样本的核验是原生通过，还是锁定版本校验器的重放通过。
`trace_status` 为 unresolved 时，项目里还没有准确的实验关联，试用停在这里。

只有 checkout 产物时，入口是外部结果登记：PET 文件和 summary 各自有 `rreg_` / `rrev_`，
`origin.kind` 为 `external`，修订上没有 `run_id`、artifact generation 或 publication。
方法写明 tracPre2 与上面的上游提交；参考写已知身份，或显式 `unknown`。
只登记了其中一个文件时，停下来。本机目录、绝对路径和开发者口头的文件位置都不是入口。
在 #228 的读取能力进入所登录的产品之前，这条路不可用，试用停在「还不能登记」。

核验结论随入口走。`pair_support_missing` 的样本在当前拒绝记录里，不作为这一版的目标 PET。
summary 里有一行数字，本身不是通过结论。原生通过和重放通过分两行写，不并成整批已发布。
数据主人若之后另行发布这些样本，或发布转换修正后的新版本，同事打开新的 generation
或新的 `rrev_`，记录只写打开的那一版。试用单不改判旧结论，也不把旧拒绝套到新版本上。
`verify_pairs`（`src/encode_pipeline/adapters/hitrac_preprocess/pairs.py`）是发布前的
双端来源检查：同一读名要有 R1 和 R2 两条比对，端点在染色体、起止和链向上与 BEDPE 相符，
两端次序可以互换。同事不重跑这个检查。

## 试用动作

每一步都先写下将要打开的版本，再看文件。版本对不上就停，不换一份继续。

### 1. 登记版本

记录这些实际打开的身份，缺一项就留空并标明缺的是哪一项：

- 所登录产品的 git SHA，以及本试用单的 git SHA
- 平台 run：`run_id`、result generation、所选 artifact id、参考修订 id 与 identity SHA-256
- 外部登记：PET 的 `rrev_`、summary 的 `rrev_`、各自内容 SHA-256、方法版本、参考状态
- 上游 tracPre2 提交 `de6cc732fa00b408551b9f4272933640c08447f1`
- WashU 页面实际加载的构建 SHA。自选行为以打开时的产品为准，并注明它是否包含 #237

### 2. 从项目找到目标 PET

平台 run：打开该实验绑定的运行，在产物里点选一份 PET。目标类型只在
`hitrac_bedpe_all`（all，MAPQ 过滤后的 BEDPE 行）和 `hitrac_bedpe_no_bg`
（去背景再按 BEDPE 前六列去重后的行）之中二选一。记下类型、artifact id 和 generation。
同一 generation 的 `hitrac_summary` 是配套 QC。BAM 和裁剪 FASTQ 保持私有，本试用不下载。

外部登记：从项目打开上面两个 `rrev_`。显示名只作辨认，身份是修订 id。
两个修订的方法版本和参考状态要能对上同一次资格；对不上就停。

### 3. 说明 QC 分母

打开 `/runs/<run_id>?view=qc`，或打开已登记的 summary。只解释与第 2 步同一 generation
或同一对 `rrev_` 的那一行。同事对着列名说出下表。符号 R、L、T、U、C、K、M、D、N、V、Cn
的定义在契约第 3 节，此处不改写。

| 原列 | `metric_key` | 同事要说出的分母 |
| --- | --- | --- |
| `total raw sequences` | `raw_pairs` | 配对数 R。一条是一对，不是 R1 与 R2 记录数之和 |
| `after linker removing sequences` | `trimmed_pairs` | 裁剪后仍成对的配对数 L |
| `mapping ratio` | `mapping_ratio` | Bowtie2 overall alignment rate 的百分数除以 100。分母不是 T/R |
| `total mapped PETs (mapq>=N)` | `all_pets` | 这列本身是分子 T：all 文件的 PET 行数。标题里的 N 是当次 MAPQ |
| `total mapped PETs redundancy` | `all_redundancy` | 1−U/T。分母是 T。U 是 all 按 QC 键的独特数，不是 noBg 行数 |
| `total mapped PETs intra-chromosomal ratio` | `all_cis_ratio` | C/U |
| `total mapped PETs close / middle / distal ratio` | `all_close_ratio`、`all_middle_ratio`、`all_distal_ratio` | K/C、M/C、D/C。分母是该 all 集合的 cis 数 C，不是 T，也不是 R |
| `unique noBg PETs` | `no_bg_pets` | 分子 N：noBg 文件的 PET 行数 |
| `Yield` | `yield` | N/R。分母是原始配对数 R |
| `unique noBg … intra-chromosomal ratio` | `no_bg_cis_ratio` | Cn/V |
| `unique noBg … close / middle / distal ratio` | `no_bg_close_ratio`、`no_bg_middle_ratio`、`no_bg_distal_ratio` | Kn/Cn、Mn/Cn、Dn/Cn。分母是 noBg 集合里的 cis 数 Cn |

另外说清两条边界，它们都写在契约第 3 节：

- 去背景看的是浮点中点距离 `<1000`，并且两端都没有 linker；跨染色体且两端都没有 linker 的行也丢。
- QC 的 close 是整数中点距离 `≤1000`，分母是对应集合的 cis 数。`<1000` 与 `≤1000` 是两条边界。

政策 A：任一样本 all 或 noBg 为空，或需要 cis 分母的那组 QC 没有 cis，整次 attempt 不发布成功。
上游在无 cis 时可能把分母写成 1；这种比例不当作已发布的计数。空结果也不表示该实验没有科学意义。

这一步的观察结果写成：列名、分子、分母、读到的数是否与契约的单位一致。
MAPQ 默认 10 与 threads 默认 2 的科学地位沿用 #231：它们是当时的软件默认，本试用不把它们冻结成实验室阈值。

### 4. 在 WashU 中只打开所选轨道

这一步要求登录的产品能从第 2 步的那次 run 打开 WashU，并且自选行为与 #227 的用例一致。
外部登记没有 `run_id`。#228 不提供 WashU 会话。在出现不冒充 run 的查看入口之前，
checkout 产物做到第 3 步为止，WashU 一行写「入口未提供」，试用未完成。
不要把本机文件路径填进查看器。

平台 run 上按打开的那个构建做：

1. 从当前 run 打开 WashU。直接打开时选择为空；勾选之前，轨道数据请求为零。
2. 只勾选第 2 步选定的那一份 PET。未勾选的文件留在列表里，不发起读取。
3. 对着画面说明：类型是 `longrange`；源文件是 PET；每条 PET 记 1 对；两条锚点的分数 `1` 是显示标志，不是 MAPQ，两行不能相加；`pair_count` 在衍生物来源摘要里。显著 loop 和 `.hic` 接触矩阵是别的产物类型。
4. 区域摘要超过 WashU 自身的 10 万行绘制上限时，页面会说明绘制不完整。记录这一句。画面不代替完整 PET 计数。
5. 去掉最后一条轨道后，选择回到空，后续不再读这些轨道。
6. 翻页时若代次变化，会话返回 409 `generation_changed`。先刷新并重选，再读新版文件。

若打开的构建在进入时已经勾选轨道，或未勾选的文件仍有数据请求，记下该构建 SHA，
这一步记为产品差距，试用未完成。

### 5. 留下脱敏记录

记录放在工单评论或数据主人指定的受权位置。字段如下。合成例子只用契约第 6 节的
手推基线（all 8 行、noBg 5 行、Yield 5/8），用来看格式；它不是某次真实文库，也不是通过。

| 字段 | 合成例子（不要照抄为真实结果） | 真实试用 |
| --- | --- | --- |
| 日期与执行人角色 | 设计示例；执行人未参与开发 |  |
| 入口 | 无。示例不登录 | `platform run` 或 `external registration` |
| 版本 | 契约第 6 节；上游 `de6cc732fa00b408551b9f4272933640c08447f1` | 第 1 步的 SHA、`run_id` 或两个 `rrev_` |
| 核验 | 合成基线，无批次裁决 | `native pass` 或 `replay pass` |
| 目标 PET | 类型 `hitrac_bedpe_no_bg`；契约预期 5 行 | 类型与 artifact id 或 `rrev_` |
| 分母 | `Yield` = N/R = 5/8；close 的分母是 cis 数 | 列名、分子、分母、与契约是否一致 |
| WashU | 未打开 | 初始选择为空；只读了所选 `longrange`；绘制是否完整 |
| 结论 | 只演示表格 | 找到 / 说清分母 / 打开轨道，或停在哪一条前置 |

记录里不写样本显示名、供方编号、读名、基因组坐标、FASTQ 或 BAM 内容、本机绝对路径、
凭据。内部 token 与用户样本名的对应留在平台 metadata，不抄进这份脱敏记录。

## 未执行

- 没有重跑该预处理批次，没有打开、复制、移动或改写原始 FASTQ、参考、BAM 或 BEDPE。
- 没有把 checkout 目录登记成平台 run，也没有代签用户验收。
- 没有在浏览器里完成第 4 步。#227 / #237 的自选实现、#228 的登记读取、#235 的目录页都不在本提交里。
- 核对命令只读本仓库文件。在仓库根目录执行后应打印 `hitrac-colleague-trial ok`：

```bash
python3 - <<'PY'
from pathlib import Path
root = Path(".")
trial = (root / "docs/development/hitrac-colleague-pet-qc-trial.md").read_text(encoding="utf-8")
prose = trial.split("```bash", 1)[0]
contract = (root / "docs/development/hitrac-preprocess-contract.md").read_text(encoding="utf-8")
washu = (root / "docs/development/washu-browser-read-contract.md").read_text(encoding="utf-8")
ops = (root / "docs/development/hitrac-preprocess-operations.md").read_text(encoding="utf-8")
results = (root / "src/encode_pipeline/adapters/hitrac_preprocess/results.py").read_text(encoding="utf-8")
assert "## 3. PET、输出和 QC 单位" in contract
assert "| `Yield` | N/R |" in contract
assert "1−U/T" in contract and "K/C" in contract
assert "hitrac-colleague-pet-qc-trial.md" in ops
assert "## 科学语义" in washu
assert "hitrac_bedpe_all" in washu and "hitrac_bedpe_no_bg" in washu
keys = (
    "raw_pairs", "trimmed_pairs", "mapping_ratio", "all_pets", "all_redundancy",
    "all_cis_ratio", "all_close_ratio", "all_middle_ratio", "all_distal_ratio",
    "no_bg_pets", "yield", "no_bg_cis_ratio", "no_bg_close_ratio",
    "no_bg_middle_ratio", "no_bg_distal_ratio",
)
for key in keys:
    assert f'"{key}"' in results
    assert f"`{key}`" in prose
for phrase in ("N/R", "1−U/T", "C/U", "K/C", "pair_support_missing", "rreg_", "rrev_", "experiment_snapshot", "generation_changed"):
    assert phrase in prose
for banned in ("s000", "/home/", "11_17", "helixweave-hitrac-preprocess."):
    assert banned not in prose
print("hitrac-colleague-trial ok")
PY
```

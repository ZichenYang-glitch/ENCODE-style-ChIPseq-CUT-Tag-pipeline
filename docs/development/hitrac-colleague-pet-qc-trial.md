# 同事试用：从已核验样本找到 PET 和 QC

状态：试用单设计。写给一位没有参与本预处理实现、本文件、以及该批次计算或核验的同事。
完成记录只说明这位同事从受权入口找到了获准样本的 PET、说清了该行的 QC 分母，并按当时构建
记下 WashU 的结果。轨道可以打开，也可以停在「已定位 PET，轨道未打开，试用未完成」。
它不是用户签字，不是整批发布，也不是生物学确认。

设计写于 2026-11-19，跟着它合并进去的 main；复审时以该 PR 的 base SHA 为准。
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
| checkout 产物的身份 | [外部结果读取契约](../architecture/registered-result-read-contract.md) 的 `rreg_` / `rrev_`。该文写明表、API 和客户端尚未实现 |
| 在 WashU 里勾选、移除轨道 | [WashU 契约「在 WashU 内选择轨道文件（IRENE-21）」](washu-browser-read-contract.md)。直接打开时选择为空，只读取勾选的文件 |
| 软件行为与科学结论分开 | [科学工作流证据台账](scientific-workflow-evidence-ledger.md) 的 Hi-TrAC 核对步骤。MAPQ 10 与 threads 2 在台账里仍是软件默认 |

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
读取契约写明相应表、API 和客户端尚未存在。登录的产品还不能按 `rrev_` 打开文件时，这条路不可用，试用停在「还不能登记」。

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
- WashU 页面实际加载的构建 SHA。自选行为对照契约「在 WashU 内选择轨道文件（IRENE-21）」

### 2. 从项目找到获准样本的 PET

一代可以有多个样本。all 与 noBg 按样本各一份，整批共用一份 `hitrac_summary`。
只按产物类型和 generation 点选，会读到另一个样本的行。

平台 run 在页面上核对三处，三处相同才继续：

1. 数据主人点名的获准样本。
2. 所选 PET 的 `metadata.sample_id`。产物表的 Sample 列就是它，是发布时的显示名，不是内部 token。类型仍只在 `hitrac_bedpe_all` 与 `hitrac_bedpe_no_bg` 里二选一。
3. QC 页 `/runs/<run_id>?view=qc` 上准备解释的那一行的 `sample_id`。

再看这一行的 `source_artifact_id`。点 View source，打开的必须是本代那份
`hitrac_summary`，且它的 `metadata.scope` 为 `run`。summary 文件本身不带样本身份。
三处有任一不一致、打开的不是本代 summary、或对应关系说不清，就停在这一步。
不按行序、文件名或 Yield 高低挑选。

外部登记没有平台 QC 行。受权核验记录必须把精确的 PET `rrev_` 指到 summary `rrev_`
里的具体一行。记录没有这层对应时停在这一步。显示名和行序都不是对应。

下面用同一代的两个合成样本演示排除错行。页面上的 `syn-a`、`syn-b` 只存在于这个例子，
不是真实文库，也不是通过。获准的是 `syn-a`。

| 页面对象 | 同代里能看见的 | 获准 `syn-a` 时留下什么 |
| --- | --- | --- |
| PET | 两份 `hitrac_bedpe_no_bg`。`metadata.sample_id` 分别是 `syn-a` 和 `syn-b` | 只留 `syn-a`。`syn-b` 的类型和 generation 相同，仍然排除 |
| summary | 一份 `hitrac_summary`，`metadata.scope` 为 `run` | 两行 QC 都指向它 |
| QC | 两行 `sample_id` 分别是 `syn-a` 和 `syn-b`，`source_artifact_id` 都是这份 summary 的 artifact id | 只解释 `syn-a` 那一行。另一行即使数字更好看也排除 |

脱敏回执不写 `syn-a` 或 `syn-b`，只写匹配结果：PET 的 `metadata.sample_id`、QC 的
`sample_id` 与获准样本一致；`source_artifact_id` 指向本代 summary；已排除同代另一行。
三处对不上时，回执写停在第 2 步，不解释那一行的分母。BAM 和裁剪 FASTQ 保持私有，本试用不下载。

### 3. 说明 QC 分母

只解释第 2 步已经对齐的那一行。同事对着列名说出下表。符号 R、L、T、U、C、K、M、D、N、V、Cn
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
当次 MAPQ 读该运行 options 里的 `mapq`，缺省 10，并与 summary 列标题中的 N 一致。
MAPQ 缺省 10 与 threads 缺省 2 的科学地位沿用证据台账：它们是软件默认，本试用不把它们冻结成实验室阈值。

### 4. 在 WashU 中只打开所选轨道

从第 2 步已经对齐的那份平台 run PET 打开 WashU。勾选、移除和代次变化按契约
「在 WashU 内选择轨道文件（IRENE-21）」。外部登记没有 `run_id`。外部结果契约不提供
WashU 会话。在出现不把登记结果写成 run 的查看入口之前，checkout 产物做到第 3 步为止，
WashU 一行写「入口未提供」，试用未完成。不要把本机文件路径填进查看器。

转换器 `washu-display-derivative-1` 对每份源文件最多接收 2,000,000 条记录
（`MAX_SOURCE_RECORDS`，`src/encode_pipeline/services/browser_derivatives.py`）。
超过就拒绝，原因是 `source_too_large`。页面写 “This file is too large to index here
and is not shown as empty data.” 空白画面不是空数据，也不是轨道已打开。
`test/services/test_browser_derivatives.py` 把源记录上限降到 1 时，两行合成 PET 得到
同一个 `source_too_large`。

IRENE-9 的脱敏审计里，10 个通过核验的样本中，最小的 all 文件有 25,799,851 条，最小的
noBg 文件有 15,829,716 条。两者都高于 2,000,000。当前构建不能为这批完整 PET 生成轨道。
即使以后取得访问授权，只要打开的仍是这些完整文件、上限仍是 2,000,000，记录就是
「已定位 PET，轨道未打开，试用未完成」。

平台 run 上按打开的那个构建做：

1. 从当前 run 打开 WashU。直接打开时选择为空；勾选之前，轨道数据请求为零。
2. 只勾选第 2 步对齐的那一份 PET。未勾选的文件留在列表里，不发起读取。
3. 文件还在上限以内、索引尚未准备时，文案是 “An index is not prepared yet. Choose this file, then prepare it.” 勾选后只点这一份的 `Prepare`。准备完成之前不把画面当成轨道。
4. 页面报 `source_too_large` 时不点准备，也不改读另一个样本。记录「已定位 PET，轨道未打开，试用未完成」。这一句就是本步结论。
5. 只有索引成功、轨道已经打开之后，才说明：类型是 `longrange`；源文件是 PET；每条 PET 记 1 对；两条锚点的分数 `1` 是显示标志，不是 MAPQ，两行不能相加；`pair_count` 在衍生物来源摘要里。显著 loop 和 `.hic` 接触矩阵是别的产物类型。区域摘要超过 WashU 自身的 10 万行绘制上限时，页面会说明绘制不完整。那句话发生在索引成功之后，不能用来描述超限而未打开的文件。画面不代替完整 PET 计数。
6. 去掉最后一条轨道后，选择回到空，后续不再读这些轨道。
7. 翻页时若代次变化，会话返回 409 `generation_changed`。先刷新并重选，再读新版文件。重新核对第 2 步的三处身份，不沿用上一份 PET。

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
| 样本对应 | 三者一致；已排除同代另一行。不写页面上的 sample id | 匹配成立，或停在第 2 步 |
| 目标 PET | 类型 `hitrac_bedpe_no_bg`；与获准样本的 `metadata.sample_id` 一致 | 类型与 artifact id 或 `rrev_`。回执只写是否匹配 |
| 分母 | `Yield` = N/R = 5/8；close 的分母是 cis 数 | 列名、分子、分母、与契约是否一致。只解释对齐的那一行 |
| WashU | 上限以内的小文件先点 `Prepare`。本批完整 PET 不走这条 | 「已定位 PET，轨道未打开，试用未完成」，或索引成功后的轨道观察 |
| 结论 | 只演示如何排除错行 | 找到且说清分母；轨道已打开，或已定位但轨道未打开 |

记录里不写样本显示名、供方编号、读名、基因组坐标、FASTQ 或 BAM 内容、本机绝对路径、
凭据。内部 token 与用户样本名的对应留在平台 metadata，不抄进这份脱敏记录。

## 未执行

- 没有重跑该预处理批次，没有打开、复制、移动或改写原始 FASTQ、参考、BAM 或 BEDPE。
- 没有把 checkout 目录登记成平台 run，也没有代签用户验收。
- 没有在浏览器里用真实样本完成第 4 步。自选行为已在本基线的 WashU 契约里。目录页仍在 [PR #235](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/235)。外部登记的读取 API 尚未实现，本提交不登记真实文件。
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
derivatives = (root / "src/encode_pipeline/services/browser_derivatives.py").read_text(encoding="utf-8")
derivative_test = (root / "test/services/test_browser_derivatives.py").read_text(encoding="utf-8")
assert "## 3. PET、输出和 QC 单位" in contract
assert "| `Yield` | N/R |" in contract
assert "1−U/T" in contract and "K/C" in contract
registered = (root / "docs/architecture/registered-result-read-contract.md").read_text(encoding="utf-8")
ledger = (root / "docs/development/scientific-workflow-evidence-ledger.md").read_text(encoding="utf-8")
assert "hitrac-colleague-pet-qc-trial.md" in ops
assert "## 科学语义" in washu
assert "hitrac_bedpe_all" in washu and "hitrac_bedpe_no_bg" in washu
assert "## 在 WashU 内选择轨道文件（IRENE-21）" in washu
assert "在 WashU 内选择轨道文件（IRENE-21）" in prose
assert "不表示相应表、API 或客户端已经存在" in registered
assert "registered-result-read-contract.md" in prose
assert "## Hi-TrAC" in ledger
assert "scientific-workflow-evidence-ledger.md" in prose
assert "MAX_SOURCE_RECORDS = 2_000_000" in derivatives
assert "def test_source_record_cap_stops_before_sorting_extra_pairs" in derivative_test
assert '"sample_id": sample["display_id"]' in results
assert "source_artifact_id=document.source.artifact_id" in results
keys = (
    "raw_pairs", "trimmed_pairs", "mapping_ratio", "all_pets", "all_redundancy",
    "all_cis_ratio", "all_close_ratio", "all_middle_ratio", "all_distal_ratio",
    "no_bg_pets", "yield", "no_bg_cis_ratio", "no_bg_close_ratio",
    "no_bg_middle_ratio", "no_bg_distal_ratio",
)
for key in keys:
    assert f'"{key}"' in results
    assert f"`{key}`" in prose
for phrase in (
    "N/R", "1−U/T", "C/U", "K/C", "pair_support_missing", "rreg_", "rrev_",
    "experiment_snapshot", "generation_changed", "metadata.sample_id",
    "source_artifact_id", "source_too_large", "2,000,000",
    "已定位 PET，轨道未打开，试用未完成", "已排除同代另一行",
):
    assert phrase in prose, phrase
assert "只按产物类型和 generation 点选，会读到另一个样本的行" in prose
for banned in ("s000", "/home/", "11_17", "helixweave-hitrac-preprocess."):
    assert banned not in prose
print("hitrac-colleague-trial ok")
PY
```

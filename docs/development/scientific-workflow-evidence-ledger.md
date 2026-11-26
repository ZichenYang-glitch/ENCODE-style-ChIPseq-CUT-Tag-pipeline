# 科学工作流证据台账

本文件是 IRENE-35 的证据台账。它对照仓库里已经版本化的运行资格和 SCI-7 记录，分开写下 ENCODE、Bulk RNA-seq、Hi-TrAC 现在能证明的软件行为，以及仍然没有科学结论的默认值。同事可以按文末步骤核对结果和 QC。软件测试通过、资格文件生成、流程退出码 0 和 Gate 成功，都只说明软件按当时的身份跑完；它们本身不是生物学确认，也不冻结生产默认。

历史科学证据的阅读基线是 main `d85c3b6a50a9de0f108691bd60efc678d851a2be`。当时的 Bulk 资格文件和实现清单是历史证据身份。当前检出包含集成基线 main `65d99776a96830fd5a20092daf4a5c5b99461aeb` 的这两份文件，那是当前软件身份。两组摘要不同。历史公共数据记录和 Bulk Gate 没有在当前软件身份上重跑；不能把历史身份下的数字换成当前摘要后，声称那些记录用了新版本。下面的 SHA256 仍是文件内容摘要，文档提交之后可以用文末命令复核。本台账不改默认值，不下载数据，不启动流程。WashU 轨道选择由 IRENE-21 另行验收。11_17 的计算与核验跟踪在 IRENE-9，本单不打开、不重跑该批次。试用验收在 IRENE-34：由同事从已核验样本定位 PET 和 QC。

## 三种说法

| 说法 | 在本台账里指什么 |
| --- | --- |
| 软件行为 | 当前代码、契约和产品页面实际登记、计算或展示的内容 |
| 软件通过 | 测试、CI、Gate、退出码 0，或资格文件已经按闭包生成 |
| 科学结论 | 按预先登记的规则，用来源固定的真实文库决定某个默认能否冻结 |

## ENCODE

工作流 ID 是 `encode-style-chipseq-cuttag-atac-mnase`。科学政策修订号是 `SCIENTIFIC_POLICY_REVISION = 7`（`src/encode_pipeline/config/scientific_policy.py`）。修订号把一次运行的过滤、重复、细胞器和 QC 语义与旧输出目录隔开。它记录软件行为的版本，不表示 SCI-7 已经冻结默认。

### 已经在跑的默认

`config/config.yaml` 把 `duplicate_analysis` 写成 `auto`，`mapq` 写成 `30`。`auto` 在 `workflow/rules/metadata.smk` 里按样本 assay 解析：

```text
return "retain" if assay == "cuttag" else "exclude"
```

显式的 `exclude` 或 `retain` 作用于整次运行。Picard 只标记重复，分析 BAM 执行上述政策，MACS3 使用 `--keep-dup all`。配对末端的 Bowtie2 `-X` 在 `auto` 时取 `BOWTIE2_MAX_INSERT_BY_ASSAY`：ChIP-seq 2000、ATAC-seq 2000、CUT&Tag 700、MNase-seq 500。

QC 说明在 `docs/qc-interpretation.md`：QC 是描述性的，不是统一的通过/失败门。运行页的 QC 徽章如果出现，要回到产生该指标的工具去读，不能直接当成科学验收。

### 显式启用的 2/3 重复混合策略

`test/workflow/test_mixed_replicate_fixtures.py` 从公开校验、任务图和产物清单检验每个实验组的最终峰归属。合成窄峰的执行用例给三份独立生物重复输入同一位点、仅第二和第三份共享的位点、仅同一份重复中重叠的位点。默认阈值下前两个位点保留，最后一个剔除。设为三份支持时仅保留第一处。二重复 ChIP 组在同次运行由原 IDR 汇总规则生成 conservative 和 optimal 两份文件。

合成记录用于验证软件行为。历史严格模式的运行身份仍用修订 7，混合模式的身份另含显式策略。实验室数据上的表现需单独测量。[配置和产物说明](../reproducibility-policy.md#same-run-with-two--and-three-replicate-groups)。

### 版本化的科学记录

后续测量以这两份记录为准：

- `docs/development/sci7-public-data/README.md`
- `docs/development/sci7-public-data/sci7-batch-summary.md`

批次摘要的状态句是：`No production default is frozen.` 已测量的比较包括 CUT&Tag H3K27me3、组蛋白 H3K36me3、MNase Rep1、TF NR2F2、ATAC Treg，以及 JUND 确认对 ENCSR785RQR / ENCSR194GTV。摘要对每一项的结论都是 inconclusive、no default 或 no conclusion。退出码 0 在该文件里被明确排除在科学验收之外。

仍未测量的确认集：CUT&Tag、组蛋白、MNase、ATAC。MNase 预先登记的 500-arm 短模板指标没有评估。JUND 是另一个靶点和细胞背景，不是 NR2F2 的生物学重复。

`docs/development/sci7-public-data/reference-identity.md` 是 2026-10-28 的准备记录。它固定拟用参考的来源和官方压缩 MD5，并在自己的状态句里把大 FASTA、GTF 和 Bowtie2 索引写成当时尚未下载、尚未构建。`docs/development/sci7-public-data/README.md` 的参考段仍复述这句准备语句。后来的执行记录另行记下 FASTA 和索引的实测；GTF 和尚未执行的确认文库仍停在准备记录。

| 物种 | 2026-10-28 准备记录 | 后来的实测 | 仍缺 |
| --- | --- | --- | --- |
| 人 GRCh38，GENCODE v44 | 压缩 FASTA `GRCh38.primary_assembly.genome.fa.gz`，844,691,642 字节，官方 MD5 `9c3fc2ca260a767530dddb0f26721a6b`。该文件当时把载荷标为未核验 | `minimal-smoke-results.md`：四份 FASTQ 与这份 GRCh38 FASTA 的字节数、官方 MD5 一致，并由该 FASTA 构建 FAI、染色体长度表和六个 Bowtie2 索引。`tf-jund-confirmation.md` 在比对前再次核对同一压缩 FASTA，SHA256 `e9a2d5a5cd225293646ae298998ad4ea8c11e8f22729bc2d6f5c3dcfefbf8ef8` | GTF `gencode.v44.primary_assembly.annotation.gtf.gz`，官方 MD5 `b182a9f3b134b9cc2da566a2b1692557`。准备记录写明两份 GTF 都未下载；`minimal-smoke.md` 写明该次 smoke 不下载 GTF |
| 小鼠 GRCm38，GENCODE vM25 | 压缩 FASTA `GRCm38.primary_assembly.genome.fa.gz`，773,507,376 字节，官方 MD5 `3bc591be24b77f710b6ba5d41022fc5a`。该文件当时把载荷标为未核验 | `atac-treg-discovery.md`：同一文件，官方 MD5 同上，实测 SHA256 `617b10dc7ef90354c3b6af986e45d6d9621242b64ed3a94c9abeac3e45f18c17`。FAI、`chrom_sizes` 和六个 `.bt2` 由该解压 FASTA 用 samtools 1.23.1 与 Bowtie2 2.5.5 生成 | GTF `gencode.vM25.primary_assembly.annotation.gtf.gz`，官方 MD5 `c5125258a0a2c5250ddb4c192abbf4e8`，准备记录中未下载。确认文库 ENCSR213YKJ（组蛋白）和 ENCSR639PPN（ATAC）未下载、未评估 |

可选 GEO/SRA 的 CUT&Tag 与 MNase 确认集也未用于解释已完成的无结论记录，批次摘要没有把它们算作已评估。人的 JUND 确认文库已经测过，见 `tf-jund-confirmation.md`；那次测量不冻结默认。

`docs/development/sci7-public-data/resource-budget.md` 开头的 “No large FASTQ has been downloaded” 是 2026-10-28 的计划句；后来哪些文库已经测过，以批次摘要为准，不用预算文件的开头句覆盖它。预算文件里的约 1.2 TiB 也不是准入门槛。

Roadmap 对 #218 `0fef7a2…` 的记录是：workflow `43628729299` 成功，真实执行 169 passed、零失败、零跳过。同一段写明这不是生物学确认，也不是 Gate。`docs/development/workflow-platform-agent-roadmap.md` 把科学确认列在 A1 交付之外，并写明本计划不启动新数据批次。

### 仍缺科学证据的 ENCODE 默认

这些值已经写进软件，SCI-7 没有冻结它们：

| 软件默认 | 已有测量 | 科学结论 |
| --- | --- | --- |
| ChIP-seq / ATAC-seq / MNase-seq `auto` → exclude | 各有一批，触发条件未满足 | 不冻结 exclude |
| CUT&Tag `auto` → retain | Rep1 浅、Rep2 未建立为独立生物学重复 | 不冻结 retain |
| ChIP-seq / ATAC-seq `max_insert` 2000 | NR2F2、ATAC、JUND 均无结论 | 不冻结 2000 |
| CUT&Tag `max_insert` 700 | retain 对 exclude、500 对 700 均无结论 | 不冻结 700 |
| MNase-seq `max_insert` 500 | Rep1 无结论；登记的短模板指标未评估 | 不冻结 500 |
| `mapq` 30 | 不在上述单变量问题里 | 没有一份确认记录冻结 30 |

### 在已完成的 ENCODE 运行上核对

使用一次已经成功、且你有权阅读的运行。本步骤不创建运行。

1. 打开 `/runs/<run_id>`。标题下的链接应是 `encode-style-chipseq-cuttag-atac-mnase`。
2. 有参考绑定时，定义列表 `run-reference-profile` 显示 Reference、Organism / assembly、Revision、Identity digest。没有这块时，只说明这次响应没有 `reference_profile`，不能据此推断运行没用参考。
3. 打开 QC 标签，或直接访问 `/runs/<run_id>?view=qc`。成功之后才会出现已索引指标；每行有显示名和 `metric_key`。这些数字是描述性 QC。运行状态成功，不把 `duplicate_analysis: auto` 或 `max_insert: auto` 变成已冻结默认。
4. 核对参考时按物种读后续实测。`reference-identity.md` 的状态句和 README 参考段仍把大 FASTA 写成未构建，那是准备记录，不用它覆盖后来的实测。人 GRCh38 看 `minimal-smoke-results.md` 的官方字节数、MD5 和六个索引，以及 `tf-jund-confirmation.md` 比对前再次核对的压缩 FASTA：844,691,642 字节，官方 MD5 `9c3fc2ca260a767530dddb0f26721a6b`，SHA256 `e9a2d5a5cd225293646ae298998ad4ea8c11e8f22729bc2d6f5c3dcfefbf8ef8`。小鼠 GRCm38 看 `atac-treg-discovery.md`：官方 MD5 `3bc591be24b77f710b6ba5d41022fc5a`，SHA256 `617b10dc7ef90354c3b6af986e45d6d9621242b64ed3a94c9abeac3e45f18c17`，以及由此生成的 FAI、染色体长度表和六个 Bowtie2 索引。两份 GTF 仍只有准备记录中的官方 MD5。页面上的 Identity digest 是这一次运行修订的身份；它与上述 SHA256 相同，只说明用了同一压缩 FASTA，不把默认变成科学结论。组蛋白确认 ENCSR213YKJ、ATAC 确认 ENCSR639PPN，以及可选 CUT&Tag/MNase 确认文库，仍按批次摘要视为未评估。

## Bulk RNA-seq

工作流 ID 是 `bulk-rnaseq`。上游是 nf-core/rnaseq 3.26.0，提交 `e7ca46272c8f9d5ceee3f71759f4ba551d3217a4`，路线 `star_salmon`。结果契约是 `src/encode_pipeline/contracts/nfcore_rnaseq/results-contract-3.26.0.json`，文件 SHA256 `a099f495de027a385580cf14e2316b7e96ea6d67b971af02f3180cbf751e859f`。

### 版本化资格

资格文件是 `src/encode_pipeline/contracts/nfcore_rnaseq/default-execution-qualification-1.1.0.json`，实现清单是 `src/encode_pipeline/contracts/nfcore_rnaseq/execution-implementation-manifest-1.0.0.json`。`candidate_status` 仍是 `source-enabled`，闭包仍是 119 个文件，运行时 `identity_sha256` 仍是 `4c2f9c99bfe07141db39c0974d652bf006a5d596351e11caacf426da62699f9e`。当前软件身份与历史证据身份不同，是因为实现清单里若干已有源码文件换了字节；路径数没有增加，结果契约也没有变。这是软件闭包变化，不是 SCI-7 或 Bulk Gate 的重跑。

历史证据身份，阅读基线 `d85c3b6a50a9de0f108691bd60efc678d851a2be`：

historical-qualification-file: 5f050c81b28cefe8bacd6b499939f46210d166b4f4fd8630a593017a4637603b
historical-record-sha256: aa0ca90e70052d8cc336a06f9028963c00f1f516a6379333a94650117e009202
historical-aggregate-sha256: 68a3798931941bd8e99d8eec34a770ad4d453faa30e4f0b7017db93840566fc2
historical-manifest-file: b3862eeec44cedc890ae11e61297e99452f494a095c1bbc656dd9c1933d20a89

当前软件身份，本检出中的文件字节。main `65d99776a96830fd5a20092daf4a5c5b99461aeb` 记录的是当时的资格文件和实现清单；目录页、WashU 打包和长期盘迁入同一棵树后，119 个闭包文件里已有源码换了字节，所以上面的历史四行不变，下面四行按本检出重算。这仍是软件闭包变化，不是 SCI-7 或 Bulk Gate 的重跑：

current-qualification-file: 1d7e8deffd1e89eeeab1d9eda213438ff642560dd7c1fcde06732b9a5598b4db
current-record-sha256: 53e90ae4011d9c548ab7f28ac0f7181d39d64305152ceb9d48ee66d8a0d64b75
current-aggregate-sha256: 7aeeeba423ddd42800f600c7900ae1e141c4cb895f9b47a0d342ef064cbb62a0
current-manifest-file: d15aad702ede68f779beb8621312adf023596de36d487dbdaf7d13afb9e84b68

历史科学证据没有在当前软件身份上重跑。

`candidate_status: source-enabled` 表示这份源码候选已绑定闭包、持久化契约、命令、参考协议和结果契约。`docs/architecture/bulk-rnaseq-adapter-decision.md` 写明：资格字节匹配的是实时核对过的实现；PR #152 的产物发现和 QC 提取、PR #153 的合成 `rapid_quant` 资格以及本地 STAR+Salmon 验收，证明的是固定离线执行和平台机制。`rapid_quant` 是服务器持有的资格模式，不是产品默认分析。该决定同时写明这些夹具不构成生物学有效性。

`docs/development/workflow-platform-agent-roadmap.md` 在 2026-11-16 的入口里记录：Bulk Gate 跑在 `e112223…`，GitHub Actions run `43792685991` attempt 2，job `133061548271` 通过 4 项执行测试和 5 项产品浏览器测试。那是该提交上的软件 Gate。同一入口当时的 main 是 `1b1db421ba3172bc39cf31ce0d03976511379fda`；roadmap 写明不能把这次 Gate 改记成在 merge SHA 上重新跑过。该次合入 push CI `43816767856` 里，`bulk-rnaseq-real-execution` 为 skipped，skipped 不是执行通过。历史证据身份读自 `d85c3b6a50a9de0f108691bd60efc678d851a2be`，Gate 没有在这个 SHA 上重跑。main `65d99776a96830fd5a20092daf4a5c5b99461aeb` 上的资格文件是更早的当前软件身份；本检出的当前四行摘要已经换过，Gate 也没有在这组摘要上重跑。Gate 成功不改写 `source-enabled`。SCI-7 记录和这份 Bulk 资格都没有把 nf-core 默认冻结成生物学结论。

### 同事如何定位定量、QC 和参考

仍使用一次已完成、且你有权阅读的 `bulk-rnaseq` 运行。不新建运行。

参考在运行页同一块 `run-reference-profile`：Reference、Organism / assembly、Revision、Identity digest。资格文件里的参考协议只规定绑定如何哈希，不指定这一次的物种或组装。组装以页面上这一次修订为准。

定量在 Artifacts 标签，或 `/artifacts` 的 Artifact publication filters。筛选字段是 Workflow、Run、Artifact type。Artifact type 一次填一个完整 `output_type`：

| 要找的结果 | Artifact type |
| --- | --- |
| 样本转录本定量 | `bulk_rnaseq.salmon.quant_transcript` |
| 样本基因定量 | `bulk_rnaseq.salmon.quant_gene` |
| Salmon 样本元信息 | `bulk_rnaseq.salmon.meta_info` |
| 合并基因计数 | `bulk_rnaseq.salmon.merged.gene_counts` |
| 合并基因长度缩放计数 | `bulk_rnaseq.salmon.merged.gene_counts_length_scaled` |
| 合并基因缩放计数 | `bulk_rnaseq.salmon.merged.gene_counts_scaled` |
| 合并基因长度 | `bulk_rnaseq.salmon.merged.gene_lengths` |
| 合并基因 TPM | `bulk_rnaseq.salmon.merged.gene_tpm` |
| 合并转录本计数 | `bulk_rnaseq.salmon.merged.transcript_counts` |
| 合并转录本长度 | `bulk_rnaseq.salmon.merged.transcript_lengths` |
| 合并转录本 TPM | `bulk_rnaseq.salmon.merged.transcript_tpm` |
| tx2gene | `bulk_rnaseq.salmon.merged.tx2gene` |
| 增补 tx2gene | `bulk_rnaseq.salmon.merged.tx2gene_augmented` |

结果契约里的 `bulk_rnaseq.salmon.merged.*` 是家族写法。页面筛选要用上一表的完整类型。样本定量文件对应上游 `star_salmon/<sample>/quant.sf` 与 `quant.genes.sf`；合并表对应 `star_salmon/salmon.merged.<suffix>.tsv`。能打开这些类型，说明产品把定量登记成了可筛选产物。它不说明这次数值已经过生物学验收。

QC 在 `/runs/<run_id>?view=qc`。成功完成后列出已索引行，每行有显示名和 `metric_key`，并可打开来源产物。Salmon 定量旁边的机器可读 QC 至少包括：

- `salmon.processed_fragments`
- `salmon.mapped_fragments`
- `salmon.mapping_fraction`

比对与文库 QC 使用同一页上的 `star.*`、`fastqc.*`、`featurecounts.*`、`rseqc.*`。部分 FastQC 行的徽章把 FastQC 自己的 PASS、WARN、FAIL 映射成 pass、warning、fail。那是 FastQC 模块状态，不是 HelixWeave 的科学结论。没有徽章时页面写 Not reported。

### 仍缺科学证据的 Bulk 默认

上游 3.26.0 的 STAR+Salmon 参数、比对过滤和定量列，由结果契约和资格闭包固定为软件行为。本仓库没有一份来源固定的真实队列记录，用来决定这些默认是否回答某个实验室问题。因此当前可以确认的是：用户能按上表找到定量、QC 和参考身份。不能确认的是：某一次 Bulk 结果在生物学上成立，或 nf-core 默认已经冻结为本实验室的科学政策。

## Hi-TrAC

工作流 ID 是 `hitrac-preprocess`，名称是 Hi-TrAC preprocessing。上游身份是 cLoops2/tracPre2，Git 提交 `de6cc732fa00b408551b9f4272933640c08447f1`。准入要求的脚本 SHA256 是 `c3c4ef4e6287fa4ade97a6f5980d20c81ea345b8e8e13c88ae3b7c4bd3a67aec`，与 `config/hitrac_preprocess/tools.lock.json` 的 `upstream.script_sha256` 相同。锁文件 SHA256 是 `392e4a4ff93f2b14b380f16572d74a77dbca699112800bcafe50992d692090da`，共 124 个 conda 包。

### 预处理和下游分析是两件事

平台实际启动的是这份锁定的 `tracPre2.py`。`qualification.py` 的参数是 `-fqd`、`-o`、`-ref`、`-n 1`、`-p`、`-mapq`。调用计划里与 cLoops2 命令相关的阶段是 `cLoops2 qc`，用来生成预处理 QC，不是 loop、domain 或 peak 调用。

公开产物只有三类：`hitrac_bedpe_all`、`hitrac_bedpe_no_bg`、`hitrac_summary`。QC 只从 `hitrac_summary` 提取，文件名是 `tracPre_summary.txt`，相对路径形如 `results/hitrac/<64 hex>/tracPre_summary.txt`。十五个指标键是 `raw_pairs`、`trimmed_pairs`、`mapping_ratio`、`all_pets`、`all_redundancy`、`all_cis_ratio`、`all_close_ratio`、`all_middle_ratio`、`all_distal_ratio`、`no_bg_pets`、`yield`、`no_bg_cis_ratio`、`no_bg_close_ratio`、`no_bg_middle_ratio`、`no_bg_distal_ratio`。QC 行的 assay 字段是字符串 `Hi-TrAC`；区分预处理和下游时，以运行页上的工作流 ID `hitrac-preprocess` 和上述三个产物类型为准。

作者可填的选项只有 `threads`（整数 2–8，默认 2）和 `mapq`（整数 0–255，默认 10）。配置对象不接受科学覆盖。没有 loop、domain、peak 或差异分析字段。

运行时锁的 `package_sources` 含有 `callCisLoops.py`、`callTransLoops.py`、`callDiffLoops.py`、`callDomains.py`、`callPeaks.py` 的源码摘要。那是整包 cLoops2 0.0.5 源码闭包的一部分。这些文件出现在锁里，只说明包字节被固定；它不表示平台执行了这些入口。

### 软件资格和产品状态

没有 `HELIXWEAVE_HITRAC_RUNTIME_BINDING` 与 `HELIXWEAVE_HITRAC_RUNTIME_SHA256`，或二者核验失败时，`load_default_hitrac_adapter` 注册只写作的适配器，`execution_availability` 返回 `not_configured`。两处坐标都通过准入后，结果适配器才把执行报成 `available`，并发布上面三类产物。页面上的 `not_configured` 表示服务器还没有这份准入，不是某次实验的科学失败。

`docs/development/hitrac-preprocess-contract.md` 记录了合成 tiny 基线：设计参考上 all=8、noBg=5，以及随后的平台链。那些原始日志写在 `/tmp/helix-hitrac-*`，不在本 Git 树里。本台账没有重跑 H5，也不能从仓库重新打开那些日志。契约里的 all=8、noBg=5 是当时对合成基线的记录。

`docs/development/adapter-conformance.md` 仍写着 Hi-TrAC 本地执行是 `not_configured` until qualified artifact/QC publication is implemented。结果发布代码已经在 `results.py`。同事判断一次运行能否执行时，以该环境的准入和环境变量为准，不用这句表格代替代码。

### 在已完成的 Hi-TrAC 运行上核对

1. 运行页链接是 `hitrac-preprocess`。Reference 块同样给出这次参考修订的名称、物种、组装、修订号和 Identity digest。
2. Artifacts 里应能按三个类型找到 BEDPE 和 `tracPre_summary.txt`。出现 loop、domain 或 peak 文件时，它不属于这条预处理工作流的公开契约，应单独追溯来源。
3. QC 标签的 `metric_key` 应落在上面十五个键里，来源产物类型是 `hitrac_summary`。`all_pets` 是 MAPQ 过滤后的 all PET 数；`no_bg_pets` 是 noBg 集合。两个数不是同一集合，不能把 all 的冗余率当成 unique/noBg 文件的计数。
4. 政策 A 拒绝空 all、空 noBg 或缺少 cis 分母的成功发布。这是接入时的成功语义。契约写明，不能据此说该实验没有科学意义。

### 仍缺科学证据的 Hi-TrAC 默认

`mapq` 默认 10 是原脚本的筛选值，本仓库没有一份真实 Hi-TrAC 文库记录把它冻结成实验室阈值。`threads` 默认 2 是资格入口的线程数，不是生物学最适值。H5 使用的是种子 `20260924` 的人工 100 kb 参考，不是基因组尺度的科学结果。本仓库没有注册第二条 Hi-TrAC 下游工作流，因此预处理通过不表示 loop、domain 或 peak 分析已经完成或已经确认。

## 不启动新的真实计算

本单没有新的科学问题，也没有新的资源预算。SCI-7 里尚未下载的确认集保持未测：组蛋白确认是不同物种和不同标记，ATAC 确认与发现集共享供体，可选的 CUT&Tag / MNase 确认没有被用来改阈值。预算文件不授权按 1.2 TiB 情景扩跑。Hi-TrAC 下游分析和 Bulk 生物学队列同样没有在本单获得问题和预算。因此这里不提议、也不启动任何新的真实文库计算。

## 未确认项

- ENCODE 的重复政策和 `max_insert` 候选仍是软件默认；SCI-7 结论保持 no conclusion。
- ENCODE `mapq` 30 没有对应的冻结记录。
- 人 GRCh38 与小鼠 GRCm38 的压缩 FASTA，以及由它们生成的 FAI、染色体长度表和六个 Bowtie2 索引，已有 `minimal-smoke-results.md`、`tf-jund-confirmation.md` 和 `atac-treg-discovery.md` 的实测。两份 GTF 仍未下载。确认文库 ENCSR213YKJ、ENCSR639PPN 和可选 CUT&Tag/MNase 确认集仍未评估。`reference-identity.md` 的状态句和 README 参考段仍把大 FASTA、GTF 和索引写成未构建，那是 2026-10-28 的准备语句。
- Bulk 的 `source-enabled` 资格和 `e112223…` Gate 是软件身份与软件测试。历史证据身份与当前软件身份不同，历史记录没有用当前实现清单重跑。没有 Bulk 公共数据科学确认记录。
- Hi-TrAC 的 H5 原始日志不在 Git 中；本单没有复测 all=8、noBg=5。
- `adapter-conformance.md` 关于 “publication is implemented” 的句子落后于 `results.py`。
- `resource-budget.md` 开头的下载状态句落后于 SCI-7 批次摘要。
- 11_17 的计算与核验跟踪在 IRENE-9，试用验收在 IRENE-34。这两项都不在本单范围。
- WashU 选择、加载和移除不在本单范围。

## 验收

在仓库根目录执行：

```bash
python3 scripts/check_scientific_workflow_evidence.py
```

命令只读本仓库文件，退出码 0 表示本台账引用的内容摘要、资格字段和关键源码句子仍与文件一致，并且历史证据身份与当前软件身份分开核对：当前四行摘要必须等于本检出的资格文件和实现清单，历史四行摘要必须仍是阅读基线 `d85c3b6a50a9de0f108691bd60efc678d851a2be` 的那一组。它不连接平台，不证明任何一次真实运行的生物学结论。对真实运行，按上面三节的页面步骤，用已有且已授权的运行 ID 核对；不要为了核对而创建或重跑批次。

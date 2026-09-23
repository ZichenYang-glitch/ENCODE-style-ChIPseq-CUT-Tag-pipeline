# 上游耦合与升级核对台账

状态：PR-5 文档返修已验收；2026-09-23 同步 fastp 字段修复（待独立复验）。
其余条目保持 2026-09-22 的核对范围。
本台账记录已有消费者的依赖和升级义务，不批准修改科学口径，不替代缺陷审计。
设计背景沿用 [bulk adapter 决策](../architecture/bulk-rnaseq-adapter-decision.md)，
交付顺序沿用 [roadmap](workflow-platform-agent-roadmap.md#已确认的维护顺序与审核边界2026-09-22)。
没有重复建立版本锁或输出契约；下文是它们的维护索引。

## 版本与证据边界

- **固定版本**：bulk 使用 nf-core/rnaseq **3.26.0**，commit
  `e7ca46272c8f9d5ceee3f71759f4ba551d3217a4`（下文 NF）。已读取对应
  [schema][NF-schema]、[入口配置][NF-config]、[主工作流][NF-main]和实际模块源码；
  下载的 NF 文件与本地 [source manifest](../../src/encode_pipeline/contracts/nfcore_rnaseq/source-manifest-3.26.0.json)
  中同路径的大小、SHA-256 核对。下文 bulk 条目均限定本地支持的 `star_salmon` 路径。
- **工具声明与运行身份**：[results contract](../../src/encode_pipeline/contracts/nfcore_rnaseq/results-contract-3.26.0.json)
  记录工具版本/格式；[container inventory](../../src/encode_pipeline/contracts/nfcore_rnaseq/container-inventory-3.26.0.json)
  记录上游镜像坐标和模块源码身份。镜像 tag 不等于 OCI digest；本轮没有读取部署侧
  availability lock、加载镜像或查询本机工具版本，**实际执行版本和部署 OCI digest 均未核实**。
  inventory 所列的可得坐标见文末，不能用它们冒充实际容器身份。
- MultiQC **1.33** tag 对应 commit `5953b5417ccb70bf4a2309562d43015fced8b585`；
  ENCODE 的 **1.35** tag 对应 `87e504bb77e3d94687dd980638724b7ef8e1e92f`。
  本轮通过官方 GitHub tag API 核对，文中链接使用 commit。
  bulk 的 1.33 defaults/base_module/table_object 字节分别与 results contract 的记录吻合。
- **源码核对（S）**：已读本地调用者及相应官方版本源码，没有据“存在函数”推断完整流程执行成功。
  **合成契约测试（T）**：初次交付实际通过 canonical bootstrap 执行的已有测试，输入为内存构造的
  schema、ZIP、JSON、文本或 TSV；这不等于运行上游科学工具。
  下文列出的“已有测试”还有未运行的项；初次执行清单与结果保留在原证据目录。
  PR-5 文档返修未重跑行为测试；UC-06 保留该次独立审核的失败对照，
  并单列 2026-09-23 修复后的合成契约验证，不能与真实工具运行混称。
- PR-5 不改行为、测试、依赖或身份，不重建 manifest/qualification，不运行 Bulk Gate。
  PR-3、PR-4 和 G0 的 Gate 欠验分别保留；本台账不能补足任何一项欠验。

## 本地接线与命令约定

基类 [adapter.py:70](../../src/encode_pipeline/adapters/bulk_rnaseq/adapter.py#L70) 的 `BulkRnaSeqWorkflowAdapter`
承担校验、规划及命令相关委托；其 [extract_artifacts:254](../../src/encode_pipeline/adapters/bulk_rnaseq/adapter.py#L254)
返回 unsupported，不提供这里描述的结果提取能力。
结果子类 [adapter.py:275](../../src/encode_pipeline/adapters/bulk_rnaseq/adapter.py#L275) 的 `BulkRnaSeqResultsWorkflowAdapter`
提供产物发现及 QC 提取：[extract_artifacts:293](../../src/encode_pipeline/adapters/bulk_rnaseq/adapter.py#L293)
委托 [artifacts.py:269](../../src/encode_pipeline/adapters/bulk_rnaseq/artifacts.py#L269) 的 `discover_bulk_rnaseq_artifacts`，
[extract_qc_metrics:305](../../src/encode_pipeline/adapters/bulk_rnaseq/adapter.py#L305)
委托 [qc.py:163](../../src/encode_pipeline/adapters/bulk_rnaseq/qc.py#L163) 的 `extract_bulk_rnaseq_qc_metrics`。
后者按封闭的 `BULK_RNASEQ_QC_SOURCE_TYPES` 分派，而不是任意解析报告目录。
机器表的公开下载审计与数字 QC 提取是两个消费者；表被允许发布，不代表每个字段都成为平台 QC 指标。

下列升级命令中的“目标”接在这一公共前缀后执行；每次使用新临时目录：

```bash
CHECK_DIR=$(mktemp -d /tmp/helix-upstream-check.XXXXXX)
TMPDIR="$CHECK_DIR" XDG_CACHE_HOME="$CHECK_DIR/cache" PYTHONDONTWRITEBYTECODE=1 \
  ./.local/envs/ci-fast/bin/python -I -S scripts/checkout_bootstrap.py \
  --repository-root . pytest <目标及可选 -k 表达式> -q --basetemp "$CHECK_DIR/pytest"
```

这是升级时的核对说明，不授权安装工具或运行真实工作流；真实工具测试另需按
[real-execution harness](real-execution-harness.md) 准备锁定工具，并重新检查 fixture、配置加载和写入位置。

## UC-01 — bulk MultiQC 样本清洗、替换与 mate 所有权

**上游**：NF + MultiQC 1.33。[NF replacement/grouping 接线][NF-multiqc]、
[NF 清洗配置][NF-multiqc-config]、[MultiQC 清洗函数][MQ-base]及 [defaults][MQ-defaults]。
`multiqcNameReplacements` 从 FASTQ `simpleName` 生成全局替换表；R1 的 simpleName 已等于样本 ID 时，
该行的 R1/R2 替换都被省略。清洗后才应用替换；不能仅比较文件 basename 或直接假定所有 `_1/_2` 都是 mate。

**本地依赖**：[authoring.py:46](../../src/encode_pipeline/adapters/bulk_rnaseq/authoring.py#L46) 的 `MULTIQC_SAMPLE_CLEAN_TOKENS` 的 159 个清洗字面量保留策略；
[validation.py:604](../../src/encode_pipeline/adapters/bulk_rnaseq/validation.py#L604) 的 `_validate_normalized_result_identities` 在有效 MultiQC 路径调用
[validation.py:635](../../src/encode_pipeline/adapters/bulk_rnaseq/validation.py#L635) 的 `_validate_multiqc_sample_identities`，并使用
[validation.py:787](../../src/encode_pipeline/adapters/bulk_rnaseq/validation.py#L787) 的 `_nextflow_simple_name` / [validation.py:794](../../src/encode_pipeline/adapters/bulk_rnaseq/validation.py#L794) 的 `_multiqc_clean_identity`。
所有权同时区分 canonical sample 和 mate，纳入 authored PE 与 UMI 丢弃一端后的布局。
当前次序为 NF extra truncate/remove →默认 clean → trim →全局 exact replacement；
不同所有者重合即拒绝，同一所有者的多 lane 可接受。此处已有保护，不能把潜在的非单射清洗直接写成未修 bug。

**测试**：[test_bulk_rnaseq_adapter.py:1201](../../test/adapters/test_bulk_rnaseq_adapter.py#L1201) 的 `test_multiqc_identity_collision_is_rejected_but_safe_underscores_remain_valid`；
[test_bulk_rnaseq_adapter.py:1411](../../test/adapters/test_bulk_rnaseq_adapter.py#L1411) 的 `test_multiqc_cleaning_before_exact_replacement_cannot_change_owner`；
[test_bulk_rnaseq_adapter.py:1495](../../test/adapters/test_bulk_rnaseq_adapter.py#L1495) 的 `test_multiqc_fastq_simplename_cannot_map_both_mates_from_one_exact_key`。
另有 `test_sample_ids_reserve_every_pinned_multiqc_cleanup_literal` 覆盖保留表，
`test_multiqc_identity_graph_deduplicates_lanes_but_not_biological_samples` 覆盖 lane。
这些测试核对本地拒绝/接受，不会启动 MultiQC 来穷举等价性。

**升级复核**：逐项 diff defaults、NF extra 清洗规则、替换先后次序、R1 no-op 条件与 UMI 布局；
检查 global replacement 是否仍覆盖所有来源，而非只检查 canonical ID 唯一。
目标：`test/adapters/test_bulk_rnaseq_adapter.py -k multiqc`。
证据：S + T；没有实际 MultiQC 1.33 运行。

## UC-02 — bulk 公开 MultiQC 机器表与 General Stats 分组

**上游**：NF + MultiQC 1.33；[NF 的 table_sample_merge 构造][NF-multiqc]、
[MultiQC 分组][MQ-base]及 [table object][MQ-table]。
General Stats 可使用 `S Read 1` / `S Read 2`，不能用解析前的 `S_1` / `S_2` 猜测它们。

**本地依赖**：[artifacts.py:1973](../../src/encode_pipeline/adapters/bulk_rnaseq/artifacts.py#L1973) 的 `_add_multiqc_specs` 只登记指定机器表；
[artifacts.py:655](../../src/encode_pipeline/adapters/bulk_rnaseq/artifacts.py#L655) 的 `_audit_published_multiqc_tables` → [artifacts.py:741](../../src/encode_pipeline/adapters/bulk_rnaseq/artifacts.py#L741) 的 `_multiqc_row_policy` →
[artifacts.py:810](../../src/encode_pipeline/adapters/bulk_rnaseq/artifacts.py#L810) 的 `_general_stats_owners` / [artifacts.py:867](../../src/encode_pipeline/adapters/bulk_rnaseq/artifacts.py#L867) 的 `_validate_multiqc_sample_rows`。
表首列严格为 `Sample`，UTF-8 TSV，禁止 BOM/NUL/CR、空行、重复行、列数不匹配和外来样本；
不同表使用 authored/effective/阈值后样本集合，General Stats 的允许分组集合单独生成。
部分表允许缺席（例如 alignment 路径可缺 MultiQC Salmon 表），不能统一补零或强制所有表出现。
HTML、`multiqc_data.json`、sources/log/parquet 不在该公开机器表清单中，不能推断平台可下载它们。

**测试**：[test_bulk_rnaseq_artifacts.py:1497](../../test/adapters/test_bulk_rnaseq_artifacts.py#L1497) 的 `test_general_stats_uses_exact_fixed_multiqc_grouped_identities`；
[test_bulk_rnaseq_artifacts.py:1408](../../test/adapters/test_bulk_rnaseq_artifacts.py#L1408) 的 `test_multiqc_table_duplicate_and_malformed_rows_fail_closed`；
已有 `test_every_published_multiqc_table_policy_rejects_a_foreign_row`、
`test_star_salmon_accepts_absent_optional_multiqc_salmon_table` 核对各表所有权和缺席策略。
这些测试不验证完整 HTML 中所有指标的可见性。

**升级复核**：核对输出目录 `multiqc/star_salmon/multiqc_report_data`、表名、首列、
模块 anchor 与 `Read 1/2` 分组形状，并逐表复核 optional/exact rows/required owners。
目标：`test/adapters/test_bulk_rnaseq_artifacts.py -k 'multiqc or general_stats'`。
证据：S + T；未运行上游报告生成或浏览器。

## UC-03 — nf-core 原生参数白名单、转换及提交入口

**上游**：NF 的 [nextflow_schema.json][NF-schema] 与 [samplesheet schema][NF-samples]。
vendored 参数 schema：60,620 bytes，SHA-256
`8f2f84a25c0aec65a18234cf01acdd74f2385e8dfac8417e4bad23a70bfb4388`；
samplesheet schema：3,218 bytes，SHA-256
`013b669b1a3d38709f548a2548350ecdabb3f2ba578f2b7de843105e2ce87a7d`。
均与本轮官方下载字节相同；schema 内 `$id` 的 master URL 不改变本地按固定哈希加载的事实。

**本地依赖**：[upstream.py:336](../../src/encode_pipeline/adapters/bulk_rnaseq/upstream.py#L336) 的 `upstream_parameter_properties` 要求 133 个上游键；
五个互斥集合合计 133：Standard 43、Advanced 16、platform-owned 57、raw args 7、unsupported 10。
[validation.py:297](../../src/encode_pipeline/adapters/bulk_rnaseq/validation.py#L297) 的 `_validate_advanced` 先区分危险运行键、未知键、Standard 冲突、平台所有权、raw args 和 unsupported；
允许值再经原 schema、公开 schema 和 [validation.py:365](../../src/encode_pipeline/adapters/bulk_rnaseq/validation.py#L365) 的 `_validate_advanced_semantics` 检查。
UI 的 [upstream.py:352](../../src/encode_pipeline/adapters/bulk_rnaseq/upstream.py#L352) 的 `projected_advanced_properties` 去掉上游 default/help 等注解，不把 NF 默认值自动变为作者输入。

Advanced 的 16 键为：`bam_csi_index`、`deseq2_vst`、`featurecounts_feature_type`、
`featurecounts_group_type`、`gffread_transcript_fasta`、`gtf_extra_attributes`、
`gtf_group_features`、`min_mapped_reads`、`min_trimmed_reads`、`rseqc_modules`、
`skip_gtf_filter`、`skip_gtf_transcript_filter`、`star_ignore_sjdbgtf`、
`stranded_threshold`、`stringtie_ignore_gtf`、`unstranded_threshold`。
其余完整集合以 [upstream.py:36](../../src/encode_pipeline/adapters/bulk_rnaseq/upstream.py#L36) 的 `STANDARD_NATIVE_PARAMETERS`、
[upstream.py:109](../../src/encode_pipeline/adapters/bulk_rnaseq/upstream.py#L109) 的 `PLATFORM_OWNED_NATIVE_PARAMETERS`、[upstream.py:173](../../src/encode_pipeline/adapters/bulk_rnaseq/upstream.py#L173) 的 `RAW_ARGUMENT_NATIVE_PARAMETERS` 为准；
例如 `aligner` 是 Standard 冲突，`input`/`outdir` 属平台，`extra_star_align_args` 属 raw args。
10 个已知但不支持的键是 `bracken_precision`、`contaminant_screening`、
`contaminant_screening_input`、`kallisto_quant_fraglen`、`kallisto_quant_fraglen_sd`、
`pseudo_aligner_kmer_size`、`save_bbsplit_reads`、`save_kraken_assignments`、
`save_kraken_unassigned`、`use_rustqc`。额外危险执行键如 `-c`/`profile`/`workDir` 也被拒绝，
不能只在上述 133 键中寻找运行时逃逸开关。未暴露上游能力是本地产品限制，不自动成为缺陷。

[validation.py:446](../../src/encode_pipeline/adapters/bulk_rnaseq/validation.py#L446) 的 `_normalize_inputs` 将 QC 主/子开关转为 `skip_*`，组合 `star_salmon`，
按 enabled 条件生成 UMI/rRNA 参数，并合并已验证 Advanced。
[execution.py:778](../../src/encode_pipeline/adapters/bulk_rnaseq/execution.py#L778) 的 `_runtime_params` 再注入受控 `input`/`outdir`、`validate_params=true` 等；
[execution.py:811](../../src/encode_pipeline/adapters/bulk_rnaseq/execution.py#L811) 的 `_samplesheet_bytes` 按 sample/library/lane 排序，实际 CSV 五列为
`sample,fastq_1,fastq_2,strandedness,seq_platform`，library/lane 本身不输出为列。
[execution.py:329](../../src/encode_pipeline/adapters/bulk_rnaseq/execution.py#L329) 的 `build_bulk_rnaseq_command` 以固定源码和 `-params-file` 读取 JSON，不把 Advanced 拼为任意 CLI。
函数存在/规范化成功不代表完成 Nextflow 执行。

**测试**：[test_bulk_rnaseq_adapter.py:385](../../test/adapters/test_bulk_rnaseq_adapter.py#L385) 的 `test_upstream_schema_is_exact_and_parameter_policy_is_closed`、
[test_bulk_rnaseq_adapter.py:735](../../test/adapters/test_bulk_rnaseq_adapter.py#L735) 的 `test_valid_advanced_allowlist_is_exactly_validated_and_classified`、
[test_bulk_rnaseq_adapter.py:825](../../test/adapters/test_bulk_rnaseq_adapter.py#L825) 的 `test_advanced_conflicts_and_dangerous_parameters_fail_closed`，另有
`test_advanced_wrong_types_and_stricter_semantics_fail_closed`。
本轮没有调用构建/物化运行时 fixture；CSV/argv 接线为源码核对。

**升级复核**：比较新 schema 的增删键、类型/范围/default 与五类分区，逐一评估已接受参数对文件名、样本集和路由的影响；
复核 native skip 的意义及 CSV 字段顺序，不能直接扩大 allowlist。
目标：`test/adapters/test_bulk_rnaseq_adapter.py -k 'upstream_schema or advanced or normalize'`。
证据：S + T；未运行上游参数校验器/Nextflow。

## UC-04 — FastQC ZIP 内部结构、文件名与整数计数

**上游**：FastQC **0.12.1** [BasicStats 输出][FQC-basic]，以及 NF [FASTQC 模块][NF-fastqc]。
**本地**：[qc.py:680](../../src/encode_pipeline/adapters/bulk_rnaseq/qc.py#L680) 的 `_parse_fastqc`、[qc.py:854](../../src/encode_pipeline/adapters/bulk_rnaseq/qc.py#L854) 的 `_fastqc_identity`、[qc.py:951](../../src/encode_pipeline/adapters/bulk_rnaseq/qc.py#L951) 的 `_fastqc_data`、
[status_evidence.py:478](../../src/encode_pipeline/adapters/bulk_rnaseq/status_evidence.py#L478) 的 `parse_fastqc_total_sequences`。
依赖受控 ZIP 根目录、`summary.txt` / `fastqc_data.txt`，版本标题、Filename、模块状态和
`Total Sequences` 整数；raw/trimmed/filtered 与 single/read1/read2 身份必须和期望路径对应。
PE post-trim 两端 retained reads 必须相等。禁止重复成员、逃逸路径、超限解压及 summary/data 状态冲突。
`Total Bases` 的简写不是精确最终保留碱基数，不能据此补一个 exact metric。

**测试**：[test_bulk_rnaseq_qc.py:443](../../test/adapters/test_bulk_rnaseq_qc.py#L443) 的 `test_fastqc_zip_extracts_decimal_metrics_flags_and_source_binding`、
[test_bulk_rnaseq_qc.py:517](../../test/adapters/test_bulk_rnaseq_qc.py#L517) 的 `test_fastqc_identity_and_closed_archive_tree_fail_closed`、
`test_trimmed_fastqc_rejects_mismatched_pe_retained_read_counts`。
**升级复核**：比较 ZIP 成员和 BasicStats 字段、版本串、raw `.gz` Filename、PE 两端计数与阶段输出路径；
目标：`test/adapters/test_bulk_rnaseq_qc.py -k fastqc`。
证据：S + T；本轮 ZIP 为合成输入，不是新跑 FastQC。

## UC-05 — Trim Galore → MultiQC Cutadapt 兼容表

**上游**：Trim Galore **2.1.0**，commit `3f6be57a7da52b0b91a2641c6121bff6e34eb6a4`，
[report.rs][TG-report] / [main.rs][TG-main]；NF [模块][NF-trimgalore] / [配置][NF-trim-config]；MultiQC 1.33 [Cutadapt parser][MQ-cutadapt]。
报告的 `cutadapt 4.0` 是 Trim Galore 主动输出的兼容标识，**不是本轮核实的 Cutadapt 安装版本**。
NF config 传入 `--fastqc_args`；Trim Galore 的 SE/PE 分支据此开启 bundled FastQC，不能把 `skip_fastqc` 当作此分支开关。

**本地**：[results_contract.py:86](../../src/encode_pipeline/adapters/bulk_rnaseq/results_contract.py#L86) 的 `trimmed_fastqc_enabled`、[qc.py:1018](../../src/encode_pipeline/adapters/bulk_rnaseq/qc.py#L1018) 的 `_parse_multiqc_cutadapt` →
[status_evidence.py:188](../../src/encode_pipeline/adapters/bulk_rnaseq/status_evidence.py#L188) 的 `parse_cutadapt_rows`。固定九列表头为
`Sample, cutadapt_version, r_processed, r_with_adapters, r_written, bp_processed, quality_trimmed, bp_written, percent_trimmed`
（实际分隔符是 tab）。`cutadapt_version` 必须为 `4.0`，样本键是单端 S 或有效 PE 的 S_1/S_2；
计数为有界整数，禁止重复/外来/缺失行，当前产品禁用 discard-untrimmed，故要求 `r_written == r_processed`。
碱基比例是 adapter/quality 后、长度/配对过滤前的统计；`percent_trimmed` 与精确基数比容差为
`0.000000001` 个百分点。最终 retained reads 另由 trimmed FastQC 提供，不能混淆阶段。

**测试**：[test_bulk_rnaseq_qc.py:1039](../../test/adapters/test_bulk_rnaseq_qc.py#L1039) 的 `test_multiqc_cutadapt_pe_has_explicit_pre_filter_semantics`、
[test_bulk_rnaseq_qc.py:598](../../test/adapters/test_bulk_rnaseq_qc.py#L598) 的 `test_trimgalore_bundled_fastqc_metrics_do_not_require_raw_fastqc`、
`test_multiqc_cutadapt_corrupt_or_unknown_rows_fail_closed`。
**升级复核**：同时核对 Trim Galore 兼容输出和 MultiQC parser，不只比较工具版本字符串；
检查九列、兼容版本、SE/PE 命名、FastQC 内嵌条件和过滤前后计数语义。
目标：`test/adapters/test_bulk_rnaseq_qc.py -k 'cutadapt or trimgalore'`。
证据：S + T；本轮未执行 Trim Galore/MultiQC。

## UC-06 — fastp JSON 的版本字段兼容差异与计数约定

**上游**：fastp **1.0.1** [JsonReporter:76–77][FASTP-json] 将版本写在
**`summary.fastp_version`**。NF **3.26.0**（固定 commit 见上文）的
[FASTP 模块:18、45、61、78][NF-fastp] 直接输出 `--json ${prefix}.fastp.json`，
再由 [fastp 子工作流:112][NF-fastp-workflow] 的 `trim_json = FASTP.out.json` 传递；
已核路径未见将版本键提升至顶层的转换。2026-09-23 在线核对 fastp 固定 tag 源码；
NF 在线获取失败，重新核对已有两份官方源码归档，大小与 SHA 均符合固定 source manifest。

**本地最终字段契约**：[status_evidence.py:400](../../src/encode_pipeline/adapters/bulk_rnaseq/status_evidence.py#L400)
的 `parse_fastp_summary` 先检查 payload 和 summary 为对象，再要求
`summary.fastp_version == "1.0.1"`。没有顶层 fallback：只有顶层版本、
嵌套版本缺失/错误/类型错误均拒绝，顶层值不能掩盖嵌套错误。
不新增未知字段限制；嵌套版本及其他校验通过时，额外顶层字段不参与版本判定。
[qc.py:1130](../../src/encode_pipeline/adapters/bulk_rnaseq/qc.py#L1130) 的 `_parse_fastp` 调用此解析器；
artifact 的 [artifacts.py:1244](../../src/encode_pipeline/adapters/bulk_rnaseq/artifacts.py#L1244)
经 [status_evidence.py:395](../../src/encode_pipeline/adapters/bulk_rnaseq/status_evidence.py#L395)
的 `parse_fastp_retained_reads` 消费相同解析结果。
计数仍消费 `summary.before_filtering/after_filtering.total_reads/total_bases`
与 `filtering_result.passed_filter_reads`：有界非负整数、输入非零、after 不大于 before、
passed 等于 after reads；重复 JSON 键、非有限值、超出原大小/深度限制的输入继续拒绝。
SE/PE 都按 JSON 的 read records 聚合，不额外乘二或当作 pair 数；不从 HTML 提取数字。

**历史缺陷与修复**：PR-5 独立审核只移动版本键、保留其他字段，合法 fastp 配置
在顶层布局下通过；上游嵌套布局使旧解析器抛 `StatusEvidenceError`，公开 QC 返回
`BULK_RNASEQ_QC_INVALID / source_content_invalid`。该探针 **1 passed / 1 failed**；
当时既有 **51 passed** 包含采用顶层版本的 fastp 合成夹具，存在布局盲区，不能否定缺陷。
原证据保留在 `/tmp/helix-pr5-review-6b0n09oq/`；本轮只阅读，未执行会覆盖旧输出的探针。

2026-09-23 按授权修正共享解析器及两份测试文件中的合成 JSON 布局，保留原行为断言。
固定同一测试字节后，修复前 **21 failed / 34 passed**，修复后 **55 passed**；
两份完整消费者测试及执行身份测试合计 **332 passed**。这是当前固定版本、合法 fastp
配置路径的兼容修复；[默认 Trim Galore 配置](../../src/encode_pipeline/adapters/bulk_rnaseq/authoring.py#L212)
不因该字段差异触发。实现与定向验证完成，待独立复验。

**正式测试**：

- [test_bulk_rnaseq_artifacts.py:407](../../test/adapters/test_bulk_rnaseq_artifacts.py#L407)
  起的 `test_fastp_shared_parser_*`：嵌套版本、顶层不能掩盖错误、非对象容器、非法计数、
  重复键、非有限值及原 JSON 限额；接受路径检查四项计数和 retained reads。
- [test_bulk_rnaseq_qc.py:1195](../../test/adapters/test_bulk_rnaseq_qc.py#L1195)
  的 `test_fastp_qc_requires_the_pinned_nested_version`：合法配置经过原校验与公开 QC 提取，
  嵌套布局的 250/150 read records、比例 0.6，以及非法布局的准确错误码。
  [1156 行](../../test/adapters/test_bulk_rnaseq_qc.py#L1156) 的 mixed-layout 测试及
  `test_fastp_and_trimmed_fastqc_retained_reads_reconcile_exactly` 保留 SE/PE、UMI 丢一端的聚合断言。
- [test_bulk_rnaseq_artifacts.py:2830](../../test/adapters/test_bulk_rnaseq_artifacts.py#L2830)
  的 `test_fastp_low_trim_sample_does_not_require_post_filter_fastqc` 实际调用公开 artifact
  发现、读取临时 JSON 并核对状态后的产物集合；
  [2569 行](../../test/adapters/test_bulk_rnaseq_artifacts.py#L2569) 的
  `test_fastp_status_evidence_rejects_impossible_fixed_report` 保留计数矛盾拒绝，
  补顶层独有/顶层掩盖错误嵌套版本的 `BULK_RNASEQ_ARTIFACT_STATUS_INVALID` 断言。
- [test_bulk_rnaseq_qc.py:1784](../../test/adapters/test_bulk_rnaseq_qc.py#L1784)
  的 `test_fastp_status_membership_uses_float_threshold_but_filter_uses_long` 保留阈值边界。

**升级复核**：核对固定版本 JsonReporter 的版本键位置、PE read-record 汇总、passed 与 after
关系，并检查 NF JSON channel 是否仍直接传递。目标：
`test/adapters/test_bulk_rnaseq_artifacts.py test/adapters/test_bulk_rnaseq_qc.py -k fastp`。
不能只改版本字面量或用本地合成夹具代替新版本生产源码。

**证据与边界**：源码核对 + 合成输入执行原解析器、公开 QC 与 artifact 发现。
artifact 路径已实际验证，输入 JSON、状态表及占位产物仍为合成数据；不等于真实科学产物或下载验证。
检查的 PATH、`.local/envs/*/bin` 与 `.snakemake/conda/*/bin` 中未找到 fastp，
没有安装或执行 fastp 二进制、容器、nf-core 完整流程，部署版本仍未核实。
本轮证据 `/tmp/helix-fastp-version-s51wlpg9/`；执行闭包和资格已按正式工具同步，
完整 Gate 留在 roadmap 的集成收尾单元，新身份待验证，不能追认任何历史身份。

## UC-07 — MultiQC failure tables 与 nf-core 阈值接线

**上游**：NF [preprocessing][NF-preprocess]、[fastp 子工作流][NF-fastp-workflow]、
[Trim Galore 子工作流][NF-trimgalore-workflow]、[STAR 子工作流][NF-star-workflow]及
[MultiQC failure tables][NF-multiqc]。这是 NF 的 Groovy Long/Float 比较约定，不是泛化 QC 标准。

**本地**：[artifacts.py:1009](../../src/encode_pipeline/adapters/bulk_rnaseq/artifacts.py#L1009) 的 `_read_sample_status_tables`、[qc.py:522](../../src/encode_pipeline/adapters/bulk_rnaseq/qc.py#L522) 的 `_sample_status` 共用
[status_evidence.py:285](../../src/encode_pipeline/adapters/bulk_rnaseq/status_evidence.py#L285) 的 `parse_status_table` / [status_evidence.py:341](../../src/encode_pipeline/adapters/bulk_rnaseq/status_evidence.py#L341) 的 `reconcile_sample_status`。
固定表头 `Sample\tReads after trimming`、`Sample\tSTAR uniquely mapped reads (%)`；
缺表与空表区分，无 native evidence 时不能凭表自行隐藏样本。
trim 表成员边界为 `<= min_trimmed_reads`，实际下游保留边界为 `>=`，因此恰好阈值的状态行不隐藏分析产物；
mapped 表为 `< min_mapped_reads`。本地还复刻二进制 float32 舍入，Trim Galore 路径为 Float(total)-Float(removed)。
fastp JSON、trimmed FastQC/Cutadapt、完整 STAR log 分别提供可核对的证据，不将缺失或非法值补成零。
fastp native evidence 的版本布局已按 [UC-06](#uc-06--fastp-json-的版本字段兼容差异与计数约定)
修正为 `summary.fastp_version`；原公开 QC 与 artifact 消费路径均通过合成输入执行验证，
既有阈值和 read-record 口径未变，真实工具与部署边界仍保留。

**测试**：[test_bulk_rnaseq_qc.py:1719](../../test/adapters/test_bulk_rnaseq_qc.py#L1719) 的 `test_exact_trim_threshold_status_row_preserves_existing_analysis_metrics`、
[test_bulk_rnaseq_qc.py:1699](../../test/adapters/test_bulk_rnaseq_qc.py#L1699) 的 `test_mapped_status_membership_uses_the_fixed_binary32_threshold`；另有
`test_trimgalore_status_reconciles_the_fixed_binary32_count_route` 和同路径替换 race 检查。
**升级复核**：比较上游 filter/report 两条分支及显式 Long/Float 转换，检查表名、精确标题、等号边界、
native evidence 及禁用路线的旧表处理；目标：
`test/adapters/test_bulk_rnaseq_qc.py -k 'status or threshold'`。
证据：S + T；没有在真实 NF JVM 流程中重测浮点边界。

## UC-08 — STAR Log.final.out 与 template 单位

**上游**：STAR **2.7.11b** [Stats.cpp][STAR-stats] / [ReadAlign_oneRead.cpp][STAR-read]、NF [STAR 模块][NF-star]。
**本地**：[status_evidence.py:562](../../src/encode_pipeline/adapters/bulk_rnaseq/status_evidence.py#L562) 的 `parse_star_log_final` → [qc.py:1219](../../src/encode_pipeline/adapters/bulk_rnaseq/qc.py#L1219) 的 `_parse_star`。
依赖固定完整标题/字段语法、时间/长度/计数/百分比；six-way partition 为 unique、accepted multimapped、
too-many-loci 及三类 unmapped，合计必须为输入 templates。PE 一个 template 是一对，SE 是一个 read。
百分比从整数比率生成，STAR 两位小数百分比仅作 `0.00005` fraction 容差的一致性检查。
`too_many_loci` 不并入 `pure_unmapped_template_fraction`，不能直接沿用其他工具的 mapping fraction 定义。

**测试**：[test_bulk_rnaseq_qc.py:1411](../../test/adapters/test_bulk_rnaseq_qc.py#L1411) 的 `test_star_rounded_percentages_use_the_exact_template_count_partition`，另有
`test_star_fixed_layout_rejects_invalid_non_mapping_fields`、`test_star_qc_accepts_zero_second_infinite_mapping_speed`。
**升级复核**：核对标题/新字段、模板单位、六类分区及零秒速度表示；不得把所有非唯一读段都叫 unmapped。
目标：`test/adapters/test_bulk_rnaseq_qc.py -k star`。
证据：S + T；本轮不验证真实比对或生物学代表性。

## UC-09 — Salmon meta_info.json

**上游**：Salmon **1.10.3** [GZipWriter.cpp][SALMON-writer]、NF [SALMON_QUANT][NF-salmon]。
writer 的 `num_processed/num_mapped` 来自 fragment counters；不等同 BAM record 数。
**本地**：[qc.py:1322](../../src/encode_pipeline/adapters/bulk_rnaseq/qc.py#L1322) 的 `_parse_salmon` 要求 `salmon_version=1.10.3`、`mapping_type=alignment`、
`num_libraries=1`；核对 `num_processed`、`num_mapped`、`percent_mapped`，输出 fragments/count 和 fraction。
处理数为正、mapped 不超 processed，报告百分比与整数比在 `1e-12` fraction 包络内一致。
路径由 [artifacts.py:1341](../../src/encode_pipeline/adapters/bulk_rnaseq/artifacts.py#L1341) 的 `_expected_specs` 绑定 `star_salmon/<sample>/aux_info/meta_info.json`；不消费 cmd_info/log 中的私有内容。

**测试**：[test_bulk_rnaseq_qc.py:1470](../../test/adapters/test_bulk_rnaseq_qc.py#L1470) 的 `test_salmon_meta_info_counts_fragments_not_records_reads_or_alignments`，另有
`test_salmon_rejects_a_percent_that_disagrees_with_fragment_counters`、duplicate/non-finite/depth 检查。
**升级复核**：检查 alignment 而非 pseudoalignment 路径、字段版本、library 数与计数单位；目标：
`test/adapters/test_bulk_rnaseq_qc.py -k salmon`。
证据：S + T；没有实跑 Salmon 或复核其他 mode。

## UC-10 — featureCounts summary 的状态目录与 BAM 名

**上游**：Subread/featureCounts **2.0.6**；已核对 NF [模块][NF-featurecounts]与 [config][NF-featurecounts-config]：
PE 加 `-p`，未加 `--countReadPairs`。本轮尝试固定 2.0.6 `readSummary.c` 的官方 GitHub URL 返回 404，
**未独立复核该版本底层计数实现**；不能把本地契约测试当作上游 PE 计数实测。

**本地**：[qc.py:1371](../../src/encode_pipeline/adapters/bulk_rnaseq/qc.py#L1371) 的 `_parse_featurecounts`、[qc.py:1423](../../src/encode_pipeline/adapters/bulk_rnaseq/qc.py#L1423) 的 `_final_bam_name` 与封闭 `_FEATURECOUNTS_STATUSES`。
TSV 必须两列，首列标题 `Status`，第二列为当前路径的 `<sample>.sorted.bam`、
`.markdup.sorted.bam` 或 `.umi_dedup.sorted.bam`；要求完整、唯一、非负状态目录。
本地公开契约将 SE/PE 数值都标记为 reads/count，assigned fraction 由 Assigned/全部状态计数生成，不能自行改成 fragments。

**测试**：[test_bulk_rnaseq_qc.py:1990](../../test/adapters/test_bulk_rnaseq_qc.py#L1990) 的 `test_featurecounts_fixed_route_counts_individual_reads_for_both_layouts`、
`test_featurecounts_requires_the_fixed_status_catalog`。
**升级复核**：先取得可核对的目标版本实现/手册，验证 `-p` 与 `--countReadPairs` 的语义；再核对状态名、
BAM 列名及 UMI/markdup 分支。目标：`test/adapters/test_bulk_rnaseq_qc.py -k featurecounts`。
证据：本地 S + T、NF wrapper S；底层 Subread 固定版本源码和真实输出未核实。

## UC-11 — Picard 经 MultiQC 转出的 duplication 表

**上游**：Picard **3.4.0** [DuplicationMetrics][PICARD-metrics]，MultiQC **1.33** [MarkDuplicates parser][MQ-picard]，
NF [模块][NF-picard]。MultiQC 有 multiple-library 合并/重算分支，机器表不是任意原始 Picard 文本的透传。
**本地**：[qc.py:1431](../../src/encode_pipeline/adapters/bulk_rnaseq/qc.py#L1431) 的 `_parse_picard_multiqc` 要求精确 11 列（`Sample`、`LIBRARY`、七项计数、
`PERCENT_DUPLICATION`、`ESTIMATED_LIBRARY_SIZE`），每个期望样本恰一行；fraction 原值在 [0,1]，不再次除以 100。
以 `(unpaired duplicates + 2*pair duplicates)/(unpaired examined + 2*pairs examined)` 核对，
兼容六位小数误差上限 `0.0000005`；library size 为 `''` 或 `-` 时不生成数值，不补零。

**测试**：[test_bulk_rnaseq_qc.py:2022](../../test/adapters/test_bulk_rnaseq_qc.py#L2022) 的 `test_picard_duplication_reconciles_counts_with_six_decimal_serialization`、
[test_bulk_rnaseq_qc.py:1882](../../test/adapters/test_bulk_rnaseq_qc.py#L1882) 的 `test_featurecounts_picard_and_rseqc_use_closed_machine_sources`。
**升级复核**：复核列顺序、library 合并是否改变行键/重算精度、零 reads 的上游省略及缺失值表示；目标：
`test/adapters/test_bulk_rnaseq_qc.py -k picard`。
证据：S + T；当前列举测试未实际执行多 library 的 MultiQC 合并，不宣称这条假设已端到端覆盖。

## UC-12 — RSeQC native 文本、TIN 分数和文件名

**上游**：RSeQC **5.0.4** 官方 [源码发行包][RSEQC-source]，SHA-256
`b7f3996f3de0b0b0a09eec949281a8f3e665a20827fcb3cbbd7546b94574a088`（与 results contract 相同）。
本轮只读取 `scripts/bam_stat.py`、`infer_experiment.py`、`read_distribution.py`、`tin.py` 及 `src/qcmodule/SAM.py`，没有安装。
对应 NF [bamstat][NF-bamstat] / [inferexperiment][NF-infer] / [readdistribution][NF-distribution] / [tin][NF-tin] 接线已核对。

**本地**：[qc.py:1619](../../src/encode_pipeline/adapters/bulk_rnaseq/qc.py#L1619) 的 `_rseqc_bam_stat_counts`、[qc.py:1702](../../src/encode_pipeline/adapters/bulk_rnaseq/qc.py#L1702) 的 `_parse_rseqc_infer`、
[qc.py:1770](../../src/encode_pipeline/adapters/bulk_rnaseq/qc.py#L1770) 的 `_parse_rseqc_distribution`、[qc.py:1906](../../src/encode_pipeline/adapters/bulk_rnaseq/qc.py#L1906) 的 `_parse_rseqc_tin`：

| 消费物 | 依赖约定与升级点 |
| --- | --- |
| bam_stat `.bam_stat.txt` | 按固定标签和 read-count 小计解析。上游 `%-40s` 是最小宽度，长标签 `Proper-paired reads map to different chrom:` 后可能直接跟数字，本地该项用 `\s*`，其余多数用 `\s+`；不得恢复成要求填充空格。核对总记录与 unique、mate/strand/splice 小计。 |
| infer_experiment 文本 | 区分 SingleEnd/PairEnd 与固定 orientation 字符串；原值为 fraction，未知 orientation 不接受。升级核对标签和布局，不按列号猜样本。 |
| read_distribution 文本 | 按固定 group/Total Tags 等语法及区间小计关系核对；tag 与 read 指标不得混称。升级复核计数重叠/互斥关系，而非只保留总和。 |
| TIN `.summary.txt` | 精确四列 `Bam_file, TIN(mean), TIN(median), TIN(stdev)`（tab）、一数据行，BAM 名精确匹配。mean/median 0–100，population stdev ≤50，公共单位 `score`；只容纳既定数值量化包络，不改为百分比。 |

`tin.py` 用全局、区分大小写的 basename `.replace('bam','')`，NF wrapper 期望 bam.baseName；
[validation.py:604](../../src/encode_pipeline/adapters/bulk_rnaseq/validation.py#L604) 的 `_validate_normalized_result_identities` 已在有效 TIN 路径拒绝含小写 `bam` 的样本 ID，保留原输入不改名。
[results_contract.py:97](../../src/encode_pipeline/adapters/bulk_rnaseq/results_contract.py#L97) 的 `effective_rseqc_modules` 的 CSI 路径会移除 inner_distance/read_distribution/tin；这是本地已声明的固定路线策略，
本轮未重新下载 NF utility 实现核对 CSI 接线。inner-distance 等仅作产物、不提取上述 QC 数字的模块，不假称已有数字 parser。

**测试**：[test_bulk_rnaseq_qc.py:1882](../../test/adapters/test_bulk_rnaseq_qc.py#L1882) 的 `test_featurecounts_picard_and_rseqc_use_closed_machine_sources`、
[test_bulk_rnaseq_qc.py:2161](../../test/adapters/test_bulk_rnaseq_qc.py#L2161) 的 `test_rseqc_se_bam_stat_accepts_fixed_zero_paired_fields`（含长标签无空格）、
[test_bulk_rnaseq_qc.py:2117](../../test/adapters/test_bulk_rnaseq_qc.py#L2117) 的 `test_rseqc_tin_population_standard_deviation_allows_exact_global_bound`；
另有 `test_rseqc_infer_rejects_unknown_orientation_labels`、`test_rseqc_read_distribution_reconciles_outer_windows_and_total` 和 validator 的 TIN 命名测试。
**升级复核**：用新版本原始四类输出核对标签、空白、PE/SE、单位和缺值，再复核 TIN basename 与 CSI 门控；目标：
`test/adapters/test_bulk_rnaseq_qc.py -k rseqc` 及
`test/adapters/test_bulk_rnaseq_adapter.py -k rseqc_tin`。
证据：S + T；未运行 BAM 分析，也未完整覆盖上游所有 RSeQC 模块。

## UC-13 — ENCODE 的 MultiQC 1.35 阶段与文件名配置

**上游与锁**：此项属于本仓库 ENCODE Snakemake workflow，不是 nf-core/rnaseq。
[MultiQC 1.35 cleaner][MQ135-base] / [defaults][MQ135-defaults] / [FastQC 模块][MQ135-fastqc]。
[workflow/envs/multiqc.yml:8](../../workflow/envs/multiqc.yml#L8) 只声明 `multiqc >=1.20,<2`，
而 [multiqc.lock:122](../../workflow/envs/multiqc.lock#L122) 锁定 `multiqc-1.35-pyhdfd78af_1.conda`，
URL 片段 `cdb20309681ba3ce8f52c110e214d4f3` 是包锁摘要，不是容器 digest；本机安装/实际 worker 使用版本本轮未核实。

**本地消费者**：[report.smk:144](../../workflow/rules/report.smk#L144) 的 `rule multiqc` 使用
[workflow/multiqc_config.yaml:6](../../workflow/multiqc_config.yaml#L6)：额外清理 `.final/.sorted/.blacklist_filtered`；
复制 1.35 清洗列表但排除会截断 `.trimmed.R1` 的 `.trim`；仅 `fastqc/zip` 用文件名作样本名。
[config:174](../../workflow/multiqc_config.yaml#L174) 以 path_filters/anchor 分开 raw/trimmed FastQC 和 pre-filter/final samtools；
General Stats 列名也分阶段。这里的消费者是报告生成，不是 bulk 的 S_1→S Read 1 图或其公开下载政策。
样本后缀清洗理论上可改变名称；本轮未证明一个当前可接受的 ENCODE 输入必然碰撞，不据此登记新缺陷。

**已有测试**：[test_staged_fastqc_and_multiqc_135_report_contract:29](../../test/real_execution/test_staged_fastqc_multiqc.py#L29)
是实际工具测试，检查 HTML anchors、`sample.one/two` 的 raw/trimmed R1/R2 机器表及两套 samtools 输出。
其中 FastQC/Trim Galore/MultiQC 为真实工具调用，samtools flagstat/idxstats 是测试写入的合成文本，
只验证报告分组，不证明 samtools 实际执行。它不穷举所有样本名，也不证明任意配置都无冲突。
**本轮未运行**，不把旧通过结果作为新实测。

**升级复核**：比较新 defaults 与本地完整复制列表，确认只保留有意差异；核对
`use_filename_as_sample_name`、模块 anchor/path_filters、机器表文件名及 General Stats 同名列行为。
具备锁定工具与隔离 fixture 条件后，使用公共 bootstrap 前缀的目标
`test/real_execution/test_staged_fastqc_multiqc.py::test_staged_fastqc_and_multiqc_135_report_contract -m real_execution`，
显式覆盖 pytest 默认排除真实执行测试的 marker 选择。
本轮不执行该命令，以遵守“不启动真实工作流”。证据：S，真实工具验证待未来升级时执行。

## 可得的上游容器坐标

以下来自 NF 固定 commit 的模块，并与本地 inventory 核对；其中 `biocontainers/...` 简写
按 [nextflow.config][NF-config] 的 `docker.registry = 'quay.io'` 补齐。它们是上游标签坐标，
不是本次已检查的运行镜像 digest。实际 digest 应在另一次获授权的运行时/升级验证中从受控部署锁获取。

| process | 固定的 image_coordinate |
| --- | --- |
| CUSTOM_MULTIQCCUSTOMBIOTYPE | `quay.io/biocontainers/python:3.12.12` |
| FASTP | `community.wave.seqera.io/library/fastp:1.0.1--c8b87fe62dcc103c` |
| FASTQC | `quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0` |
| MULTIQC | `community.wave.seqera.io/library/multiqc:1.33--ee7739d47738383b` |
| PICARD_MARKDUPLICATES | `community.wave.seqera.io/library/picard:3.4.0--e9963040df0a9bf6` |
| RSEQC_BAMSTAT | `community.wave.seqera.io/library/rseqc_r-base:2e29d2dfda9cef15` |
| RSEQC_INFEREXPERIMENT | `community.wave.seqera.io/library/rseqc_r-base:2e29d2dfda9cef15` |
| RSEQC_READDISTRIBUTION | `community.wave.seqera.io/library/rseqc_r-base:2e29d2dfda9cef15` |
| RSEQC_TIN | `community.wave.seqera.io/library/rseqc_r-base:2e29d2dfda9cef15` |
| SALMON_QUANT | `quay.io/biocontainers/salmon:1.10.3--h6dccd9a_2` |
| STAR_ALIGN | `community.wave.seqera.io/library/htslib_samtools_star_gawk:ae438e9a604351a4` |
| SUBREAD_FEATURECOUNTS | `quay.io/biocontainers/subread:2.0.6--he4a0461_2` |
| TRIMGALORE | `community.wave.seqera.io/library/trim-galore:2.1.0--27e6376b8f6c1872` |

## 交付与返修证据、剩余边界

初次交付证据根目录：`/tmp/helix-pr5-ledger-si4g0_12/`，保留原样。
`logs/upstream-fetch*.json` 保存 URL、HTTP 结果、大小和 SHA-256；上游文件只保存于该目录。
官方链接核对不等于工具执行；初始受限网络失败和 Subread 404 也留有记录。
`logs/source-identity-check.json` 记录 NF 下载件与 source manifest、vendored schema、
MultiQC 1.33 三份清洗/分组源码、RSeQC 源码发行包及 13 个容器坐标的核对结果。
`selected-tests.json` / `logs/targeted-contracts-command.json` 精确记录初次交付实际选择的 23 个测试函数，
展开 **51 passed / 0 failed / 0 skipped**。根 fixture 与所选 helper 已检查；访问审计没有保护文件、数据库或网络访问。
所选测试只消费合成数据或固定契约资源，没有调用 runtime fixture、真实工作流或邮件。

独立审核证据：`/tmp/helix-pr5-review-6b0n09oq/`；文档返修证据：
`/tmp/helix-pr5-doc-revision-3e8r8g68/`。返修纠正结果子类归属，并记录 UC-06 的当前 fastp 兼容差异；
撤回原交付报告中“没有确认需要新增修复条目的生产问题”的笼统结论。
是否及如何纳入既定修复单元待用户裁决，不据此新增 PR 或插队。本次返修仅文档检查，未运行行为测试。

尚未核实：部署 OCI digest、实际 worker/本机工具版本；Subread 2.0.6 底层实现；NF CSI utility 本轮未重核；
真实 MultiQC 多 library/完整 HTML 和新版本兼容性。没有因这些边界擅自新增修复项；科学政策与 PR 队列不变。
本轮只核对上述实际消费者，不承诺覆盖所有上游工具，也不把测试存在写成对全部假设的证明。

[NF-schema]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/nextflow_schema.json
[NF-samples]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/assets/schema_input.json
[NF-config]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/nextflow.config
[NF-main]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/workflows/rnaseq/main.nf
[NF-multiqc]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/subworkflows/local/multiqc_rnaseq/main.nf
[NF-multiqc-config]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/workflows/rnaseq/assets/multiqc/multiqc_config.yml
[NF-fastqc]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/modules/nf-core/fastqc/main.nf
[NF-fastp]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/modules/nf-core/fastp/main.nf
[NF-trimgalore]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/modules/nf-core/trimgalore/main.nf
[NF-trim-config]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/conf/modules/trimgalore.config
[NF-preprocess]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/subworkflows/nf-core/fastq_qc_trim_filter_setstrandedness/main.nf
[NF-fastp-workflow]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/subworkflows/nf-core/fastq_fastqc_umitools_fastp/main.nf
[NF-trimgalore-workflow]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/subworkflows/nf-core/fastq_fastqc_umitools_trimgalore/main.nf
[NF-star-workflow]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/subworkflows/local/align_star/main.nf
[NF-star]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/modules/nf-core/star/align/main.nf
[NF-salmon]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/modules/nf-core/salmon/quant/main.nf
[NF-featurecounts]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/modules/nf-core/subread/featurecounts/main.nf
[NF-featurecounts-config]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/conf/modules/featurecounts.config
[NF-picard]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/modules/nf-core/picard/markduplicates/main.nf
[NF-bamstat]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/modules/nf-core/rseqc/bamstat/main.nf
[NF-infer]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/modules/nf-core/rseqc/inferexperiment/main.nf
[NF-distribution]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/modules/nf-core/rseqc/readdistribution/main.nf
[NF-tin]: https://github.com/nf-core/rnaseq/blob/e7ca46272c8f9d5ceee3f71759f4ba551d3217a4/modules/nf-core/rseqc/tin/main.nf
[MQ-base]: https://github.com/MultiQC/MultiQC/blob/5953b5417ccb70bf4a2309562d43015fced8b585/multiqc/base_module.py
[MQ-defaults]: https://github.com/MultiQC/MultiQC/blob/5953b5417ccb70bf4a2309562d43015fced8b585/multiqc/config_defaults.yaml
[MQ-table]: https://github.com/MultiQC/MultiQC/blob/5953b5417ccb70bf4a2309562d43015fced8b585/multiqc/plots/table_object.py
[MQ-cutadapt]: https://github.com/MultiQC/MultiQC/blob/5953b5417ccb70bf4a2309562d43015fced8b585/multiqc/modules/cutadapt/cutadapt.py
[MQ-picard]: https://github.com/MultiQC/MultiQC/blob/5953b5417ccb70bf4a2309562d43015fced8b585/multiqc/modules/picard/MarkDuplicates.py
[MQ135-base]: https://github.com/MultiQC/MultiQC/blob/87e504bb77e3d94687dd980638724b7ef8e1e92f/multiqc/base_module.py
[MQ135-defaults]: https://github.com/MultiQC/MultiQC/blob/87e504bb77e3d94687dd980638724b7ef8e1e92f/multiqc/config_defaults.yaml
[MQ135-fastqc]: https://github.com/MultiQC/MultiQC/blob/87e504bb77e3d94687dd980638724b7ef8e1e92f/multiqc/modules/fastqc/fastqc.py
[FQC-basic]: https://github.com/s-andrews/FastQC/blob/v0.12.1/uk/ac/babraham/FastQC/Modules/BasicStats.java
[TG-report]: https://github.com/FelixKrueger/TrimGalore/blob/3f6be57a7da52b0b91a2641c6121bff6e34eb6a4/src/report.rs
[TG-main]: https://github.com/FelixKrueger/TrimGalore/blob/3f6be57a7da52b0b91a2641c6121bff6e34eb6a4/src/main.rs
[FASTP-json]: https://github.com/OpenGene/fastp/blob/v1.0.1/src/jsonreporter.cpp
[STAR-stats]: https://github.com/alexdobin/STAR/blob/2.7.11b/source/Stats.cpp
[STAR-read]: https://github.com/alexdobin/STAR/blob/2.7.11b/source/ReadAlign_oneRead.cpp
[SALMON-writer]: https://github.com/COMBINE-lab/salmon/blob/v1.10.3/src/GZipWriter.cpp
[PICARD-metrics]: https://github.com/broadinstitute/picard/blob/3.4.0/src/main/java/picard/sam/DuplicationMetrics.java
[RSEQC-source]: https://files.pythonhosted.org/packages/f8/04/f32be86c3e14fa77b0dae3a9cdea9e5795642960ade79e3d463664101a4d/rseqc-5.0.4.tar.gz

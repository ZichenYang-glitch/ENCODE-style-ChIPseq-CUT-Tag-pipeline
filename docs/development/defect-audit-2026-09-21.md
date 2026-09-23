# HelixWeave 缺陷审计：2026-09-21 原记录与后续处置

本文件保留原审计的 S1–S20、P1–P6、D1–D2 编号和证据归属，
并于 2026-09-23 按已接受的复核及 PR-4、PR-8～PR-16 修复更新。
它不是“当前仍有 28 个缺陷”的清单，也不是所有运行环境均已验收的声明。
PR-17 只同步文档，不修改科学行为。当前顺序与验收边界以
[roadmap](workflow-platform-agent-roadmap.md#已确认的维护顺序与审核边界2026-09-22) 为准。

## 证据、统计与默认值

原正文共 **28 条**：V 13、R 13，D1/D2 正文没有 V/R 标签。
原汇总表只有 27 条，漏了 D2，且把正文无标签的 D1 写成 V。
下表已补齐，标签仍指**原作者当时的证据归属**：

- V：原作者实际复现或读过具体源码；包含源码阅读，不等于全部实测。
- R：原作者引用并行审查者的报告，未自行复现；不代表后续一直未验证。
- 后续证据区分 A（真实工具/原脚本或真实临时 HTTP）、B（dry-run、隔离调度；
  如有替身另外注明）、C（仅源码核对）。不能把 dry-run 当作 shell 成功。
- N1/N2/N3 是后续补充，不计入原 28 条。部分确证可能同时含默认值、机制、
  影响范围或政策判断错误，不能全部计作正确/误报，也不据此生成单一准确率。

原审计记载 Snakemake 8.30.0、MACS3 3.0.4 和本地 conda 中的 `run_spp.R`。
这些是 **2026-09-21 实验环境记录**，不是 PR-17 实测版本或部署证明。
原命令与输出补充见 `/tmp/helixweave-defect-review-Nobias/{review,commands}.md`；
补强复核见 `/tmp/helix-independent-corrections.IEHFno/`。修复前的失败输入、
输出与工具身份由各实施/独立报告保存，不能继续称为当前版本可复现配方。
PR-17 主要做源码与已有记录核对；本轮小型 consensus 原脚本验证另见 roadmap。
这些临时材料是补充，以下条目同时提供仓库源码、契约和正式测试入口。
原审计的“只读”描述属于当时过程，不适用于本轮文档编辑；本文件也不背书其他
会话的全过程保护文件访问声明。本轮未读取用户指定的五项保护内容。

默认真源与调度不能混为一谈：

| 配置 | 当前默认与资格 | 源码依据 |
| :-- | :-- | :-- |
| QC 11 项 | 7 true：blacklist_filter、frip、library_complexity、nrf_pbc、signal_tracks、summary、cuttag_fragment_size；4 false：cross_correlation、preseq_complexity、picard_metrics、tss_enrichment | [`config/qc.py`](../../src/encode_pipeline/config/qc.py):11 |
| reproducibility | 父级 enabled=false；consensus 子级 enabled=true 不能越过父级 | [`config/reproducibility.py`](../../src/encode_pipeline/config/reproducibility.py):57 |
| consensus 默认目标 | 还需 replicate_analysis、适用 assay/mode 和至少两个 treatment bioreps；min_replicates 控制保留支持度 | [`metadata.smk`](../../workflow/rules/metadata.smk):380、[`targets.smk`](../../workflow/rules/targets.smk) 的 `_consensus_targets` |
| ChIP narrow IDR | chipseq_idr=false；开启需 replicate_analysis 且恰好两个适用 bioreps；不由 reproducibility 父级替代 | [`config/validator.py`](../../src/encode_pipeline/config/validator.py):203、[`metadata.smk`](../../workflow/rules/metadata.smk) |
| 扩展 IDR | reproducibility 父级与各 IDR 子开关、replicate_analysis 和模式资格共同控制 | [`workflow/Snakefile`](../../workflow/Snakefile):75、[`config/reproducibility.py`](../../src/encode_pipeline/config/reproducibility.py) |
| MultiQC | multiqc=true，独立于上述 QC 表；按启用规则与样本生成报告 | [`config/validator.py`](../../src/encode_pipeline/config/validator.py):170、[`report.smk`](../../workflow/rules/report.smk):143 |

表中描述默认目标选择。显式请求某规则产物是另一条可达路径，不能用默认未调度
推断规则永远不可达。通用 summary 还会请求自己的组件指标。

## 28 条处置总表

“已接受”限于 roadmap 中实现和定向验证范围，不包含尚欠的 Gate 或完整科学流程。

| ID | 原正文标签 | 原主张简记 | 复核后当前处置 |
| :-- | :-- | :-- | :-- |
| S1 | V | mixed MNase 项目默认 DAG 失败 | 合法 mixed treatment 条件下成立；PR-8 已修并接受 |
| S2 | V | consensus 默认开启且失败 | 默认值错误；opt-in 空参数缺陷由 PR-10 修复并接受 |
| S3 | V | 传递性合并是算法缺陷 | 连通分量是显式设计；政策待决，文档澄清 |
| S4 | V | speak=0 使整个互相关计算失效 | 覆盖选峰，不停止曲线；PR-11 已修并接受 |
| S5 | V | blacklist 未用于交付物 | raw/filtered 并存、消费者不同；PR-17 澄清，不据此改变科学输入 |
| S6 | V | 去重后 NRF/PBC 恒退化 | 非恒等；输入阶段缺陷由 PR-9 修复并接受 |
| S7 | R | preseq 去重输入、PE 缺 -P | 两个问题分开确证；PR-9 已修并接受 |
| S8 | V | FRiP 未过滤、未平移必偏低 | 条件满足时双过滤，偏差可双向；PR-17 澄清 |
| S9 | R | AND/OR 差异意味着结果必更少 | 算法/候选集也不同，不能推出包含关系；政策待决 |
| S10 | R | 缺 borderline 是实现缺陷 | 本地二态严格 <2 契约；政策待决，与 N1 分开 |
| S11 | R | 缺 oracle 必错且添加可保证嵌套 | 兼容差异；无真实 IDR 拟合证据支持保证，政策待决 |
| S12 | V | p/q=-1 的下游用途及全支持例中的恒定 score | 撤回实现缺陷；score 与支持比例相关，-1 合法 |
| S13 | R | 多候选 estFragLen 解析为 NA | 非必然但已由真实 R 产物确证；PR-11 已修并接受 |
| S14 | V | 本地窗口必须随 S4 修复 | ENCODE wrapper/WDL 也设窗口；政策待决，未随 PR-11 修改 |
| S15 | R | 担忧线粒体 reads 影响指标，并建议先定政策 | 损害证据不足；不同 assay/版本不能混称统一标准 |
| S16 | R | 未声明 log 导致永不重建 | 撤回永不重建；首轮依赖缺失由 PR-12 修复并接受 |
| S17 | R | MNase insert-size 数值常为 NA | 实为路径字段、Picard opt-in；PR-14 已修依赖并接受 |
| S18 | R | ATAC 完全无分布生产路径 | 撤回中心指控；Picard histogram 路径存在，默认关闭 |
| S19 | V | 项目表未接 MultiQC，相关指标全不可见 | 特定表未接入成立；完整 HTML 可见性证据不足 |
| S20 | R | fraction_lt_120 是死指标 | 工作区 TSV 有该字段；部分消费者未接入，平台可下载未证实 |
| P1 | V | trusted 分支差异即越界漏洞 | 未证明不可信输入可达攻击路径；信任边界待决 |
| P2 | V | events/logs limit 无上限 | 真实 SQLite/HTTP 溢出确证；PR-15 已修并接受，Gate 待补 |
| P3 | R | 500 envelope 未声明 | PR-16 补实际 operation 声明并接受，Gate/流式成功边界保留 |
| P4 | R | 响应未显式区分本次取消与既有终态 | 撤回；幂等返回 canonical status |
| P5 | R | 所有通知失败无痕 | 普通失败已有事件；PR-4 补内外静默边界及 doctor 日志并接受 |
| P6 | V | 四类读取无认证 | 撤回；父 router 鉴权，匿名请求 401 |
| D1 | — | BigWig planned 文字过时 | PR-17 按已有条件转换功能修正文案 |
| D2 | — | blacklist 文档证明全局过滤违约 | 不足以证明该承诺；PR-17 明确副本及消费者 |

## 科学条目：原主张、纠正与当前依据

### S1 — 混合项目的通用 summary 依赖（PR-8 已接受）

原 V 记录称无需 FASTQ 就能触发 MNase MACS3 guard。复核纠正：缺 FASTQ 可先
触发 `MissingInputException`；合法输入、mixed treatment（存在 peak assay）且
summary 开启时，旧项目汇总遍历所有 treatment 才会将 MNase 拉入
`qc_summary → peak_counts → macs3_callpeak`。不能泛称任何混合项目同样失败。

当前 [`qc.smk`](../../workflow/rules/qc.smk):924 的 `project_qc_summary` 遍历
`PEAK_SAMPLE_IDS`，与 [`targets.smk`](../../workflow/rules/targets.smk):118 的
通用 summary 目标一致；MNase guard 和专用 summary 保留。
正式证据：[`test_project_qc_summary_dag.py`](../../test/workflow/test_project_qc_summary_dag.py)
的 `test_project_summary_assay_dependencies`，覆盖合法 mixed 默认/true/false、
MNase-only、peak-only 及缺 FASTQ。PR-8 独立 123 项通过含六组 DAG；这是 B，
不代表测序或 shell 端到端执行。

### S2 — narrow consensus 空参数（PR-10 已接受）

原 V 把 `reproducibility` 写成默认开启，且混淆 DAG 与执行失败。父级实际默认
false；启用父级、consensus、replicate_analysis 且样本合格后，ChIP/ATAC narrow
可返回空 final_output。旧 shell `--final-output {params.final_output:q}` 少一个参数
token，使 argparse 报 `argument --final-output: expected one argument`。

当前 [`consensus.smk`](../../workflow/rules/consensus.smk):300 使用
`--final-output={params.final_output:q}`；非空路径仍有引号保护。
`final_output` 是 summary 字段，不能补造路径或声称该脚本写了第二个产物。
[`test_consensus_execution.py`](../../test/workflow/test_consensus_execution.py)
实际执行原 narrow 规则和原脚本，含 ChIP/ATAC 空值、CUT&Tag 非空及空格路径、
非法峰反例和默认门控。PR-10 独立 174 项通过；A 使用预置合法峰，不证明峰生产者。

### S3 — 连通分量与支持语义（政策待决）

原 V 的传递性合并示例可得到 400 bp 区间，但不是违反现有契约的证据。
[`compute_consensus.py`](../../scripts/compute_consensus.py):354、401、434、499
按双向重叠阈值建图，取连通分量，按不同 biorep 数筛选，输出最小 start 到最大 end。
两重复时也不等同逐碱基 intersection；合并区间不保证每个碱基获 N 个重复支持。
这是 C 及已有原脚本测试支持的显式算法，是否更改由科学政策决定。
见 `test_overlap_components_respect_interval_topology` 和
[修正后的政策说明](../reproducibility-policy.md#22-consensus)。

### S4 — 自动选峰被 speak 覆盖（PR-11 已接受）

原 V 将 `-speak=0` 描述成停止整条曲线计算，机制错误。PR-11 核对的已安装 `run_spp.R`
先计算曲线并按相关系数排列候选，再由 speak 覆盖；主候选用于 NSC/RSC。
PR-11 只从 [`qc.smk`](../../workflow/rules/qc.smk):762 移除固定 speak，保留
BAM、`-x=-500:15`、`-rf`、线程及失败语义。移除后不保证 NSC/RSC 必改善或多峰。

[`test_cross_correlation.py`](../../test/real_execution/test_cross_correlation.py)
保存同 BAM、其他参数一致的真实 R 对照和原规则/summary 执行。原审计同时改
窗口的数字不是该单变量实验，不能混作同一次输出。PR-11 独立 199 项回归及
4 项真实工具测试通过；完整锁定 R 环境、远端 CI 等边界仍见 roadmap。S14 未改。

### S5 — blacklist 副本与消费者（PR-17 文档澄清）

原 V 将“未全局替换交付物”与“没有过滤”混同，并误称 peak_counts 消费过滤 BAM。
当前 [`qc.smk`](../../workflow/rules/qc.smk):33、72、97、160、210 源码核对（C）：
`_has_blacklist_qc` 同时要求开关及 `BLACKLIST_SAMPLES` 资格；后者是配置了
非空 blacklist 资源的 treatment 样本（[`metadata.smk`](../../workflow/rules/metadata.smk):470）。
过滤规则生成额外 BAM/峰副本，原始产物保留。

`_frip_inputs` 满足条件选两个过滤副本，否则选 final BAM 与原峰目录；选择建立
DAG 生产者依赖，不是文件未生成时静默回退。`_peak_counts_inputs` 只选原峰与
可选过滤峰。默认 CPM coverage（[`common.smk`](../../workflow/rules/common.smk):396）
用 final BAM；MACS3 峰调用及 FE/ppois 分别消费其 BAM、pileup/background，
不使用这些 QC 过滤副本。final BAM 是否去重取决于策略。
见 [配置说明](../configuration.md#qc-block)；不据文案改变科学输入。

### S6 — NRF/PBC 的输入阶段（PR-9 已接受）

原 V 用已无重复的输入演示 1/1/NA，却推成去重后恒等。脚本片段键与
samtools/Picard 重复键不同（包括 clipping），真实去重后存在非恒等反例；
无重复文库的 `distinct == total` 本身合法，PBC2 零分母保持 NA。

当前 [`qc.smk`](../../workflow/rules/qc.smk):326 消费配置对应的 `{MAPQ_TAG}.bam`，
位于过滤后、duplicate_handling 前；不改变指标键、过滤、FRiP 或峰输入。
[`test_library_complexity.py`](../../test/real_execution/test_library_complexity.py)
的真实 samtools 去重链、clipping/空输入/合法 1/1/NA 与原规则执行（A），加上
[`test_complexity_inputs.py`](../../test/workflow/test_complexity_inputs.py) 的
SE/PE、remove_dup 与配置 MAPQ 对照（B），支持输入修复。不能由此规定所有库应有重复。

### S7 — preseq 输入与 PE 模式（PR-9 已接受）

原 R（输入前提 V）把去重影响泛化成必输出低曲线。频数信息不足可直接导致
非零退出，不能保证外推成功；旧 SE BAM 加 `-P` 探针不是合法 PE 证据。
后续合法双端 BAM 的同输入单变量对照补强了 `-P` 统计差异。

当前 [`qc.smk`](../../workflow/rules/qc.smk):791 使用过滤后保留重复 BAM，按已验证
layout 给 PE `-P`，SE 不加。preseq 默认关闭。生产参数没有为夹具加入 `-Q`，
失败仍失败，不伪造曲线。上述真实工具测试的
`test_preseq_pair_mode_is_the_only_difference`、
`test_uninformative_preseq_input_fails_without_fabricating_a_curve` 分别覆盖模式和失败。
历史 `-Q/-e/-s` 诊断实验不等于原生产规则执行；PR-9 报告分别列出两者。

### S8 — FRiP 计数与过滤（文档澄清，兼容政策未改）

原 V 源码主张“blacklist-inclusive 且未平移必低估”均过强。条件成立时输入双过滤，
否则 raw/final；[`calc_frip.py`](../../scripts/calc_frip.py):24、37、93 用
`samtools view -c` 记录数和 `bedtools intersect -u`，不平移，不按 PE 片段计数。
复核 A 的相同 read/两种峰边界对照可使平移后 FRiP 分别升高或降低，不代表真实误差分布。
特定 ENCODE ChIP 版本的 tagAlign/shift 接线与本地不同，不足以判定必然科学错误。
见 [FRiP 说明](../qc-interpretation.md#frip-fraction-of-reads-in-peaks)。

### S9 — overlap AND/OR（政策待决）

原 R 从本地 AND、上游 OR 推出最终峰必更少，推论不成立。
本地用双向条件建图；ENCODE-DCC **ChIP v2.2.2** 的
[overlap 实现](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/v2.2.2/src/encode_task_overlap.py)
是 pooled 候选依次与两重复相交，单次判定使用任一方向比例。候选空间及算法层次均不同，
不能推出最终集合包含关系。本轮核对官方版本源码（C），未做上游流程执行。
原复核归档 `upstream/sources.json` 的六项**不含 overlap**；本轮在线核对不补成历史归档已完整。

### S10 — 二态判级（政策待决；N1 已修）

原 R 将无 borderline 视为实现错误。本地两个 summary 明确两比率有效且严格 <2
才 pass；否则 fail。ENCODE-DCC ChIP v2.2.2
[reproducibility 源码](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/v2.2.2/src/encode_task_reproducibility.py)
使用 >2 边界及 borderline/fail 分层，属于特定版本兼容政策，不是统一标准。
本轮核对已有版本源码归档（C），未执行上游。本地舍入错误 N1 已由 PR-13 修，
不引入 borderline、不改等号；原八列/十五列与复制行为保持。

### S11 — IDR oracle（兼容政策待决）

原 R 称缺 oracle 导致错误，加入后可保证集合嵌套，证据不足。
[`idr.smk`](../../workflow/rules/idr.smk) 与
[`idr_reproducibility.smk`](../../workflow/rules/idr_reproducibility.smk) 的调用
未使用上游 ChIP v2.2.2 [IDR wrapper](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/v2.2.2/src/encode_task_idr.py)
的 `--peak-list`。仅此源码差异不能证明两次独立拟合等价或集合嵌套；没有真实 IDR
拟合反例/等价实验支持原保证。本轮不修改候选空间或拟合算法。

### S12 — consensus score 与 p/q（撤回实现缺陷）

原 V（部分核验）将 p/q=-1 的下游用途及全支持例中的恒定 score 列为缺陷，
没有直接断言 -1 格式不合法。复核确认 -1 是合法未赋值标记，score 随支持比例
变化，因此撤回实现缺陷。
[`compute_consensus.py`](../../scripts/compute_consensus.py):514 的 score 为
`int(1000 * support_count / n_bioreps + 0.5)`，3/3=1000、2/3=667；signalValue
取支持峰最大信号；:565、:581 的 p/q=-1 表示未赋值，符合
[UCSC 格式](https://genome.ucsc.edu/FAQ/FAQformat.html#format12)。
已有原脚本格式/支持度测试，PR-17 也以合法三 biorep 小输入核对输出（A）。
需要修的是 N3 文档，不新增输出字段或统计。

### S13 — 多候选片段长度（PR-11 已接受）

原 R 的缺陷成立但并非所有输出触发：单候选原可解析，人工列表不足以单独证明真实工具可达。
后续真实 R 多候选 `.cc.qc` 已补强（A）。当前
[`parse_cross_correlation.py`](../../scripts/parse_cross_correlation.py) 的
`_primary_fragment_length` 在 header/headerless 路径按**原序列第一候选**输出标量，
对应工具主峰，不排序/平均/跳过坏首项；空/NA/畸形列表为缺失。NSC/RSC 的通用
标量解析未放宽，quality_flag 仍只依赖 NSC/RSC，原 `.cc.qc` 保留全部候选。
见 [`test_parse_cross_correlation.py`](../../test/test_parse_cross_correlation.py) 与真实 R 测试。
移除 speak 不保证所有输入多峰或必变 NA。

### S14 — 排除窗口（政策待决，未纳入机械修复）

原 V（文本）称 ENCODE 保持 R 底层默认窗口、当前必须随 S4 修改，应纠正。
本地仍传 `-x=-500:15`；ENCODE-DCC ChIP v2.2.2
[WDL](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/v2.2.2/chip.wdl) 设置下界 -500，
[wrapper](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/v2.2.2/src/encode_task_xcor.py)
按 TF/histone 与 read length 推导上界并传 `-x`，不是裸 R 默认值。
本轮在线核对 WDL、核对已归档 wrapper 及其 SHA（在线 wrapper 获取失败）；证据 C。
一组 BAM 在两窗口结果相同不能证明窗口普遍无影响。PR-11 保留窗口，政策继续待决。

### S15 — 线粒体过滤（科学损害证据不足）

原 R 以 ENCODE xcor 的 chrM 过滤为参照，担忧线粒体 reads 对 NSC/RSC/FRiP
的影响，并明确建议先定政策；没有声称所有 assay 必错。本地各 assay 的损害
及统一过滤必要性仍缺证据，继续保留政策决策。
[`common.smk`](../../workflow/rules/common.smk):213 的过滤与
[`qc.smk`](../../workflow/rules/qc.smk):762 的互相关接线可源码核对；特定 ChIP
上游 wrapper 的过滤不代表 ChIP/ATAC/CUT&Tag/MNase 的统一基线。
本轮未重跑原来的宽范围 grep，也不评价保护配置。缺规则本身不能量化科学损害，
不借文档任务新增过滤或改变 MAPQ/primary/proper-pair 政策。

### S16 — SE broad CUT&Tag 首轮延伸（PR-12 已接受）

原 R（分支前提 V）把缺 log input 写成永不重建，撤回该泛化。
合法 SE broad CUT&Tag treatment、extend_reads=auto/yes 可因消费者先跑使用
回退 200，而生产者先跑读取模型长度。下一轮 Snakemake params 检测可触发重跑，
但不能纠正已交付的首次结果。

当前 [`common.smk`](../../workflow/rules/common.smk):373 的 `_bamcoverage_inputs`
增加对应 MACS3 **声明输出目录**依赖，不把未声明的 log 伪装成输出。
真实 metadata dispatch、CUT&Tag 委派与 ChIP helper 保留；PE/no/整数/control/
narrow --nomodel 不因修复新增该依赖。无有效预测仍回退 200。
[`test_bamcoverage_dependencies.py`](../../test/workflow/test_bamcoverage_dependencies.py)
覆盖首次优先级、单独 bigWig、对照和增量 params；独立 171 项通过。
证据 B 使用假 MACS3/bamCoverage，证明调度/argv，不证明科学 bigWig；旧 v2 手写
ChIP dispatch 的探针不等于该真实调用链。真实科学产物/远端 CI 等边界保留。

### S17 — MNase Picard 路径依赖（PR-14 已接受）

原 R 的顺序问题成立，但 `insert_size_metrics` 是**路径或 NA**，不是数值。
Picard 默认 false，MNase 专用 summary 不依赖通用 qc.summary 开关。
当前 [`mnase.smk`](../../workflow/rules/mnase.smk):447 在规范化 Picard 开关为真时
依赖 [`qc.smk`](../../workflow/rules/qc.smk):822 声明的 insert_size 输出。
[`mnase_qc_summary.py`](../../scripts/mnase_qc_summary.py):112 保持存在性判断：
关闭 Picard 但旧文件存在仍报告路径；开启缺参考仍拒绝，不静默跳过。

[`test_mnase_summary.py`](../../test/real_execution/test_mnase_summary.py)
保留原规则/helpers/脚本，用真实 samtools 和合法 PE BAM，Picard 外部工具为替身。
独立 14 项专项＋213 项回归通过：两优先级、单独 summary、失败阻断、关闭/旧文件、
23 列及 read-record 计数保持。A 原脚本/计数与 B 替身调度分开，未验证真实 Picard
分布、科学 bigWig 或完整流程。旧 v2 的手写 helpers 不是原链路证据。

### S18 — 插入长度分布的生产路径（撤回中心指控）

原 R 称完全没有 ATAC 分布生产路径，忽略了 opt-in Picard。
[`qc.smk`](../../workflow/rules/qc.smk):819 声明 `insert_size_histogram.pdf`，
[`targets.smk`](../../workflow/rules/targets.smk):259 的目标覆盖相应 treatment。
MNase summary 已有 sub/mono/di 记录计数，不应将计数称作 PE 片段数，也不等同
自动产生所有希望的比例/图。源码证明路径存在（C），不证明真实数据图形质量。
本轮只澄清默认/开启条件，不引用或评价排除在范围外的本地脚本。

### S19 — MultiQC 项目表（缩窄结论；完整 HTML 未验证）

原 V 由项目 TSV 未接入推出所有相关指标不进 HTML，推论过强。
[`report.smk`](../../workflow/rules/report.smk):180 搜索 active sample 目录及条件性的
cross-correlation summary；[`multiqc_config.yaml`](../../workflow/multiqc_config.yaml):223
有 cross-correlation/MNase 自定义表，没有 project_qc_summary 自定义接线。
这不排除标准模块解析相关原始报告。特定表缺消费者是 C；完整 HTML 内容、尤其
duplication 展示仍依输入和环境，证据不足，不能宣称全不可见或全部已可见。
是否新增汇总展示留决策，不以工作区文件存在证明平台列表/下载。

### S20 — CUT&Tag fraction_lt_120（缩窄结论）

原 R 的“死指标”不准确：[`calc_cuttag_fragment_size.py`](../../scripts/calc_cuttag_fragment_size.py):64
实际写该字段，[`targets.smk`](../../workflow/rules/targets.smk) 的
`_cuttag_fragment_targets` 选择 active CUT&Tag，SE 也有 TSV 路径但指标为 NA。
[`artifact-inventory.yaml`](../architecture/artifact-inventory.yaml) 的该项
`manifest_output_type: null` 不证明平台下载可达；其旧“PE treatment only”注记
也不能覆盖实际 layout/role 门控。本轮不改这一运行时契约文件。
部分汇总未消费字段与工作区存在字段可同时成立（C）；额外展示继续待决。

## 平台条目

### P1 — trusted adapter 授权边界（证据不足，决策保留）

原 V 阅读到 [`command_builder.py`](../../src/encode_pipeline/services/command_builder.py):149、169、275
的 trusted exact-instance 与 capability 分支检查不同；前者允许既有 authority 下
cwd=None。差异本身没有证明不可信输入能够越界。已接受复核中的四项授权/命令
测试通过也不等于穷尽安全证明。需要攻击可达性和权限前提证据，未并入 PR-2/6，
不能按平台目录自动豁免其将来可能影响的执行闭包。

### P2 — HTTP 页大小与 SQLite 前瞻溢出（PR-15 已接受）

原 V 只记无上限，后续 A 使用原应用、临时 SQLite、认证 HTTP 确证 `limit+1`
越界：旧 `2^63−2` 加一仍可绑定；`2^63−1` 及更大可返回 500，并非所有大值都溢出。
当前 [`routes/runs.py`](../../src/encode_pipeline/api/routes/runs.py):781、837
是 default=50、ge=1、le=100，超限按既有 400/API_REQUEST_INVALID 拒绝。
内部仍可查询 101 条前瞻，未改 repository 成 100 上限。
[`test_run_progress_pagination.py`](../../test/api/test_run_progress_pagination.py)
同组 HTTP 断言覆盖内存/SQLite、游标、顺序、stream 隔离、认证、404 和 OpenAPI。
执行闭包变更已同步身份，PR-15 新身份 Gate 仍欠验；原“not a gate item”撤回。

### P3 — 500 响应声明（PR-16 已接受，验证边界保留）

原 R 的声明遗漏经真实临时应用故障注入确证，不能写成正常请求必然失败或已发现脱敏失效。
PR-16 为 30 个 operation 中缺失的 22 项补声明，保留已有 8 项专用结构。
createRun 用 RunResponse；getRun 的显式引用证据失败用 RunResponse，其他全局
异常用 ValidationResponse，当前 :569 声明联合模型。见
[`runs.py`](../../src/encode_pipeline/api/routes/runs.py)、
[`main.py`](../../src/encode_pipeline/api/main.py) 及
[`test_internal_server_error_contract.py`](../../test/api/test_internal_server_error_contract.py)。

独立 155 项 Python、60 项前端通过，声明/生成客户端/资产/身份核对通过（已有记录，
PR-17 未重跑）。成功流式下载仍超时，原因未确定，无证据认定为 PR-16 引入回归；
新身份 Protected Bulk Gate 仍待补验，不能把定向验收写成全验证通过。

### P4 — 终态取消（撤回缺陷指控）

原 R 担忧响应未显式区分“本次取消”与“此前已结束”，没有声称 run.status 丢失。
复核认为返回 canonical 状态的现有幂等契约可接受，不认定代码缺陷。
[`run_cancellation.py`](../../src/encode_pipeline/services/run_cancellation.py):76
及 RunService 返回 canonical run，HTTP body 仍含真实终态。例如 failed run
重复取消可返回 200/ok=true/run.status=failed。复核 A 的临时 HTTP 与 C 源码
支持保留该契约，但不等于证明响应显式标记了本次是否发生状态变更；
不新增状态或 issue 改变生命周期。

### P5 — 两个通知静默边界与 doctor（PR-4 已接受）

原 R 称所有通知失败无痕，错误。普通 transport 失败已有 `terminal_email_failed`
持久事件，原因码 `TERMINAL_EMAIL_DELIVERY_FAILED`。
[`terminal_notifications.py`](../../src/encode_pipeline/services/terminal_notifications.py):226
的 `_record_outcome` **内部**吞掉事件写入异常，不会因此传到外层；
[`runs.py`](../../src/encode_pipeline/services/runs.py):1121 的 notifier 保护是另一个边界。
PR-4 分别添加固定 component/phase/reason_code 日志，保留发送次数、事件、终态及吞吐语义，
另覆盖 worker 成功通知和超时保护，未补发/重试或重复创建普通失败事件。

当前 doctor 没有 `_safely`。实际覆盖
[`DeploymentDoctor._run_probe`](../../src/encode_pipeline/deployment/doctor.py):185
与 [`_run_json_doctor`](../../src/encode_pipeline/cli/local_platform.py):594 的
environment/workflows/recovery 降级；保留 ProbeResult/JSON fallback/退出码。
故障注入验证格式化日志及 extra 不含载荷、环境值、私有路径、原异常或 traceback；
均使用假发送器。证据入口为 roadmap PR-4 与正式通知/doctor 测试，未发送真实邮件。
日志自身仍 best effort，PR-4 精确身份 Gate 欠验独立保留。

### P6 — 读取鉴权（撤回原事实判断）

原 V 只看各路由是否单独声明 Depends，忽略
[`api/routes/__init__.py`](../../src/encode_pipeline/api/routes/__init__.py):22 的父 router
统一注入 require_principal/enforce_csrf。runs/run/events/logs 匿名读取复核 A 均为 401。
原“无认证但与可信 LAN 范围一致”不是当前事实，撤回，不需要放宽或另加认证政策。

## 文档条目与后续补充

### D1 — BigWig 文案（PR-17 修正，待验收）

原正文无 V/R 标签，表中误记 V。源码已具备 FE/ppois bedGraph 到 BigWig 的
条件转换（`qc.smk:421、446` 及 pooled 对应规则；`metadata.smk:496`）。
[QC 文档](../qc-interpretation.md#macs3-signal-tracks) 删除 conversion planned，
写清 signal_tracks、资格及有效 chrom_sizes 条件。文案修订不代表本轮新增转换功能。

### D2 — blacklist 文档（PR-17 修正，待验收）

原正文无 V/R 标签、原总表遗漏。原“Filter BAMs and peaks”不足以证明承诺全局
替换所有产物；实现有过滤副本。已与 S5/S8 同步[配置说明](../configuration.md#qc-block)，
明确消费者及 DAG 选择，不为迎合旧指控修改科学输入。

### 后续 N1 / N2 / N3（不计原 28 条）

- **N1，PR-13 已接受：** 原三列重复峰只能证明计数/舍入代码路径，后续唯一合法
  十列峰重新确证 9998/5000=1.9996 被三位格式化误判。当前
  [`chipseq_idr_summary.py`](../../scripts/chipseq_idr_summary.py):73 与
  [`idr_reproducibility_summary.py`](../../scripts/idr_reproducibility_summary.py):114
  用正分母与整数交叉比较严格 <2；显示仍为 2.000，status 可 pass。
  恰好 2、NA、inf 仍 fail。独立 50 项专项及 151 项相关回归通过，八列/十五列及
  复制字节保持；见 [`test_idr_reproducibility_summary.py`](../../test/scripts/test_idr_reproducibility_summary.py)。
- **N2，政策待决：** [本地 IDR 契约](../idr-contract.md#final-peak-sets) 规定
  true→conservative、pooled→optimal，即使 Nt>Np 也不取较大集合；统一 summary
  复制 true peaks。PR-13 未改变这一兼容选择，不能重复立项或偷偷调整。
- **N3，PR-17 文档修正待验收：** score 是支持比例，signalValue 是最大信号，
  p/q=-1 未赋值；一行 13 列 summary 含 support_distribution，但没有声称的逐峰
  支持字段。保留正确表头，修正[政策 §4.6](../reproducibility-policy.md#46-consensus-output-format-caveats)
  及“两重复等于 intersection”的过强描述，不新增生产字段。

## 历史实验的引用边界

下列是修复前的已接受复核记录摘记，**本轮未复跑旧版本**；命令入口及原始 stdout/
stderr 在对应历史证据目录，不能把旧数字当作当前部署结果：

| 历史记录 | 当时操作与结果 | 可证明/不可证明 |
| :-- | :-- | :-- |
| S1/S2，补强目录 `dag/commands.json` | 原 Snakefile dry-run：合法 mixed 默认 summary 触发 MNase MACS3 guard；缺 FASTQ 先 MissingInput。预置十列峰实际执行 narrow consensus，argparse 缺 final-output 值 | 前者 B 调度，后者 A 原消费者；不是完整科学执行 |
| S4/S13，`root-confirm/qc/commands.json`、`multipeak-commands.json` | 同 BAM 保留 `-x=-500:15`，只去掉 speak；另真实 R 输出 `200,300,1115`，旧 parser 输出片段长度 NA | A：固定零位移覆盖及合法列表可达；不保证所有样本多峰或 NSC/RSC 达标；未传 -p 的机制实验不是生产规则 |
| S6/S7，`root-confirm/qc/commands.json` | `samtools sort -n → fixmate -m → sort → markdup -r` 后跑原指标脚本；合法 PE 的 preseq 加减 -P，以及同库去重前后对照 | A：典型退化与 clipping/反向长度反例；含 -Q/-e/-s 的诊断不冒充生产默认参数 |
| S16/S17，`schedule/commands.json` | 旧隔离调度中 CUT&Tag argv 200/150，下一轮可检测 params；MNase 字段 NA/路径，后来文件出现未自动重跑 | B：原消费者/补强 helpers 与假生产者；不证明科学产物，早期 Nobias v2 简化 helper 边界另见条目 |
| P2/P3/P4/P6，`root-confirm/platform/api-v2-results.json` | 真实临时 SQLite/进程内 HTTP：极大 limit 500、注入异常 500、终态取消保留 failed、四类匿名 401 | A：实际 API 契约；不是生产数据库或正常请求普遍失败 |
| N1，`root-confirm/values/commands.json` | 唯一十列峰 Nt=5000/Np=9998/N1=N2=5000，旧输出 rescue=2.000、status=fail | A：原脚本判级错误；不是实际 IDR 拟合验证；PR-13 已修 |

此处补强目录指 `/tmp/helix-independent-corrections.IEHFno/`。
原上游六项归档的 URL/SHA 可分别核对，但不含 overlap；不得声称原引用全部已归档。
本轮在线确认 UCSC 格式及 ChIP v2.2.2 overlap/WDL；其他固定版本对照按已归档
一手源码核对，在线失败另列 PR-17 证据。没有拼接为跨 assay/版本的“ENCODE 标准”。

## 旧建议顺序的处置与仍开放事项

原“建议修复顺序”**已被 roadmap 取代，不再是执行计划**。其把 S2 称为默认
DAG 失败、把 S14 并入 S4/S13、把若干政策差异直接视为科学缺陷的分组已撤回。
本文件不重排 PR-2～PR-17，也不因 N1 再开修复。原“Verified clean”是当时有限
源码/探针观察，不是一揽子当前合格证；默认 CPM 可配置、脚本 wrapper 存在不等于
规则接线、DAG 覆盖不等于科学执行，不能保留无条件 clean 断言。

S3/S9/S10/S11/S14/N2 等政策、P1 信任边界、S15 损害、S19 完整 HTML、S20
额外展示仍有待决或证据缺口。fastp 固定版本字段兼容问题见
[上游耦合台账 UC-06/07](upstream-coupling-ledger.md)，按用户决定延期，修复归属未定。
PR-2 历史 MACS3 bug #3 当前未复现且无完整部署证明；PR-6 的完整在线 staging、
部署 admission 与科学小样本边界也不因文档完成而消失。

Gate 按最终执行闭包、配置、参考及 artifact/QC 契约影响判断，不按目录或“科学”
标签自动触发/豁免。PR-17 只改普通 Markdown，闭包不变，不重建身份、不新增 Gate。
G0 与 PR-3、PR-4、PR-6 各身份、PR-7、PR-15、PR-16 的欠验分别保留，后续身份
通过不能追认旧版本。PR-17 完成只表示文档交付待独立验收，不表示所有问题清零。

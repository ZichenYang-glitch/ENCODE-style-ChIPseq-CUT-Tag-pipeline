# 科学计算修复 spec / plan

状态：**实施基线。用户已安排 Kimi 通过独立分支/PR 开发，Codex 独立验收；科学代码尚未实现。**

审查基线：`fa0d8401cbc3fecb765fed6d7e01e18bb8fbc96d`。
本文记录当前会话的范围决定、实现设计、执行关卡和验收要求。工作区已有用户暂存改动；本轮只新增本文，不修改科学规则、配置、测试或运行身份。

## 1. 已确认的产品决定

1. 首轮修复默认分析路径和 QC；IDR 方法对齐另立后续工作。
2. **legacy 和 updated 中确认的错误均须修复。**legacy 是兼容配置与既有工作方式的入口，不是冻结错误行为的开关。
3. 保留历史结果、快照、版本与运行记录；在新版本上重跑不承诺得到旧错误行为的相同输出。
4. CUT&Tag 默认重复政策纳入首轮：以保留重复、报告重复率、必要时显式去重为候选，验证后确定 assay 默认，不再直接用 narrow/broad 代理。
5. 普通用户主要提供实验事实；维护者负责默认方法及其适用边界。高级选项只用于有真实需求的覆盖，不为每个工具增设一套相互独立的科学政策。
6. 执行前先讨论并定稿本文；本次新增文档不表示实现或生物学验证已经完成。
7. **blacklist 是可选功能，默认开启，允许显式关闭。**默认开启应过滤分析主链，而不只是生成事后 QC 副本；资源存在性和身份验证以是否开启为条件，不要求所有用户必须使用该功能。
8. **主片段分布 QC 跟随实际分析输入。**blacklist 未启用时使用未做 blacklist 过滤的分析 BAM；启用分析主链过滤时使用过滤后送入下游的 BAM。此处“过滤前”仅指 blacklist 之前，MAPQ、重复和细胞器政策仍独立生效；复杂度和细胞器比例保留各自过滤前的诊断口径。
9. **重复处理按 assay 统一候选默认，两个 profile 一致：ChIP-seq、ATAC-seq、MNase-seq 为 exclude，CUT&Tag 为 retain。**包括 legacy broad ChIP 的默认调整；显式覆盖保留既有作用域和冲突检查。用户已确认采用该候选矩阵，最终默认在真实数据验证后冻结。
10. **CLI 的显式重复政策可覆盖 assay 默认。**例如 ATAC 默认 exclude，但用户明确设 retain 时，分析 BAM、coverage 和所有适用的 peak caller 均保留重复；报告记录此次覆盖。两个 profile 的 CLI/config 入口均须验收。
11. **采用 BAM 层决定重复政策、MACS3 不再隐式二次去重的方案。**所有 MACS3 调用显式使用 `--keep-dup all`，覆盖普通、pooled、consensus 和为 IDR 生成候选峰的调用。IDR 自身方法仍后续处理；公共上游修复可能改变其输入与结果，不能承诺结果不变。该范围由当前会话连续讨论 MACS3 行为后确认，纳入首轮验收。
12. **PE 插入上限先对齐已核对的权威流程，选定本轮默认基线：ChIP-seq/ATAC-seq 2000 bp、CUT&Tag 700 bp、MNase-seq 500 bp。**CLI 可显式覆盖。后续验证实现、边界和代表性结果，不以寻找每个文库的最优上限作为实施前置条件；选定基线不等于已完成科学验收。
13. **细胞器采用 assay 默认，两个 profile 一致：ATAC/CUT&Tag 默认排除，ChIP/MNase 默认保留并报告。**用户可显式覆盖；只处理参考中已确认的细胞器名称，不猜测任意 accession。
14. **cross-correlation 的未知靶标及非 ChIP assay 展示诊断值和适用说明，不直接套 ChIP 阈值判合格。**该 QC 维持默认关闭；不能把无法判级写成通过。

## 2. 首轮结果目标及边界

| 编号 | 目标 | 范围 |
| --- | --- | --- |
| S1 | peak caller 忠实消费上游选择的分子集合 | 普通、pooled、consensus、IDR-ready MACS3；首轮统一重复参数，IDR 方法后续 |
| S2 | PE 片段上限明确、有依据、可追溯 | ChIP-seq、ATAC-seq、CUT&Tag、MNase-seq；SE 不受此限制 |
| S3 | CUT&Tag 重复处理按 assay 政策生效 | 两个 profile 的 BAM、coverage、MACS3、SEACR、pooling |
| S4 | 默认 blacklist 生效，参考过滤不因命名静默失效 | 可显式关闭的 blacklist 分析主链、细胞器参考绑定、实际 contig 验证、分析输入与 QC 交付 |
| S5 | fragment-size 和 MNase QC 统计阶段明确 | 合法负 TLEN、零 TLEN、重复状态、选定窗口及未分类计数 |
| S6 | cross-correlation 不因错误窗口产生无条件合格 | 现有 opt-in QC 的输入、排除窗口、结果解释与报告 |
| S7 | 文档、schema、解析、执行与结果契约一致 | 所有被 S1–S6 直接影响的入口与消费者 |

首轮不改 IDR oracle、optimal 集合选择、rescue/self-consistency 判级政策及 IDR 专用输出结构。也不引入 footprint/cut-site 产品、spike-in 功能、平台重构、数据库迁移或完整文库元数据系统。

保留已成立的合理选择：FRiP 的 read-record any-overlap 定义、当前 fragment 模型、MNase 非穷尽窗口、MACS3 hs/mm 常量。各自的输入阶段和含义须写清楚。

## 3. 建议的科学与配置契约

以下把已确认决定转为实现契约。字段设计和执行顺序由维护者负责，按分批 PR 实施；数据验收关卡不能被文档定稿替代。具体交接见 [Kimi 开发与 Codex 验收说明](scientific-computation-kimi-handoff.md)。

### 3.0 分析阶段与消费者路由

两个 profile 共用“有效分析政策 → 分析 BAM → 下游消费者”的路由。保留既有 raw/final 文件的来源语义，并对两个 profile 都提供明确的 analysis BAM；不把所有历史 final.bam 重新定义为同一重复状态。

| 阶段 | 含义 | 主要消费者 |
| --- | --- | --- |
| 去重前 MAPQ 数据 | 现有 MAPQ/基础比对过滤后的数据，尚未执行重复排除、细胞器排除和 blacklist | NRF/PBC、细胞器诊断及适用的原始文库诊断 |
| 重复标记/处理 | 后端按解析后的政策标记或处理重复；保留原始 metrics | 原始重复诊断、分析视图生成 |
| 分析 BAM | 执行重复、细胞器和可选 blacklist 主链政策后的实际分析输入 | MACS3、SEACR、coverage、FRiP、主 fragment-size、MNase 分箱及相关分析 |
| 最终峰 | 对分析 BAM 调用的峰，按启用政策进行最终峰级 blacklist 过滤 | 首轮普通/pooled/consensus 交付与对应计数/FRiP |

PE 分析视图须保持完整模板：proper-pair、同 contig、成对 R1/R2、互相一致的 mate 坐标及有效非零 TLEN；MAPQ 或其他过滤只留下单侧的模板整对排除并计数。复用现有模板处理边界，补足当前“同 QNAME 恰好两条”不等于完整 mate 验证的缺口；不引入全新比对引擎。legacy 的原始中间 BAM 可保留既有比对记录，但其下游分析视图也执行这一契约，因此 BAM/coverage 的有效集合可能变化，须列入兼容说明。

各过滤阶段以模板为单位执行，不能依赖 blacklist 开启后才修复孤儿；关闭 blacklist 也须满足模板闭合。外部 control 在 admission 边界验证相同的布局与模板前提；无法确认其原始文库重复处理时如实记录，不声称已恢复外部 BAM 中删除的信息。

既有 read-level `.CPM.bw` 与 fragment-level `.fragment.CPM.bw` 继续作为不同产品，明确前者的计数单位；本轮不以修文档为由改变前者算法。IDR 专用交付/计数的历史差异仍属第 7 节，不把首轮普通峰修复宣传为 IDR 方法对齐。

### 3.1 重复处理只有一个分析政策

- BAM 层负责执行已解析的保留/排除政策；Picard 标记数据保留用于复杂度和重复率诊断。
- 新修正版中所有 MACS3 调用显式使用 `--keep-dup all`，避免再次按坐标裁剪上游分子集合。普通、pooled、consensus 和 IDR-ready 必须共用相同的重复参数构造逻辑；原有正常/宽峰/IDR 阈值仍各自保留。
- 不新增一组彼此独立的 `macs3.keep_dup`、`idr_macs3.keep_dup`、`consensus.keep_dup` 用户开关。
- 对现有 extra_args 中的 `--keep-dup`：同义 `all` 可规范化去重；冲突设置明确拒绝并给出迁移说明，不能靠命令行参数先后顺序决定结果。
- 两个 profile 的显式 `remove_dup=yes/no` 和 `duplicate_analysis=exclude/retain` 都应产生明确、可验证的结果。冲突继续拒绝；auto 的解析应与每个 sample 的实际政策一致。
- CUT&Tag auto 的候选是 retain，narrow/broad 不决定是否排除重复。显式 exclude 仍可使用。IgG 不等同所有 `role=control`；首轮不据此自动把全部 control 去重，也不在缺少元数据时猜测控制类型。
- 同一文库跨 lane 与独立文库须区分。首轮不在所有 merged/pooled BAM 上追加坐标去重；当前 sample/LB 表达能力的限制应记录。完善 library_id 是另一个有独立消费者的需求，不夹带完整样本模型改造。

已确认的候选默认矩阵：

| assay | auto 的有效分析政策 | narrow/broad 是否改变此政策 |
| --- | --- | --- |
| ChIP-seq | exclude | 否；legacy broad 也调整为 exclude |
| ATAC-seq | exclude | 当前支持 narrow |
| CUT&Tag | retain | 否 |
| MNase-seq | exclude | 不走 peak-mode 找峰路径 |

缺省/auto 必须保留到每个样本的 assay 政策解析阶段；不能先规范化为全局 exclude，导致混合 ChIP/CUT&Tag 项目丢失不同默认。显式 exclude/retain 则仍按用户配置的实际作用域生效。执行快照与报告记录每样本解析结果；测试覆盖混合 assay、两个 profile、显式覆盖及二次规范化幂等。

CLI 覆盖的验收例是 ATAC 显式 retain。建议通过现有输入配置提供 `duplicate_analysis: retain`，并兼容 `remove_dup: "no"`；两者同时存在须一致。新版本中一致的显式选择覆盖 auto，冲突选择明确报错。当前 CLI 主要经配置文件传递此政策，不能在实现独立命令行 flag 之前宣称已有 `--duplicate-analysis` 选项。Snakemake `--config` 的覆盖需要经过相同 validator；仅修改 MACS3 extra_args 不算完成 BAM/coverage 的重复政策覆盖。

`retain` 覆盖只改变重复处理，其他明确启用的 MAPQ、proper-pair、organelle 和 blacklist 处理继续生效。现有设置是 run 配置作用域，不能把全局 retain 描述为“只影响该 run 中的 ATAC”；若需要同 run 中每 assay/每样本独立覆盖，另设计有明确消费者的作用域，不在当前示例中隐含承诺。

为 IDR 生成候选峰的 MACS3 调用也存在相同的重复处理。本轮同步上述公共重复参数，使普通/pooled/consensus/IDR-ready 均忠实消费上游 BAM；不改 IDR 的拟合、oracle、optimal 或判级方法。此项与“重做 IDR 方法”不同，但可能改变 IDR 的输入峰及最终结果，须明确告知并验证。这里确认的是 spec 范围，实施仍等待本文整体定稿。

### 3.2 PE 插入长度采用维护的 assay 默认

- 在现有 Bowtie2 tool_parameters 边界增加一个结构化上限设置，具体键名建议 `max_insert`，支持 assay 默认与正整数覆盖。
- 解析为具体每样本数值后，显式渲染 `-X`，记录到执行身份与用户可读的分析说明。SE 不渲染该参数。
- 用户已选定本轮默认基线：ChIP-seq/ATAC-seq 2000，CUT&Tag 700，MNase-seq 500。依据是权威流程的现有设置和当前支持窗口，不能仅凭合成数据宣布适用于所有文库。
- MIN/MAX 组合、过长 reads 与 Bowtie2 实际边界语义一并验证；首轮不顺带更换 local/end-to-end、dovetail、MAPQ 等政策。
- 对 extra_args 中 `-X`/`--maxins` 及等号形式制定一致迁移规则：未提供结构化值时可规范化合法旧覆盖；显式冲突时拒绝。历史覆盖不能静默消失或被新默认压过。
- 500/501/600 bp 只是证明已有边界问题；新候选还须测试其上限两侧及真实配对行为。

| assay | 本轮默认 max_insert | 依据与验证边界 |
| --- | ---: | --- |
| ChIP-seq | 2000 | ENCODE ChIP v2.2.2 PE wrapper；验证长片段和正常短片段，不声明全部 ChIP 文库必须使用该值 |
| ATAC-seq | 2000 | ENCODE ATAC PE wrapper；关注 nucleosomal 长片段、配对率及片段分布变化 |
| CUT&Tag | 700 | 方法团队教程的候选范围；验证不同目标及 retain/exclude 比较中的独立影响 |
| MNase-seq | 500 | 保持现有配对上限；nucMACC 发表流程有相同设置；不改变 sub/mono/di 选定窗口，检查当前支持范围与长尾诊断 |

MNase 参数需在产品说明中单独区分：约 147 bp 是核小体核心包裹的 DNA 长度；130–200 bp 是常见单核小体信号选择范围。仓库当前 dyad 的 `--MNase` 窗口是 130–200，mono occupancy/mono BAM 是 140–200，di 窗口是 300–400。表中的 500 仅是上游配对上限，不能描述为“单核小体取到 500 bp”。配对时若直接把上限设为 200，会先排除当前 di 产品需要的模板；若未来提供 mono-only 运行策略，应单独设计并说明输出范围，不能在当前多窗口流程里隐式收紧。

本轮暂保留 mono=140–200、dyad=130–200 的既有区别并明确报告，不因用户提出常见 130–200 窗口就直接更改 mono 定义。若需要统一为 130–200，应与 sub 窗口边界和指标解释一起决策，并单独比较 130–139 bp 模板的归属及信号变化。

上述数值已作为本轮默认基线选定，尚未执行本轮实现验证。`max_insert` 是配对约束，不等于生物片段分类窗口；提高它也不意味着所有新配对都是可靠生物信号。验收需检查新增配对的质量和下游影响，保留与上限一致的 QC 说明。

用户提出过大上限的风险，因此真实数据验收必须同时评估收益与代价：

- 上限变大扩大合法配对的搜索范围，可能改变候选配对、比对位置、MAPQ 和 concordant/discordant 分类；它不会把 read 延长，也不等同放宽单条 read 的碱基错配评分。
- `proper_pair` 是给定配对参数下的判定，不能把上调 -X 后 proper 数量增加直接视为准确率提高；高 MAPQ 也不能单独证明新增长片段全部正确。
- 搜索窗口扩大有运行时间成本，记录相同线程、参考、输入和其他参数下的 wall time 与资源使用；不能假定成本与上限线性变化。
- 真实数据至少比较旧默认与本轮权威基线：ChIP/ATAC 为 500 对 2000，CUT&Tag 为 500 对 700；MNase 保持 500 并验证当前窗口与诊断。中间值或更宽值比较可用于解释异常，不作为为每个文库寻优的必经步骤。固定输入，记录新增长片段数量、质量分布、坐标变化、原高置信短片段的稳定性及信号/背景和重复间一致性。
- 对本流程 BAMPE 和 fragment coverage，新增长片段也会改变信号形状及背景，须实际比较峰/coverage；保持 MAPQ 等其他阈值不变，避免通过同时加严其他过滤掩盖上限变化的效应。
- 默认基线与对应上游版本一起记录，报告验证的支持文库范围和资源代价；不以越大越好、保留越多越好或单一峰数作为正确性标准。若验证发现不适用的明确条件，记录例外并讨论调整，不自动静默调参。

CLI 配置表达为 `tool_parameters.bowtie2.max_insert: <正整数或 auto>`，缺省为 auto；拒绝布尔值、零、负数、小数和空字符串。auto 由维护的 assay 默认解析。现有配置作用域是全 run，不能把该整数误报成只覆盖某个 assay；每样本有效数值在快照/报告中可见。针对 SE，校验配置合法性但不渲染 -X，并记录“不适用”。

旧 extra_args 使用 shell token 解析而非子串匹配，覆盖 `-X 700`、`-X700`、`--maxins 700`、`--maxins=700`。缺省/auto 与一个合法旧值并存时采用旧显式值；显式整数与旧值相同则规范化为一个参数，不同则拒绝。多个互相冲突的旧值也拒绝；解析成功后移除旧 token，保证最终命令只有一个 -X。已有 -I/--minins 覆盖保留并校验有效 min ≤ max，不另加新的最小插入默认。

### 3.3 fragment QC 采用明确阶段，保留诊断指标

- 主 fragment-size QC 描述最终分析使用的片段，消费与分析相同的 BAM；blacklist off 时为未做 blacklist 过滤的分析 BAM，analysis 时为过滤后的分析 BAM；是否包含 duplicate 由该样本的分析政策决定，不再无条件额外排除 0x400。copies 仅生成额外副本，其主指标仍跟随未做 blacklist 过滤的分析输入。
- 只计一次完整有效 PE template。合法负 TLEN 必须保留，以 `abs(TLEN)` 求长度；TLEN=0、缺失 mate、跨 contig 或不一致模板单独计数或明确报错，不能落入 0 bp 生物片段桶。
- 分析视图生成时已排除的无效模板记入阶段计数；最终 fragment-size 表的无效输入计数另列，不能因为最终 BAM 已清理就让上游 TLEN=0/孤儿数量消失。
- 重复率和复杂度继续取适当的去重前阶段，不能为了统一输入而改成 analysis BAM。原始文库诊断与最终分析分布是不同指标，应标明阶段。
- QC 记录至少包含统计阶段、重复/细胞器/blacklist 政策、解析后的片段上限及有效/无效模板数。新增字段通过版本化 QC 契约提供，不静默扩展严格固定 header 的旧解析器。
- MNase 保留选定窗口；增加 total/invalid/classified/unclassified 的模板计数。可配置窗口可能重叠，未分类数以窗口并集计算，不能简单 total 减各类之和。
- read count 与 fragment count 分列，避免在 legacy 可能存在不完整 pair 时直接除以二。

主片段分布与诊断指标的阶段原则已获用户确认；具体统计阶段矩阵如下。这里的“最终分析 BAM”表示首轮统一路由后的实际分析输入，不依赖 legacy/updated 的历史文件名猜测状态。

| 指标/产品 | 统计输入阶段 | 必须说明的含义 |
| --- | --- | --- |
| 主 fragment-size QC | 最终分析 BAM 中的有效 PE 模板 | 长度取 abs(TLEN)，每模板一次；重复保留与否服从分析政策 |
| percent_duplication / estimated_library_size | 重复标记及其原始统计 | 描述被标记前的文库分子分布；不能从已排除重复的 analysis BAM 反推 |
| NRF/PBC | 现有去重前 MAPQ 输入 | 首轮不改其片段键定义；明确其细胞器/blacklist 阶段，避免假称与分析片段分母一致 |
| mito_fraction | 去重前 MAPQ 输入，位于重复排除、细胞器排除和 blacklist 之前 | mapped alignment 比例，含潜在重复；两个 profile 均不依赖 final.bam 的不同去重状态；命名匹配与是否应用排除另列 |
| blacklist 丢弃数量 | 进入 blacklist 模板过滤器前后 | PE 同时报告 template 与 alignment 计数；不把此前 duplicate/organelle 丢弃数混入 |
| FRiP | 最终分析 BAM 与该政策下的对应峰 | read-record、any-overlap；记录分母和峰阶段，不改为 fragment FRiP |
| MNase classified/unclassified | 最终分析 BAM 中有效 PE 模板 | invalid 单列；按窗口并集统计，保留各类实际 read/fragment 数 |

报告分别记录 `input_stage` 和 `blacklist_applied`，不能仅凭“final”或“过滤前”推断是否做过 blacklist。关闭 blacklist 时，主指标没有额外的“过滤后”副本要求；开启时可报告过滤前后计数，但不要求每项 QC 都重算两份。

`retain` 场景下主 fragment-size 含标记重复是符合分析政策的行为；需要修的是隐含阶段变化和无效模板计数，不能把 0x400 存在本身一律判错。

### 3.4 参考过滤：默认 blacklist、显式身份与过滤前指标

- 已确认：两个 profile 的 organelle_mode 缺省为 auto，ATAC/CUT&Tag 解析为 exclude，ChIP/MNase 解析为 report；显式 keep/report/exclude 优先。keep 不执行排除也不要求细胞器诊断资源，报告记录此次选择。
- 优先使用参考绑定中的 organelle_chroms，并核对实际 BAM/reference contig；fallback 仅作为已声明别名，不能从任意 `NC_*` 猜测线粒体。
- 显式 organelle_chroms 中每个名称均须存在，部分拼错也拒绝；缺省别名候选只取实际匹配项，不能要求 chrM/MT/chrMT 三者同时存在。最终记录使用显式清单还是兼容别名匹配。
- exclude 且无任何匹配时须在产出分析结果前失败，提示明确配置参考身份或显式选择 keep。参考确实不含细胞器时，可显式 keep 并记录“不执行该项过滤”；不把缺少名称或空列表自动解释为“参考生物学上没有细胞器”。这避免重释现有 reference binding 中缺省空列表的含义，也不需要本轮增加一套参考元数据模型。
- report 模式匹配失败时，记录 `unresolved_contigs`、比例 NA 和可操作提示；不能输出伪造的 0。匹配成功但 mapped alignment 分母为零时，记录 `no_mapped_reads` 和 NA。名称匹配且分母大于零、细胞器 reads 为零，才是合法的零比例。
- mito_fraction 使用 3.3 定义的去重前 MAPQ 数据；标明包含潜在重复，接入 adapter QC/artifact 与 MultiQC。不能改用排除后的 BAM 造成恒零。
- 新修正版未指定 blacklist 政策时，规范化为 `blacklist_mode=analysis`，在两个 profile 的分析主链生效；用户可显式设为 off。普通新建表单与发行模板提供相同默认。
- 开启时在运行前校验资源存在、有效区间与参考 contig 身份；缺失、损坏或完全不匹配时明确失败，提示补齐资源或显式关闭，不静默退化为 off。部分 contig 不匹配的判断需结合 assembly/reference 身份，不能仅靠存在一个同名 chr1 宣布参考正确。
- 关闭后不要求 blacklist 资源，也不调度其过滤依赖；结果说明中明确记录未过滤。
- 继续精确区分 off/copies/analysis：现有显式 copies/off 以及旧 alias 的明确选择保留原意，copies 不能被描述为全分析过滤。发行模板可迁移到 canonical analysis，但不修改用户已有配置或静默重释旧 alias。
- 当前 legacy 的 analysis-mode 验证禁令、updated-only 派生规则及消费者路由都需调整，才能使修复真正覆盖 legacy；不自动把整个 legacy 切换为 updated 或连带改变其他科学政策。
- PE 使用完整模板过滤，首轮交付范围内的最终峰再做区间过滤；MACS3、coverage、SEACR、FRiP 与汇总读取的阶段须有明确矩阵。原始重复诊断和细胞器比例继续使用各自合适的过滤前输入。关闭过滤时不改变 FRiP 的 read-record any-overlap 定义。

资源校验采用可执行的边界：BED 坐标为整数、满足 0 ≤ start < end ≤ contig 长度；非注释记录不能为空；BED 中每个 contig 都必须属于所用参考，未知名称或越界区间直接拒绝。BED 没覆盖某个参考 contig 并非错误，不要求 blacklist 覆盖全部染色体。FASTA/FAI、chrom_sizes、BAM @SQ 对同名 contig 的长度必须一致；提供的 assembly/binding 身份必须相容。本地未绑定来源的 BED 只可标注为“坐标兼容、assembly 来源未独立核实”，不能仅凭同名 contig 宣称 assembly 已认证。存在显式 blacklist_mode 与旧 qc.blacklist_filter 时，按现有契约检查一致性；不一致则报错，不能让 alias 默默压过新默认。

### 3.5 cross-correlation：计算修复与解释边界

- 维持 opt-in，不提高为默认必跑 QC。
- ChIP 根据实际输入长度和可识别的 TF/histone 类型计算排除上界；TF 为 max(read_len+10,50)，histone 为 max(read_len+10,100)，下界为 -500。不能仅凭 broad/narrow 判定目标类型。
- 首轮使用最终分析 BAM 转换的 tagAlign：PE 仅取有效模板的 R1，SE 取有效主比对；不自动裁剪至 50 bp，不另改重复政策。该 QC 输入仍随主链 blacklist 开关变化，阶段记录为 `analysis_r1` 或 `analysis_se`。这明确区别于 ENCODE 的裁剪、去重前输入，不承诺 NSC/RSC 与其数值直接可比。
- 增加一个有明确消费者的高级选项 `qc.cross_correlation_target_type: auto|tf|histone`，缺省 auto。auto 复用既有明确的 histone 靶标分类；其他名称保留 unknown，不因“不在 histone 名单”就判为 TF。显式 tf/histone 仅用于 ChIP；非 ChIP 始终诊断展示。此设置作用域为 run，不宣称解决混合靶标的逐样本覆盖。
- 长度估计与锁定 run_spp.R 对同一 tagAlign 的实际估计一致，并记录其方法、估计值与输入长度分布；当前锁定脚本为前 500 条有效记录跨度的中位数再按 R 规则取整，不能换成原始 FASTQ 长度后假定等价。变长输入增加 `variable_read_length` 状态，首轮只作诊断。unknown/非 ChIP 的诊断窗口使用 max(read_len+10,100)，说明该窗口用于避开 phantom，不作为该 assay 的最优模型。
- phantom 邻域被选中时输出解释状态，不能仅因 RSC≈1 就无条件标 ok。该状态表示解释存疑，不自动宣称实验失败。
- ATAC/CUT&Tag/MNase 不直接按 ChIP 阈值作统一生物学合格判定；展示诊断值与适用说明。
- 不把“省略 -x 使用底层默认窗口”视为与 ENCODE 的 -500 下界等价。保留真实 R 选峰与生产 parser 的联动测试。

parser 分开记录执行状态、解释状态和数值阈值提示。无数据/工具失败/格式错误不写成生物学失败；未知类型、非 ChIP、变长输入或 phantom_selected 时不输出无条件 ok。固定长度且类型明确的 ChIP 可以保留带有输入政策说明的本地 NSC/RSC 阈值提示，名称不得暗示 ENCODE 认证。phantom_selected 应依据实际 phantom 峰位置、输入长度邻域和对应主峰候选联合判定，不能仅因 RSC 接近 1 就判错；多候选长度与分数须成对解析，覆盖并列峰和缺失值。

### 3.6 报告和解析器必须解释同一份政策

- 每个样本报告有效 assay/profile、重复政策、片段上限、细胞器匹配与排除状态、blacklist 状态、参考资源身份和关键阶段计数。数值采用解析后实际生效值，不只展示用户输入的 auto 或缺省键。
- 用户可读的报告不泄露工作区绝对路径、环境值或私有配置；参考名称、绑定版本和允许公开的摘要用于追溯。
- 固定 TSV header 的新指标采用明确的 v2 QC 契约；adapter 同时精确识别旧 v1 与新 v2，新增指标进入既有通用 QC 模型，平台中立层不判断 assay。
- 新 writer、artifact inventory、adapter extraction、QC summary 和 MultiQC 在同一实施单元更新；验收包含用户实际可见的列表、数值和下载，不能仅证明目录里生成了 TSV。
- 缺测、未启用、不适用、解析失败和合法零值必须区分。未获方法学验证的阈值不输出生物学 pass/fail；尤其不能用统一 duplicate-rate 阈值给全部 CUT&Tag 文库判级。

建议的最小输出契约如下；实现时在 writer/reader 中集中定义精确 header，禁止“接受任意额外列”的宽松解析。

| 输出 | 本轮契约 |
| --- | --- |
| sample qc_summary v2 | 以 `schema_version` 开头，随后保留原 v1 列序，再追加本轮分析模板、细胞器及 cross-correlation 指标与状态列；分析模板字段使用 `analysis_` 前缀，避免覆盖原复杂度指标 `total_fragments` |
| mnase_qc_summary v2 | 以版本列和原 v1 列序为基础，追加有效/无效模板数、classified/unclassified、各窗口 fragment 数、fraction_classified 及 TLEN=0 计数；原 sub/mono/di read 数继续保留 |
| pooled_qc_summary | 首轮保留现有固定列；其中 peak count 必须对应实际交付峰，不引入 IDR 专用计数 |
| fragment-size / mito / cross-correlation 明细 | 各自带版本和明确的阶段/状态；adapter 从固定格式提取本轮需要的值进入 sample v2，不把任意明细表当成通用可解析 TSV |
| 每样本分析政策说明 | adapter-owned、版本化 JSON，固定字段记录实际政策、输入阶段、contig 匹配状态和各阶段计数；作为已登记 artifact 可下载，并由 QC 报告/MultiQC 展示相关摘要 |

sample v2 至少提供 `analysis_valid_fragments`、`analysis_invalid_templates`、`organelle_reads`、`organelle_total_mapped_reads`、`organelle_fraction`、`organelle_status`，以及 cross-correlation 的 NSC/RSC 与解释状态。未启用的 optional QC 写明确状态和 NA，不伪造数值。分析政策说明至少记录 `schema_version`、`scientific_policy_revision`、sample、assay、profile、resolved_duplicate_analysis、resolved_max_insert、organelle_mode、organelle_status、matched_organelle_contigs、blacklist_mode、blacklist_applied、input_stages 和 stage_counts；参考摘要沿用现有允许公开的身份字段。

adapter 精确支持旧 v1 与新 v2；不修改历史 artifact，不要求历史 QC 重新计算。未知版本、重复列、错误类型或未定义状态明确拒绝该 QC 源并报告提取问题，不能把提取失败当作零值。新指标进入既有 metric 模型；政策说明保持 adapter 所有，不把 assay 私有字段塞入平台中立生命周期模型。

计数必须为非负整数，fraction 必须在 [0,1]；NA 通过状态解释。NSC/RSC 按各自数值语义验证有限值，尤其不能把有限负 RSC 截成 0 或因沿用旧计数解析器而丢弃整个 QC；NaN/Inf 和分母为零的未定义结果以缺测状态处理。

### 3.7 配置解析与 CLI 覆盖

配置输入版本继续接受现有 1.0.0/2.0.0；不增加第三个 analysis_profile。科学修复版本通过第 4 节内部政策身份表达，不能由用户选择旧错误行为。authoring schema 因新增字段和默认语义改变，从当前 2.0.0 升至 3.0.0；它与 YAML 的 config_schema_version 是不同契约，不捆绑改号。

| 配置项 | 缺省/auto | 显式输入与迁移 |
| --- | --- | --- |
| duplicate_analysis | 新增 auto，按样本 assay 延后解析 | exclude/retain 优先；旧 remove_dup yes/no 转为相同政策；两者冲突拒绝，remove_dup auto 不压过显式值 |
| tool_parameters.bowtie2.max_insert | auto，按 3.2 矩阵 | 整数覆盖；旧 extra_args 按 3.2 规范化 |
| blacklist_mode | 没有 canonical/alias 显式选择时为 analysis | 显式 off/copies/analysis 保留；仅旧 alias 时保留其明确意图，迁移说明提示 copies 的有限范围 |
| organelle_mode | auto，按 1.13 解析 | keep/report/exclude 覆盖 |
| qc.cross_correlation_target_type | auto | tf/histone 是有输入范围说明的高级覆盖 |

先读取用户显式选择，再迁移 alias，最后补缺省；不能把 schema 自动填入的默认当成用户显式输入。两次 validator/authoring/materializer 往返后，canonical 配置和每样本有效政策必须稳定。保持未解析 auto 与 resolved sample policy 分开，避免混合 assay 被全局默认抹平。

原有配置版本对 profile/control 等未在本文变更的默认继续按既有兼容规则解释，规范化快照明确写出解析后的 profile；本轮修正 schema 文案，不能再写“legacy preserves existing outputs”。authoring 对当前 execution gate 不支持的组合给出准确条件约束，例如 updated + samtools 后端；本轮不为了消除文档矛盾就顺带实现新 backend 或新 control 策略。

CLI 使用现有配置文件及受 validator 约束的覆盖入口。验收必须从用户实际执行入口出发，证明 ATAC retain 同时改变 BAM/coverage/peaks，而不是只测 Python helper。未新增独立 CLI flag 前，文档示例只展示真实存在的配置入口。

**命令行注入覆盖默认是正式验收要求。**同一个 canonical 配置键的优先级为：本次 CLI 显式覆盖 > 运行配置文件显式值 > assay/workflow 默认；不能在校验或规则参数生成时重新填回默认，把覆盖吞掉。Snakemake 的 `--config KEY=VALUE` 是现有入口，合并后的 config 已进入 workflow 的 validate_config。顶层覆盖示例：

```bash
snakemake -s workflow/Snakefile --configfile run.yaml \
  --config mapq=20
```

本轮实现后须另外验收 `--config duplicate_analysis=retain` 在 ATAC 上贯穿全链；若旧配置同时显式设置了相反的 remove_dup，按上述 alias 冲突规则明确提示统一选择。嵌套工具参数可通过运行 YAML 提供；不宣称未经验证的点号键语法会更新嵌套 mapping。有效配置须保留未覆盖的兄弟字段，并把覆盖后的值写入身份与分析说明。新增 override 应使用与之匹配的新运行目录，不重写旧运行身份。

CLI 注入仅改变经校验的公开配置；同一样本的分析重复政策仍只在 BAM 层执行，MACS3 keep-dup all 作为该契约的下游实现保持一致。dry-run 和实际运行须解析出相同覆盖值，非法输入明确拒绝；不能把任意未登记的工具参数当成已经支持的结构化配置。

平台科学参数沿用已存在的配置链：用户通过 schema 表单或 YAML 编辑配置（API 接收配置对象），校验并冻结为快照后，由 adapter 生成工作区 config.yaml，再以 --configfile 执行。用户不必手写或上传 YAML 文件；表单覆盖同样是有效配置输入。平台不提供绕过已验证快照的临时 --config/任意 shell 参数注入。修改科学参数须重新校验并生成新运行快照；cores 等已有执行选项继续走独立、受校验的 options 契约，不混同科学配置。

CLI 需要修正/声明的范围到此为止，不另外设计一套科学参数体系：

1. **入口与优先级**：文档区分 Snakemake 直接运行与 helixweave 平台管理 CLI；不把 Snakemake 的 --config 冒充为所有子命令通用选项。给出真实可执行的顶层覆盖示例；复杂嵌套配置推荐运行 YAML，并说明参数作用于整个 run。
2. **覆盖确实生效**：合并完成后再校验和解析 assay 默认，最终工具命令与所有下游消费者使用同一政策。有效值与原配置矛盾时不能静默回退；alias/结构化参数冲突按本文规则报错。
3. **校验对象一致**：现有独立 `validate --config FILE` 只验证该文件，不能称为验证了之后额外注入的 --config 值。带覆盖运行应通过 workflow 的合并配置校验与 dry-run 检查；单独校验复杂配置时使用实际运行的最终 YAML。首轮不因此新建第二套覆盖语法。
4. **查看与追溯**：运行前显示用户可读的有效 assay 政策摘要（包括重复、max_insert、organelle、blacklist）；实际运行保存合并、规范化后的有效配置和每样本政策。dry-run 可显示但不写入运行目录。用户可见摘要不打印私有路径或完整敏感配置；本地运行记录按既有私有文件边界保存，不作为未经清洗的公开 artifact。
5. **回归与升级说明**：验收 YAML 显式值、CLI 覆盖同名键、两者得到相同有效值时的结果一致性，以及旧 alias 冲突、非法值、嵌套兄弟字段保留和新旧 outdir 隔离。区分“当前入口已支持”与“本轮才修复的跨 profile/全链行为”，不能先发布尚未实现的 flag 或保证。

## 4. 兼容与执行身份

- 不迁移、覆盖或重新解释既有结果和不可变快照；历史原版本行为依赖原运行时与工作流版本，而不是新版本继续保留错误。
- legacy 与 updated 的已确认修复都进入新发行的工作流身份。采用内部、非用户可覆写的 `scientific_policy_revision`；首轮修复使用新的修订号，解析后的 per-sample 政策和工具参数进入既有运行身份。工作流打包身份仍覆盖真实源码，不以政策修订号代替代码身份。
- 复用 `.helixweave/analysis_identity.json`，升级为带 marker schema version、profile、policy revision 和有效配置身份的固定格式，两个 profile 都检查。保留原有 symlink/非法文件拒绝，不自动把旧 marker 当成新格式。
- 新代码不得向已含旧科学输出的目录继续写入修正版产物。CLI 与平台均验证，不依赖普通用户自行记住版本差异。
- 更新表单默认、validator 默认、序列化/渲染、导入配置、CLI 和平台入口，并测试两次规范化不再变化。不能只修改示例 YAML，让常规新建运行继续绕过修复。
- 科学政策属于 ENCODE adapter/workflow；不向平台中立 services/API 引入 assay 私有字段或路径知识。
- 只在实际 OpenAPI 变化时重生成并检查规范和前端客户端；原生 authoring schema/QC 契约变化也须覆盖其消费者。不要手工编辑 generated 客户端。
- 身份清单、运行时打包和适用 Gate 按最终影响逐项判断。ENCODE 修改不自动触发 Bulk 全量科学执行；也不能因为本轮没有改 nf-core 源码就忽略共享闭包影响。

输出目录的明确处理规则：

| 目录/marker 状态 | 行为 |
| --- | --- |
| 新目录，或仅含允许的调度元数据且没有旧科学产物 | 正常创建本次 marker，之后才执行科学任务 |
| 新格式 marker 与 profile、政策修订及有效身份一致 | 允许同一运行继续；仍遵守现有输入/运行身份约束 |
| 已有科学产物但没有 marker | 拒绝在原目录执行，要求指定新 outdir；不删除已有内容 |
| 旧格式 marker，或任一身份不一致 | 拒绝继续写入，保留原始结果；提供新目录重跑说明 |
| marker 损坏、非普通文件或 symlink | 明确拒绝，不尝试自动修复或覆盖 |

CLI dry-run 应能发现不匹配，但 dry-run 不创建 marker 或修改科学产物。平台继续通过 adapter 验证与现有不可变快照/工作流身份传递此约束；不把路径判断逻辑复制到通用 API/services。两个 profile 旧目录的拒绝是本轮有意的升级行为，不提供强制覆盖旧科学结果的便捷开关。

## 5. 建议实施单元与验收

这里的 M0–M6 仅为本文工作单元，不重排或复用 roadmap 中已有 PR 编号。

| 单元 | 消费者结果 | 核心验收 |
| --- | --- | --- |
| M0 契约与验证准备 | 用户与维护者知道支持范围、默认政策和升级含义 | 按本文核对开发基线；固定上游 revision；准备精确 QC fixtures、规范化/身份预期和公开数据 manifest；不把依赖准备当作科学验收 |
| M1 peak caller 重复一致性 | 明确选择的分子不会被 MACS3 隐式再裁剪 | 两个 LB 各一个同坐标模板；retain 中含同坐标模板；普通/pooled/consensus/IDR-ready 的 MACS3 参数与代表性真实计数；混合 assay auto 解析、显式覆盖与幂等；参数冲突与迁移；IDR 方法未被顺带更改 |
| M2 PE 插入上限 | 支持文库中的合法模板被正确配对与消费 | 真实 Bowtie2 500/501/600 及各新上限边界；两种 R1 方向；真实 MACS3 BAMPE 读入；不只验证手工 99/147 |
| M3 assay 默认和参考过滤 | 两个 profile 的 CUT&Tag 政策与默认 blacklist 生效，参考身份可验证 | auto/显式 retain/exclude；narrow/broad；blacklist 默认 analysis、显式 off/copies 和旧 alias；off 无资源可运行、on 缺资源/不匹配明确失败、PE 无孤儿；细胞器名称正常/RefSeq/无匹配/存在但零 reads；验证发生在结果产出前 |
| M4 QC 交付一致性 | 统计阶段与分析产物对应，异常模板可见 | 负 TLEN、零 TLEN、重复、缺失 mate；MNase 边界和窗口并集；mito_fraction 可在平台列表/下载/QC 与 MultiQC 读取；严格解析器的版本兼容 |
| M5 cross-correlation | phantom 主导输入不产生无条件合格，强真实峰仍正确 | 同输入单变量窗口对照、纯 200 bp 对照、实际 R 多候选输出、变长输入与未知目标；BAM/原规则可达证据及非 ChIP 解释状态 |
| M6 集成与发行依据 | 用户正常建项目即使用所声明的修复，旧输出不混用 | CLI/authoring/materializer 一致；新旧 outdir 拒绝混用；受控 tiny 科学执行；真实数据比较；适用身份与契约检查；无遗留服务 |

每个单元对应一个主要消费者结果，按实际 diff 控制大小；必要时把纯机械契约准备与科学行为变化分开。不要为每个单元重复全套浏览器/平台/科学测试。

实施依赖按以下顺序推进，不按表中数字机械串行：

1. M0 后先落配置解析、内部政策修订和两个 profile 的 outdir guard；建立新旧契约 fixtures，防止后续科学规则产生混合结果。
2. M1 与 M3 的重复政策解析一起完成，再把共享 MACS3 参数接入全部路径。调用点逐项做运行时/调度覆盖清单，不能以 rg 未发现旧参数代替验收。
3. M2 落显式上限与旧参数迁移；M3 完成分析视图、模板闭合、参考校验及消费者接线。
4. M4 的 writer/reader/artifact/MultiQC 一起落地；M5 接入相同阶段和版本化解释状态。
5. M6 汇合 tiny 科学执行、入口兼容和真实数据单变量比较。先验证修复机制，再评估 assay 默认；按结果冻结 CUT&Tag 候选默认，记录例外与适用边界。

M0–M5 的定向实现与测试可以在真实数据准备期间推进；正式宣称默认科学方案已验收，必须完成对应 M6 证据。遇到真实数据不支持候选政策时，保留结果并回到该项决策，不靠降低测试要求通过。

真实数据验证至少选 ATAC PE、CUT&Tag PE，并选择有代表性的 ChIP PE/MNase 数据作为其默认值与指标的证据。数据的 accession/路径、参考身份、选取方式、固定子采样、资源预算和评价阈值须先记录。不得自动读取用户本地测序数据或进行大型下载。

验收不以“峰更多”作为改进证明。应分别比较配对率、有效模板数、片段分布、重复率、信号/背景、重复间一致性；小型子采样只证明受控边界，不能取代充分深度的生物学评估。

真实数据的比较按固定基线拆开变量：先仅改变 MACS3 keep-dup，再仅改变 max_insert，再比较 CUT&Tag retain/exclude，最后做组合回归。相同数据、参考、下采样和其余参数保持一致；不把多项变更混合后的峰数变化归因于其中一项。

所有默认值变更至少具备两层证据：受控真实工具数据证明实现与边界正确；来源固定的代表性数据评估实际结果是否支持该政策。插入上限的权威基线已选定，不等待每文库寻优后才实施；CUT&Tag 重复政策仍按已确认方案验证后冻结。真实数据评估未完成时，可报告基础修复与工具测试完成，不能宣布整体科学验收通过。

### 5.1 最低测试矩阵与证据交付

| 层级 | 必须证明的行为 | 现有测试边界 |
| --- | --- | --- |
| 配置/authoring | 两个 profile、混合 assay、auto/显式覆盖、alias 冲突、两次规范化、合法整数与 extra_args 迁移 | test/config/，test/services/test_input_authoring.py 及相关 adapter authoring 测试 |
| 脚本/严格解析器 | TLEN 正负/零、缺失 mate、重复保留、MNase 窗口并集、合法零与 NA、新旧精确 header | test/scripts/，test/test_parse_cross_correlation.py，test/adapters/test_encode_qc.py |
| 真实工具合成输入 | Bowtie2 距离边界及方向；MACS3 两个 LB 同坐标、retain、普通与派生路径；SPP phantom 与真实峰对照 | test/real_execution/ 中对应测试；新增真正由 Bowtie2 产生的模板，不能只手写 proper flags |
| DAG/小型执行 | 两个 profile 的 on/off/copies；各消费者实际读取预期 BAM/峰；Mito/QC 可提取；旧目录拒绝 | test/workflow/，已有 analysis-chain/blacklist/fragment tiny 执行 |
| 平台交付 | 新旧 QC 均可列出/读取、数值与下载一致、缺测状态准确、无私有路径泄露 | adapter/artifact/QC service 定向测试；仅涉及实际 UI 行为变化时扩展浏览器验证 |
| 真实文库 | 固定输入、参考和其余参数，分别比较 keep-dup、max_insert、CUT&Tag 重复政策，再组合回归 | 独立证据报告与冻结的数据 manifest；不把测试通过数等同科学有效性 |

按最终改动运行相关 validator、DAG dry-run、脚本/适配器测试及 real_execution；命令使用仓库 canonical bootstrap。科学相关 schema/快照需检查精确 diff；改 OpenAPI 才 regenerate client。通过后不无故重复全套测试；完整发行 Gate 由最终闭包影响决定。本轮文档阶段没有运行上述修复验收。

交付证据包含：代码/环境身份、输入 manifest、完整命令、定向测试报告、有效政策清单、各阶段计数、旧/新结果差异和未完成项。真实工具的边界断言必须精确；真实文库的统计对比按下节预先登记，禁止事后只挑有利指标。

### 5.2 公开数据与真实文库验收准备

默认使用公开数据准备验收，不读取用户私有样本。当前仅做元数据选择，尚未下载 FASTQ、建参考或运行全流程。候选与可用证据明确分开：

| assay | 当前选择入口 | 已知事实与尚需核验 |
| --- | --- | --- |
| ChIP PE | ENCODE 的 PE TF/histone 数据；可从 ENCSR375VXU 的对应分析记录开始核查 | 需要文件级原始 FASTQ 配对、重复/对照关系和参考身份；不凭单个分析对象推断完整实验设计 |
| ATAC PE | ENCODE ATAC 官方 ENCSR356KRQ 示例及其原始实验记录 | 官方示例可作为工具集成入口；须核验布局/文件、取样范围和真实长片段是否保留，染色体级/稀疏子集不能独立冻结全基因组默认 |
| CUT&Tag PE | GSE145187；已核实 GSM4308159（H3K4me3，SRX7713680）与 GSM4308163（H3K36me3，SRX7713684）的 PE 方法记录 | 两者为不同靶标，不能当成一对重复；在同系列中核验匹配重复，并加入不同富集形态和低输入/高重复的条件。若没有这些条件，明确未验证的适用范围 |
| MNase PE | 优先核查 nucMACC 论文数据清单中的公开 MNase 数据 | 核验真实 PE、消化条件和尺寸选择；mono-only 文库不能验证 di 分支；不以 MNase HiChIP/连接型数据替代普通 MNase-seq |

仓库旧计划中的 accession 仅作为线索。例如 ENCODE ChIP 官方 ENCSR000DYI 示例明确为 SE，不能把它当成 PE 插入上限验收数据。网页/仓库示例不能直接证明每个原始文件的布局与重复关系，实施前必须冻结以下 manifest 字段：experiment/sample/run/file accession、assay、target、layout、mate 对应关系、biological/technical replicate、control、原始 URL、校验和、read 长度、参考 assembly/资源身份、尺寸选择/输入量等已知条件、固定抽样方法/seed、输入 read-pair 数及计划资源预算。

先做小规模配对保真的工具验证，再运行有足够深度的结果比较；抽样以模板为单位且两侧同步，不从不同 FASTQ 独立抽样。下采样会改变重复碰撞和找峰统计，不允许仅凭低深度子集冻结 CUT&Tag 去重默认。至少 ATAC、CUT&Tag 要有可比较的真实重复；ChIP 覆盖 PE TF 与 histone 两类，MNase 覆盖实际支持的 mono/di 窗口。可分批验证，未覆盖的 assay/条件单列，不宣传已验证。

对比在运行前登记：

- 固定输入、参考、软件、其他参数和资源；先看单变量，再看组合。
- 报告配对率、有效模板数、MAPQ/片段长度分布、原高置信短模板的配对坐标变化、峰与 coverage 差异、信号/背景及可用重复间一致性。ATAC 关注长片段尾部和 nucleosomal 结构；CUT&Tag 分别评估 retain/exclude 的信号和背景，不能以重复率较高本身证明应排除。
- 记录 wall time/内存，解释宽配对范围的代价。工具报错、非法模板泄漏、配置不生效、计数不守恒、交付/QC 阶段不一致属于硬失败。
- 不预设缺乏依据的通用“峰更多/相关性 ≥ 某数即通过”。对真实生物学指标，先冻结对照/计算定义，再由差异报告确认是否出现无法解释的背景增加、重复一致性下降或既有高置信信号损失；存在此类问题则不冻结对应默认，追加归因比较并记录结论。该审核结论与数值必须同时交付。
- 下载前在 manifest 中记录输入规模、现有参考复用方式及磁盘/CPU/时长预算。当前未授权启动大型下载；若完整科学验收超出可用资源，保留未验收状态并提供具体数据/预算安排，不把 tiny 运行升级称为全量验证。

公开数据的文件级冻结和运行预算是执行准备任务，不是新增的科学政策选择；来源不可访问时如实记录，选择同等实验条件的可验证数据，不伪造已验证 accession。

## 6. 决策状态与执行关卡

| 决定 | 当前推荐 | 状态 |
| --- | --- | --- |
| 首轮范围 | 默认路径和 QC；IDR 方法对齐后续 | **用户已确认** |
| legacy 已确认错误 | 同样修复，历史结果原样保留 | **用户已确认** |
| CUT&Tag auto 政策 | retain 为候选，验证后冻结 | **用户已确认纳入本轮验证**，最终默认待数据 |
| 两个 profile 的 assay 默认矩阵 | ChIP/ATAC/MNase exclude，CUT&Tag retain；legacy broad ChIP 同样调整；允许显式覆盖 | **用户已确认按 assay 统一候选默认**，真实数据验证后冻结 |
| CLI 显式覆盖 | ATAC 可显式 retain，覆盖 assay 默认且贯穿 BAM/coverage/peak caller；明确全局配置作用域和冲突处理 | **用户已明确要求此能力**，纳入验收 |
| IDR 的公共参数接线 | 共同上游和 MACS3 重复政策生效；不改 oracle/选择/判级；覆盖参数回归，承认输入与结果可能变化 | **纳入首轮全 MACS3 一致性范围**；IDR 方法对齐仍后续 |
| blacklist 默认与缺资源行为 | 可选且默认开启分析主链；显式关闭不要求资源；开启时资源缺失/身份不匹配不能静默跳过；明确旧 copies/alias 迁移含义 | **用户已确认可选、默认开启**；上述资源与迁移规则为本 spec 的建议实现 |
| 主 fragment-size QC 阶段 | 跟随实际分析输入：blacklist off/copies 为未做 blacklist 过滤的分析 BAM，analysis 为过滤后的分析 BAM；复杂度/细胞器比例保留各自诊断阶段 | **用户已确认按 blacklist 是否作用于主链选择实际分析输入** |
| assay 插入上限 | ChIP/ATAC/CUT&Tag/MNase 分别为 2000/2000/700/500；尊重合法显式覆盖 | **用户已确认先对齐权威流程安排本轮默认基线**；实现与影响待验证 |
| 细胞器 assay 默认 | 两个 profile 的 ATAC/CUT&Tag exclude，ChIP/MNase report，允许覆盖 | **用户已确认** |
| cross-correlation 未知目标/非 ChIP | 展示诊断，不套 ChIP 阈值判合格；维持 opt-in | **用户已确认**；输入和窗口按 3.5 具体设计 |
| 无匹配细胞器名称 | exclude 失败；report 为 unresolved + NA；确实不执行过滤可显式 keep，不把空列表当成生物学证明 | 技术方案已写入 3.4，随整体 spec 审阅 |
| 真实验收数据 | 公开数据优先，候选见 5.2；冻结文件级 manifest、配对关系、深度和预算 | **执行准备关卡**，尚未完成或运行；不阻塞定向实现，但阻塞整体科学验收 |
| 新旧运行身份/输出契约版本 | 内部政策修订、两个 profile marker 防混用；authoring v3、QC v2、保留旧 QC 读取 | 技术方案已写入 3.6/3.7/第 4 节，待实现与兼容测试 |

主要产品决策已经齐备，用户已安排进入分支开发，不再重复要求逐项或整体授权。开发者先核对本文、相关代码及测试，再按交接说明逐批实现；若发现实质性冲突，在 PR 中明确指出。CUT&Tag 候选默认的数据冻结、具体公开文件 manifest 和最终发行 Gate 仍是有证据要求的执行关卡，不能因为本文已完整就标为完成。

## 7. IDR 后续事项

首轮之外保留：Nt/Np 选择与 conservative/optimal 名称、真实/pooled/self oracle、matched/pooled control 比较政策、IDR 专用 pre/post-blacklist 计数和判级语义。它们需要单独 spec 与真实 IDR 拟合验收，不能通过首轮通用修复宣布已经对齐 ENCODE。

## 8. 证据与一手依据

- 两轮独立报告及本会话真实工具复核证明了一组可复现机制；合成输入不能证明真实文库中的效应量。
- [ENCODE ChIP v2.2.2 Bowtie2 wrapper](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/v2.2.2/src/encode_task_bowtie2.py)：PE -X2000。
- [ENCODE ATAC Bowtie2 wrapper](https://github.com/ENCODE-DCC/atac-seq-pipeline/blob/master/src/encode_task_bowtie2.py)：PE -X2000；实现前把当前已核验内容固定到具体 revision，避免以浮动 master 作为长期身份依据。
- [ENCODE ChIP v2.2.2 MACS wrapper](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/v2.2.2/src/encode_task_macs2_chip.py)：显式 --keep-dup all。
- [ENCODE ChIP v2.2.2 输入契约](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/v2.2.2/docs/input.md)：cross-correlation 输入及 TF/histone 排除范围。
- [CUT&Tag 方法团队教程](https://yezhengstat.github.io/CUTTag_tutorial/)：-X700 示例、apparent duplicates 的解释与条件去重。该教程是候选依据，不等于所有文库的唯一政策。
- [SAM/BAM 规范](https://samtools.github.io/hts-specs/SAMv1.pdf)：TLEN 的符号与不可用零值。
- [Bowtie2 官方手册](https://bowtie-bio.sourceforge.net/bowtie2/manual.shtml)：concordant 配对的方向/距离约束、MAPQ 的解释、扩大 -I/-X 搜索窗口的运行代价。
- [deepTools bamCoverage](https://deeptools.readthedocs.io/en/3.5.6/content/tools/bamCoverage.html)：MNase 中心信号的片段长度窗口；当前仓库显式渲染 dyad=130–200，与 mono BAM=140–200 分开。
- [RCSB PDB NCP147 结构](https://www.rcsb.org/structure/1KX5)：147 bp 核小体核心 DNA 的结构依据。
- [nucMACC MNase-seq 发表流程](https://pmc.ncbi.nlm.nih.gov/articles/PMC11221511/)：Bowtie2 -X500；其他 MNase 论文也有 -X1000，不将 500 描述为唯一行业标准。
- [GSM4308159 的 CUT&Tag 原始方法记录](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4308159)、[GSM4308163 方法记录](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSM4308163)：H3K4me3/H3K36me3 的 PE 示例，不能据此宣称已完成数据下载或默认政策验证。
- [ENCODE ChIP 的 ENCSR000DYI 示例](https://github.com/ENCODE-DCC/chip-seq-pipeline2/blob/v2.2.2/example_input_json/ENCSR000DYI_subsampled_chr19_only.json)：paired_end=false，用于识别旧数据计划中的布局误判。
- [Snakemake 配置与命令行覆盖](https://snakemake.readthedocs.io/en/stable/snakefiles/configuration.html)：--config 的显式配置覆盖及配置合并；具体嵌套输入形式由本仓库锁定版本的入口测试验证。

本轮验证仅为文档检查；没有新增或重跑科学实验，没有执行实现、重建 manifest、打包或发布。

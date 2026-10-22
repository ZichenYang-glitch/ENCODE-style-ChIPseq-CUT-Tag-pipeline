# 科学计算修复：Kimi 开发 / Codex 验收交接

实施依据：[scientific-computation-remediation-spec.md](scientific-computation-remediation-spec.md)。
科学审查基线：fa0d8401cbc3fecb765fed6d7e01e18bb8fbc96d。
用户已安排 Kimi 负责实现、Codex 负责独立验收，沿用独立 branch + PR 的交付方式。本文是执行交接，不代表任何修复或科学验收已经完成。

## 工作区与分支

- 阅读 AGENTS.md、spec 和相关实现/测试后开发。根工作区 main 已有用户暂存和未跟踪文件；不得 reset、stash、提交、覆盖或清理这些既有内容。
- 从包含本 spec 的已确认基线创建独立 worktree/分支。文档 PR 未合并时，可从其 head 新建实现分支，并把实现 PR 的 base 指向文档分支，清楚声明依赖；不要把整份 spec 再复制提交进每个实现 PR。
- 每个 PR 对应一个完整消费者结果。禁止直接在 main 实现或一口气混入全部修复。依赖 PR 尚未合并时，记录准确 base/head SHA；不擅自 rebase、强推、合并或删除分支。
- 用户已授权本任务创建开发分支、提交本任务代码并开 PR；不包含自动合并、部署或发布。提交保留人类作者身份，按 AGENTS.md 记录实际参与模型的默认 co-author trailer，不虚构参与者。
- Codex 当前负责文档交接和后续验收，不代替 Kimi 提交科学实现。Kimi 未实际参与的当前文档提交不填写其虚构 co-author；用户提供的 Kimi/Peri 报告作为审查线索与来源说明，执行证据仍须复核。

## 分批 PR

以下 SCI 编号只用于本次修复，不替代既有 roadmap 编号。分支前缀建议 `agent/science-`。

| 批次 | 建议分支后缀 | 一个主要交付结果 | 范围与关键验收 |
| --- | --- | --- | --- |
| SCI-0 | remediation-plan | 维护者可复核的实施依据 | spec + 本交接；仅文档，由 Codex 准备，不含根工作区既有改动 |
| SCI-1 | execution-identity | CLI 运行参数可追溯，修正版不混写旧结果 | spec 第 4 节和 3.7 的有效配置记录；两个 profile 的 marker、同身份继续、旧/损坏/symlink marker 拒绝、CLI 合并配置校验与摘要、dry-run 不写运行目录。暂不激活后续尚未接通的科学默认，不创建无消费者的配置框架 |
| SCI-2 | duplicate-policy | 用户的重复保留/排除选择贯穿 BAM 与 peak caller | 混合 assay auto、两个 profile、legacy broad、显式覆盖、全部 MACS3 keep-dup all（含 IDR-ready）、参数冲突迁移；真实 MACS3 同坐标双 LB/retain 测试；不改 IDR 方法 |
| SCI-3 | paired-insert | PE 配对上限显式且可覆盖 | 2000/2000/700/500 默认、旧 -X 迁移、SE 不渲染；真实 Bowtie2 的边界与双方向、MACS3 BAMPE 读入和配置入口回归 |
| SCI-4 | reference-analysis | 参考过滤在两个 profile 的分析主链正确生效 | canonical analysis BAM、完整模板、assay organelle 默认、命名与 blacklist 资源校验、default analysis/显式 off/copies/alias；所有消费者接线及真实模板过滤测试 |
| SCI-5 | qc-contract | 用户看到与实际分析阶段一致的 QC | fragment-size/MNase/mito、合法负 TLEN、无效模板计数、未分类并集；QC v2 writer/reader、artifact、平台提取与 MultiQC 一起更新，旧 v1 可读 |
| SCI-6 | cross-correlation | phantom 主导输入不被无条件判合格 | 输入/R1、实际长度估计、TF/histone 窗口、unknown/非 ChIP 诊断、多候选 parser；真实 R + 生产 parser 及规则可达验证 |
| SCI-7 | scientific-validation | 有可追溯证据的默认政策与集成交付 | 文件级公开数据 manifest、固定变量比较、真实数据与 tiny 集成、CLI/平台一致性、身份/Gate 影响评估、默认冻结及升级说明 |

SCI-2 的重复政策和分析视图、SCI-4 的过滤接线可能共享现有规则。按消费者需要协调边界：不能先合入 retain 开关却仍让 legacy 上游永久删除重复，也不能先声明默认 blacklist 已生效却让消费者仍读旧 BAM。若无法保持一个可验证的完整结果，可调整相邻批次拆分并在 PR 说明，但不得扩大到 IDR 方法或通用平台重构。

每批更新必要 schema/文档/测试，使该批行为自洽；authoring schema v3 随首次实际公共契约变化落地，不能单独提前更名后宣称新政策可用。SCI-5 前已有 QC 不得被前面 PR 生成无法解析的新 header。

每批实现后提交 PR 供 Codex 验收，再推进依赖批次；未完成真实数据要求的候选默认不得标成已科学验收或直接作为发行完成。可在依赖分支推进实现与数据准备，清楚区分“代码与工具测试完成”和“默认冻结/发行关卡完成”。

## 第一轮给 Kimi 的任务

**现在只执行 SCI-1。**

1. 使用独立 worktree/branch，记录基线 SHA，确认只有本任务 diff。保持根工作区原样。
2. 阅读 workflow/Snakefile 的 identity 检查、analysis.smk 的 marker 规则、validator、authoring/materializer 和现有 execution identity 测试，确定最小改动；运行前先检查可用的锁定环境。
3. 给两个 profile 接入版本化 marker 和内部政策修订；身份覆盖本批实际有效配置。相同身份可继续，旧/不匹配目录明确拒绝，不提供覆盖旧结果的快捷绕过。
4. 沿用 Snakemake `--config`，验证其合并值进入 validator、运行身份和有效配置记录。输出安全的用户可读摘要；实际运行保存有效配置。不要新增未经设计的 `--set`/点号覆盖语法，不把 validate 文件模式说成已经校验了后续注入值。
5. 不在本批提前激活 auto 重复矩阵、MACS3 参数、max_insert、blacklist 或 cross-correlation 新策略；摘要如实显示本批实际行为。内部身份机制必须被真实规则消费，不做只有测试使用的抽象。
6. 增加有意义的正反测试：legacy/updated；相同/不同身份；旧格式/损坏/symlink；未标记旧输出；CLI 与 YAML 合并值；dry-run 不写 marker 或有效配置文件。用受控 fixture，勿读取用户测序数据或把私有路径写进公开日志。
7. 运行相关定向测试、必要 validator/DAG dry-run 和 git diff --check。缺少真实运行依赖时明确记录，不能加 skip 或伪造通过。公共 OpenAPI 真有变化才重生成 client。
8. 提交本批变更、推送本任务分支并开 PR。未完成规定验收的 PR 保持 Draft；不要合并。然后把下述交付信息给用户，交 Codex 独立验收。**不要顺带开始 SCI-2。**

公开数据只先做元数据/manifest 准备；不启动大型下载，不读取未明确授权的私有 FASTQ。后续使用真实工具不等于可以在共享环境安装临时依赖；遵守仓库 task-local locked environment 规则。

## 每个 PR 的交付信息

- PR URL、branch、base SHA、head SHA、相对 base 的改动文件和范围。
- 具体触发场景及修复前后行为，涉及的 spec 条目。
- 原样可复现的验证命令、工具/环境身份、实际通过/失败/未运行项目和证据位置。
- 科学行为变更、CLI/YAML/平台差异、兼容与升级影响；没有验证的条件明确列出。
- 本批结论分别写：实现完成情况、契约/工具测试情况、真实文库科学验收情况。后者未执行不能写成整体通过。
- PR description 围绕最终变更重写，不能把整份审计报告粘贴进正文或用测试数量代替证据。

## Codex 独立验收

以交付的准确 head SHA 为对象；先复核 base/diff 和 spec 对照，再在干净 checkout/worktree 复现关键证据。不得仅依据 Kimi 自报通过或既有测试锁定当前行为就批准。

关注两种 profile 的默认/覆盖、消费者实际输入、执行身份、产物与 QC 一致性；科学调用变化使用真实工具合成边界复核，真实文库的效应另按 spec 验收。发现问题提供路径/位置、触发条件、影响和修正建议；结论区分阻塞缺陷、建议及缺失证据。

验收通过不自动合并、不免除后续默认冻结和发行 Gate。验收期间新增修复提交后，对新的准确 head 复查相应变更；不把旧 head 的结论沿用到新提交。

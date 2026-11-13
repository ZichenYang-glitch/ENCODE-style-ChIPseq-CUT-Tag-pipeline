# HelixWeave Product Roadmap

## 当前入口：2026-11-13

本文件记录当前产品边界、优先级和交付顺序。准确 PR/SHA/CI 快照见
[集成清单](current-delivery-integration.md)，实施分工见
[session 计划](next-delivery-session-plan.md)，新设计见
[数据目录、Agent 与 Viewer 边界](../architecture/catalog-agent-viewer-boundaries.md)。
本轮交付文档、设计及定点复审，没有修改产品、运行迁移、新的科学数据批次、Gate 或发布。
复审仅运行隔离 API/文件访问和 tiny Picard 回归，范围见集成清单。

远端 main 查询于 2026-11-13 09:31 UTC，为
`5bc21a040951ba28333bcc04dfb151735089fcd6`（#215 merge）。
2026-11-13T10:26:21Z 再次查询，main 与 #212、#214、#216、#217、#213 的 head
仍是集成清单中的那些 SHA，因此那次审阅只适用于这些 SHA。
根工作区仍为 `ff303296…`，其既有未提交修改未动。
当前入口取代历史文档中“当前”“下一步”的指向；历史实验、失败、SHA 和裁决不追改。

## 产品边界与已确认决策

- 保留一个仓库、SQLite 和现有文件存储。SQLite 管身份、元数据和生命周期；
  BAM/bigWig/PET 等字节在文件存储中，Redis/RQ 只承担执行交接。
- 数据目录、科学 adapter、计算执行、Agent 接入和 Viewer 有明确接口及写入归属。
  查询和已有结果浏览不以计算 worker 或 Peri 在线为前提。
- 保留实验室成员共享项目模型，不在本轮引入项目 ACL、多租户或新数据库。
- AI 整理：成员可在选定项目内查询和生成建议；管理员可在网页或终端审阅、
  批量确认。两种入口共用 API、权限检查、版本比较和回执。Peri 不取得数据库管理权限。
  此规则只针对新增批量整理，不顺带改变已有手工编辑或 A1 草稿采用权限。
- Viewer 先用 WashU；在本仓库独立构建并嵌入平台是下一阶段目标，当前 #217
  仍与平台前端共享构建。无需另建 HelixViewer 仓库；IGV 等第二引擎等具体需求。
- 科学分析已决定收成单一路径，不保留 legacy 第二政策或旧键别名过渡期。
  auto 候选矩阵为 ChIP/ATAC/MNase exclude、CUT&Tag retain；Picard 只标记，
  分析 BAM 执行政策，MACS3 显式 `--keep-dup all`。#216 最新 Picard 返修定点复审通过，
  整体合入审阅仍待完成，不等于默认已获生物学冻结。
  完整范围和不顺带改变的参数见[原决定](workflow-platform-agent-roadmap-history-through-20261112.md#科学分析单路径收口2026-11-12已决定待-cursor-实施)。

## 已交付与尚未验收

已进入 main 的基础包括项目/样本及不可变修订、受管输入和参考身份、运行/取消、
结果发布与下载；科学实现及 SCI-7 本批记录、skills、源码身份和两批测试精简；
#210 取消竞态修复、#211 固定边界内 A0 的归档、#215 浏览器依赖镜像修复。

SCI-7 本批记录已交付，结论仍为 no conclusion。JUND 三臂及登记比较已完成，
是跨靶点检查，不是 NR2F2 的生物学重复。组蛋白、ATAC、可选 CUT&Tag/MNase
confirmation 未执行；MNase 原登记 500-arm 指标未评估。本计划不启动新数据批次。

#212 的最新边界复审通过，可供 C3 消费；C4 和整个 A1 未完成。
#214 的 `2c3c439…` 正在把目录接到真实 snapshot，该 head 的 CI 为 cancelled，
不能记为成功。#216 单路径收口、#217 WashU 仍在 Draft。
普通 CI 成功、分支实现、独立审阅和产品闭环是不同状态。
A0 通过不代表 A1、真实 provider、生产容量或从零构建通过。

## 后续交付计划

每项记录范围、依赖、下一步、验收标准和关联 PR。细节以设计及 session 计划为准。

| 工作线／范围 | 依赖 | 下一步 | 验收标准 | 关联 PR |
| --- | --- | --- | --- | --- |
| A1：受限配置建议与校验 | #212 `e00130e…` 的最新边界复审通过，可供 C3 消费 | #214 接目录到草稿/快照，再完成 C4；不加入新分类系统 | 真实 HTTP 保存、回读、snapshot、回执及固定 native/stub；拒绝路径有效；最终组合 CI 和适用 Gate | #212、#214；集成 PR 待创建 |
| 科学单路径 | 已确认政策；#216 `c44c4d6…` Picard 返修定点复审通过 | 整体合入审阅；集成时复核 A1 schema | 一套政策、旧键明确拒绝、历史身份不改；实际 Picard/分析 BAM/MACS3 契约一致 | #216 |
| WashU 首版 | #217 最新 head 独立复审 | 收口授权 Range、产物/索引身份、有限转换和真实绘图 | 已知窗口数据可见；错误、空窗口、截断和缺失区分；保留原科学产物 | #217 |
| 数据目录 M1 | 本轮设计；A1 共享 API/迁移接点稳定 | 字段契约复核后开发课题/样本/实验分类与查询 | 手工管理先可用；组合筛选/分页准确；旧修订不改写；不夸大结果归属 | 独立 PR 待创建 |
| Viewer 隔离 V2 | WashU 首版通过；版本化轨道契约 | 同仓独立构建与样式/依赖隔离，保持同源认证 | 独立构建可嵌入；不直接查 SQL/调执行；worker/Peri 停止后仍可读已有结果 | 独立 PR 待创建 |
| Peri 收口 P1／真实模型 P2 | A1 C4；固定 runtime/权限边界 | 收紧科学适配和 repository 依赖；真实 provider 独立接入 | stub/真实 provider 可区分；外发/超时/撤权边界有效；固定新身份验收 | #212 后续独立 PR |
| AI 元数据整理 M2 | M1、P2、目录投影及工具资格 | 先读与建议，再接管理员批量采用；网页/终端共用契约 | 无 SQL 工具；有来源、差异、版本、期限和回执；撤权/冲突拒绝写入 | 独立 PR 待创建 |
| Agent 运行操作 | A1、真实 provider、业务确认设计 | 设计 create/preflight/start/cancel 与既有生命周期的映射 | 不另造任务 API；确认范围、重试、取消请求/确认及审计有明确契约 | 尚未开始 |
| 文档、skills、测试与 CI 运维 | 实际反馈及已确认缺陷 | 统一入口；核对单路径后的 skills；按具体收益处理维护 | 历史不倒改；skip 不当通过；不为省测试弱化契约 | 本轮文档；#213 待协调；#194/#195 待审 |

## 顺序与并行规则

1. 现有五个 Cursor session 为 Backend、Frontend、Science、Viewer 和维护。
   前四条继续当前闭环，维护接手本轮文档并协调 #213；加 Codex 主会话共六个。
   Codex 复审准确 head；CI 成功不会自动冻结接口或批准合并。
2. #216 可独立验收；A1 最终组合必须使用拟交付的科学 schema/adapter。
   #212 与 #214 有祖先依赖，集成前检查 DAG，避免重复补丁。
3. M1 现在可准备设计/fixture；共享后端契约稳定后再开始数据库与 API 写入。
   M1、V2、P1/P2 在写入目录互不重叠且各自依赖满足后可并行。
4. M2 等 M1/P2；运行操作另设目标。不给当前 A1 无限追加能力。
5. 每条 session 使用独立 worktree；迁移/OpenAPI/生成客户端有一个明确写入者。
   Codex 负责设计与独立审阅，Cursor 负责 coding、测试及交付。
   subagent 可做独立检查或不相交的实现；主 session 负责集成和最终核验。

## 后置工作与交付纪律

不把 PostgreSQL、对象存储、多机执行、长期 Peri fork、自由文本、第二 viewer、
新的 confirmation、科学默认冻结或发布作为本批隐含工作。
它们保留为后续决策，不因本设计自动授权计算、下载、合并或发布。
早期维护条目先按当前代码复核，不把历史待办直接重开为新目标。

测试从实际变更层开始，CI 绑定准确 head，并记录 PR 合成 checkout。
Gate 按最终执行闭包和身份影响裁决；相同闭包的复用要给出证据，不能将历史
Gate 写成新 SHA 的执行。禁止降门槛、加 skip 或无依据反复重跑。

#213 保存的是较早状态快照。本轮文档与它需由文档 owner 协调一次纳入，
不并行维护两份“当前”，本轮不关闭或合并它。根工作区已有的单路径决策已逐字保全。

## 历史入口

[整理前的完整 roadmap](workflow-platform-agent-roadmap-history-through-20261112.md)
与[旧集成记录](current-delivery-integration-history-through-20261112.md)保留原字节。
其中“当前”“未开始”“待授权”只描述原记录时点。本轮未重写 A0、科学或 Gate 证据。

以下锚点保留给已有文档链接；它们是历史入口，不是当前排期。

<a id="已确认的维护顺序与审核边界2026-09-22"></a>

- [2026-09-22 维护顺序原记录](workflow-platform-agent-roadmap-history-through-20261112.md#已确认的维护顺序与审核边界2026-09-22)

<a id="后续-hi-trac-预处理-adapter五个-pr"></a>

- [Hi-TrAC 五阶段原计划](workflow-platform-agent-roadmap-history-through-20261112.md#后续-hi-trac-预处理-adapter五个-pr)

# HelixWeave Product Roadmap

## 当前入口：2026-11-14T17:40:08Z

本文件记录当前产品边界、优先级和交付顺序。准确 PR/SHA/CI 快照见
[集成清单](current-delivery-integration.md)，实施分工见
[session 计划](next-delivery-session-plan.md)，新设计见
[数据目录、Agent 与 Viewer 边界](../architecture/catalog-agent-viewer-boundaries.md)。
本次只更新当前状态，没有修改产品、运行迁移、新的科学数据批次、Gate 或发布。

2026-11-14T15:57:43Z 的检查细节保留在集成清单：当时 workflow_dispatch CI 仍 in progress，
随后于 16:12:17Z failure。2026-11-14T17:40:08Z 再次查询，main 仍为
`101983f9480bb72306da69875e79dc6ad65bca5a`，父提交为 #215 merge `5bc21a0…`
与 #216 head `85ac7a2…`。#216 的 Pulls API 仍是 open、merged=false，禁止再次合并。
该 merge SHA 仍没有 push 检查，也不把 `85ac7a2…` 的 pull_request 检查记成它的检查。
同一窗口 #212 head 未变；#214 head 未变；#217 head 已从审阅对象前进。
head 没有移动时，原审阅继续只适用于该 SHA；head 已经前进的，不继承旧裁决。
根工作区仍为 `ff303296…`，其既有未提交修改未动。
当前入口取代历史文档中“当前”“下一步”的指向；历史实验、失败、SHA 和裁决不追改。
本文件已在 Draft #213 中，该 PR 仍未合并。本次更新不把随后生成的提交 SHA 写回本文件。

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
  分析 BAM 执行政策，MACS3 显式 `--keep-dup all`。`analysis_profile` 与 `remove_dup`
  会被拒绝，不能静默改义。narrow/broad 不改变 auto 选择；显式 exclude/retain 作用于整次运行。
  外部 `control_bam` 不被接受。科学政策修订号是 7，新运行使用新的 outdir，历史结果目录不重写。
  #216 `85ac7a2…` 已通过整体独立审阅，可作为 A1 的已审阅科学 schema 目标。
  main 上已有包含该 head 的 merge 提交，但查询时 PR 对象仍 open、未标记 merged。
  这不等于 C4、A1 完成，也不等于科学默认已获生物学确认。
  完整范围和不顺带改变的参数见[原决定](workflow-platform-agent-roadmap-history-through-20261112.md#科学分析单路径收口2026-11-12已决定待-cursor-实施)。

## 已交付与尚未验收

已进入 main 的基础包括项目/样本及不可变修订、受管输入和参考身份、运行/取消、
结果发布与下载；科学实现及 SCI-7 本批记录、skills、源码身份和两批测试精简；
#210 取消竞态修复、#211 固定边界内 A0 的归档、#215 浏览器依赖镜像修复。

SCI-7 本批记录已交付，结论仍为 no conclusion。JUND 三臂及登记比较已完成，
是跨靶点检查，不是 NR2F2 的生物学重复。组蛋白、ATAC、可选 CUT&Tag/MNase
confirmation 未执行；MNase 原登记 500-arm 指标未评估。本计划不启动新数据批次。

#212 head 仍是 `e00130e…`。2026-11-13 的边界复审继续只适用于该 SHA，可供 C3 消费；
本轮没有新的源码审阅。C4 和整个 A1 未完成。
#214 head 仍是 `399c3b3…`。独立复审结论是需要返修，不是 C4 或整个 A1 验收：
PE 改为 SE 并确认保存后草稿仍 dirty；已有计划时，样本复选框和实际保存成员不一致。
原定向测试 100 passed；独立新增的 3 项回归断言失败。`2c3c439…` 的 cancelled CI 仍只属于旧 head。
#216 代码已在 main。PR 对象仍 open，禁止再次合并。
`85ac7a2…` 仍是已审阅科学 schema，但 dispatch CI 的科学真实执行失败仍是剩余阻塞。
源码核对见测试 PATH 遗漏 Picard，并仍保留旧 profile 分流和复制 BAM 的测试准备。
不能把 `--maxfail=1` 下已经跑过的数量写成整个套件通过。
下一步是独立小 PR 修复该测试层，并验证完整真实执行。
#217 对 `cbe5a0e…` 的独立审阅已完成：多成员 gzip 在首个 member 后漏掉后续数据，
PR 范围空白检查不干净。分支随后到 `414872e…`，17:40:08Z 的 head 是 `aa0ef9f…`。
这两个新提交都还没有正式交付后的复审，不继承 `cbe5a0e…` 的裁决，也不提前算修复通过。
普通 CI 成功、分支实现、独立审阅和产品闭环是不同状态。

原已接受的 Peri 运行 binary `60144bbe…` 现已缺失。新构建候选
`5cc82caf…` 的 2026-11-13 只读复核允许它进入隔离的候选资格验证，
尚未接受为 A0、A1 或 C4 runtime。它不能继承原 binary 的 native 验收，
也没有恢复原身份。宿主支持 systemd 和 cgroup v2，但用户服务尚未启动。
Backend 已有限定指令，可以启动 `user@1000.service` 并验证任务资源隔离；
本文件不记录这项操作已经执行。只有拿到实际执行证据后，才能写成环境合格或候选通过。
历史 A0 在固定版本上的验收结论保留，不改写成失败。
A0 通过不代表 A1、真实 provider、生产容量或从零构建通过。

## 后续交付计划

每项记录范围、依赖、下一步、验收标准和关联 PR。细节以设计及 session 计划为准。

| 工作线／范围 | 依赖 | 下一步 | 验收标准 | 关联 PR |
| --- | --- | --- | --- | --- |
| A1：受限配置建议与校验 | #212 `e00130e…` 的边界复审仍只适用于该 SHA；科学 schema 目标为已审阅的 #216 `85ac7a2…` | #214 先返修 PE→SE 保存基线和样本成员一致性，再复审；runtime 候选另行完成资格验证。不加入新分类系统 | 确认后的规范化草稿变为干净基线；显示选择与保存成员一致；真实 HTTP、snapshot、回执及完成资格验证的 native/stub；拒绝路径有效 | #212、#214、#216；集成 PR 待创建 |
| 科学单路径 | #216 代码已在 main，PR 仍 open；真实执行测试仍有剩余阻塞 | 独立小 PR 补上测试 PATH 中的 Picard，去掉旧 profile 分流和复制 BAM 的测试准备，再跑完整真实执行层。禁止重复合并 #216 | 完整真实执行层通过。不以 `--maxfail=1` 的已运行数量代替全套。不是生物学确认 | 修复 PR 待创建 |
| WashU 首版 | `cbe5a0e…` 审阅要求返修；当前 head `aa0ef9f…` 未复审 | 正式交付新 head 后复审 gzip 全 member 和范围空白。不把中间提交写成已通过 | 多成员 gzip 不丢后续记录；范围空白检查干净；已知窗口可见且错误、空窗口、截断和缺失可区分 | #217 |
| 数据目录 M1 | 本轮设计；A1 共享 API/迁移接点稳定 | 字段契约复核后开发课题/样本/实验分类与查询 | 手工管理先可用；组合筛选/分页准确；旧修订不改写；不夸大结果归属 | 独立 PR 待创建 |
| Viewer 隔离 V2 | WashU 首版通过；版本化轨道契约 | 同仓独立构建与样式/依赖隔离，保持同源认证 | 独立构建可嵌入；不直接查 SQL/调执行；worker/Peri 停止后仍可读已有结果 | 独立 PR 待创建 |
| Peri 收口 P1／真实模型 P2 | A1 C4；原 binary 缺失，候选未通过资格。用户服务尚未启动 | 先取得启动 `user@1000.service` 和任务资源隔离的实际证据；再按设计说明收口专用适配层。P1 尚未实现 | 环境合格和候选通过都要有执行证据；更换 runtime 不改业务 API，但新实现单独验收。不新增通用 Agent 框架 | #212 后续独立 PR |
| AI 元数据整理 M2 | M1、P2、目录投影及工具资格 | 先读与建议，再接管理员批量采用；网页/终端共用契约 | 无 SQL 工具；有来源、差异、版本、期限和回执；撤权/冲突拒绝写入 | 独立 PR 待创建 |
| Agent 运行操作 | A1、真实 provider、业务确认设计 | 设计 create/preflight/start/cancel 与既有生命周期的映射 | 不另造任务 API；确认范围、重试、取消请求/确认及审计有明确契约 | 尚未开始 |
| 文档、skills、测试与 CI 运维 | 实际反馈及已确认缺陷 | 当前入口保持一份；按具体收益处理维护 | 历史不倒改；skip 不当通过；不为省测试弱化契约 | 本状态在 Draft #213，未合并；#194/#195 待审 |

## 顺序与并行规则

1. 现有五个 Cursor session 为 Backend、Frontend、Science、Viewer 和维护。
   维护继续在 Draft #213 更新当前入口，不另开文档 PR。加 Codex 主会话共六个。
   Codex 复审准确 head；CI 成功不会自动冻结接口或批准合并。
2. #216 的代码已经在 main，禁止再次合并。A1 仍使用 `85ac7a2…` 这个已审阅 schema，
   但要等真实执行测试的小 PR 完成，并单独完成 C4。#212 与 #214 有祖先依赖，集成前检查 DAG。
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

当前入口已纳入 Draft #213，该 PR 仍未合并，也不并行维护第二份“当前”。
本次状态更新不关闭、不合并它，也不把本提交 SHA 写回文件。
根工作区已有的单路径决策已逐字保全。

## 历史入口

[整理前的完整 roadmap](workflow-platform-agent-roadmap-history-through-20261112.md)
与[旧集成记录](current-delivery-integration-history-through-20261112.md)保留原字节。
其中“当前”“未开始”“待授权”只描述原记录时点。本轮未重写 A0、科学或 Gate 证据。

以下锚点保留给已有文档链接；它们是历史入口，不是当前排期。

<a id="已确认的维护顺序与审核边界2026-09-22"></a>

- [2026-09-22 维护顺序原记录](workflow-platform-agent-roadmap-history-through-20261112.md#已确认的维护顺序与审核边界2026-09-22)

<a id="后续-hi-trac-预处理-adapter五个-pr"></a>

- [Hi-TrAC 五阶段原计划](workflow-platform-agent-roadmap-history-through-20261112.md#后续-hi-trac-预处理-adapter五个-pr)

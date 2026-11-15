# 下一批交付：session 分工与验收交接

日期：2026-11-13 的分工仍然有效。2026-11-15T07:21:49Z 只追加 Backend 与 Viewer 的最终交接。
2026-11-14T19:18:06Z 的分工留在集成清单快照，不在这里倒写。
排期入口为 [roadmap](workflow-platform-agent-roadmap.md)，
准确 PR/SHA 为 [集成清单](current-delivery-integration.md)，新增契约见
[设计](../architecture/catalog-agent-viewer-boundaries.md)。
这份计划不表示本轮已经启动 coding、下载、Gate 或合并。
它已在 Draft #213 中，该 PR 仍未合并。

## 通用约定

- Codex 负责设计、精确提交的独立审阅和集成裁决；Cursor 负责实现、测试、修复及交付。
- 优先继续现有 owner session，避免另开同一职责的 backend/frontend。
  新工作线使用下面的稳定命名，完成既有闭环后才扩大范围。
- 每个 owner 可以使用多个 subagent：先做不重叠的探索/测试分析/审查，或明确分配
  不相交文件。主 session 写公共接线、迁移、生成物并整合验证。
  不允许两个 agent 同时修改同一个 worktree/文件或各自重跑同一科学任务。
- 新工作从届时最新已接受基线建独立 worktree；既有 PR 在原分支追加修复。
  先报告 HEAD、base、工作区状态和 owner，再写入；不移动根工作区或覆盖已有修改。
- 每个 session 登记一个有限 goal，写明消费者结果、依赖、非目标和退出条件。
  goal complete 只表示该目标完成，不代表 A1/生产/科学验收全部完成。
- 交接采用准确 SHA 和 API/投影版本，不用“最新”代替身份。
  生成客户端由后端公共契约 owner 用正式入口生成，前端只消费和核对。
- 提交时保留人类作者，包含实际参与的 Cursor 工具默认署名；Codex 参与设计/审阅也
  计入，必要时使用 `Co-Authored-By: Codex <noreply@openai.com>`。不加入未参与者。
- 未明确授权的合并、push、发布不因这份计划自动获准；原会话已经获得的授权继续有效。

## A. 现有 session：现在继续

| session 名称 | 范围／写入所有权 | 依赖 | 下一步 | 验收标准／交接 | PR |
| --- | --- | --- | --- | --- | --- |
| `Cursor-A1-Backend` | 现有 C4 worktree 的 native 创建线程生命周期 | `d2a0cc…` 加未提交接线。保存、回读和真实 snapshot 已成功；提案被 SIGKILL 阻断 | 让实际 Popen 的创建线程覆盖 native 生命周期，再完成采用、拒绝、解释和重新登录回执。不提高内存，不删除 `--die-with-parent` | 原 JUnit 仍是 2 passed、1 failed、2 skipped。修复后另交集成 Draft。C4 未完成 | 尚无集成 Draft |
| `Cursor-A1-Frontend` | 不新开 C3；配合时不手改 generated | `a566943…` 继续可消费 | 没有新的页面缺陷。不另跑一套 C4 | 不把可消费写成 C4 或整个 A1 已通过 | #214 |
| `Cursor-Science-OnePath` | 不再改 #218 的已复审文档增量 | `ca40bce…` 单行文档复审通过，生产和测试代码未变 | 保留 `0fef7a2…` 的 169 项真实执行身份。不重跑科学层，不重复合并 #216 | 不把文档增量写成新的科学执行 | #218 |
| `Cursor-WashU-W1` | 取消预检后的失败状态；不扩大成 V2 | `db3d8a5…` 的生产外发 P1 已通过。新 P2 待修 | 区分主动 AbortError 和真实读取失败，使成功重选恢复。保留生产 CSP | 交错回归和最小组件测试通过；真实失败仍显示失败。不是 Viewer 全部通过 | #217 |
| `Cursor-Docs-Delivery`（现有维护 session） | roadmap、集成清单及当前 session 计划 | 当前入口已在 Draft #213 | 只做有界状态更新；不开新功能或第三批测试精简 | 历史副本逐字保留；Draft 未合并；链接、身份、空白检查通过 | #213 |
| `Codex-Delivery-Review` | 只读复审、契约裁决和交接清单 | 各 owner 的准确提交与复现 | 区分阻塞缺陷、已满足契约和后续扩展 | 每条 finding 有代码/复现/影响/最小修复边界；不把新产品目标当修复循环 | 上述 PR |

共五个 Cursor owner，加 Codex 主会话为六个。`ff4bdd5…` 只增加文案和署名，不扩大 `a08aa62…` 的摘要绑定结论。
#214 没有新的 C3 任务。C4 已经尝试，目前未完成。
冻结只表示可消费的接口版本，不代替最终实现审阅或 C4。遇到新安全/一致性缺陷继续修；
遇到新功能需求登记后续，不把它变成当前修复的无限追加项。

## B. 当前批次集成：准备好后启动

C4 仍由 `Cursor-A1-Backend` 在现有集成 worktree 继续。保留 `d2a0cc…` 和未提交接线，
不另起一份产品实现。Frontend 没有新的 C3 任务。Viewer 的 P2 不并入这个 worktree。

1. 核对 main、依赖 DAG 和来源提交，#214 已含的 #212 祖先只纳入一次。
2. 再生 OpenAPI/客户端/前端打包资产只能有一个 owner；报告合并冲突的行为处理。
3. 手工保存、刷新回读和真实 snapshot 已经成功。提案及后续步骤被 native SIGKILL 阻断。
   原 binary `60144bbe…` 的缺失和历史 A0 记录保留；`5cc82caf…` 不恢复旧身份。
   生命周期修复后，继续完成尚未执行的 C4 步骤：
   选择项目/样本/输入修订/参考，保存和回读，显式连接，读工具，提案/解释，采用/拒绝，
   真实 validator snapshot，重新登录发现回执；包含拒绝、冲突、丢响应及断连清理。
4. 没有真实 snapshot、native 或工具结果时不得用 fixture 成功文案代替。
   此层不引入真实 provider，不执行 managed run 或真实 FASTQ。
5. 最终组合的普通 CI 通过后，根据准确闭包运行需要的 Protected Bulk Gate。
   #212 已改闭包，不能拿 #210 的 Gate 代替。配置准入沿用已有授权，无授权的管理变更
   先给出具体缺口；不为了拿到绿灯降低保护。
6. #217 的生产外发 P1 已通过，取消预检后的 P2 仍独立修复。未修好时不阻塞 A1 准备，
   也不作为已验收 Viewer 进入组合。
   导航和打包资产交给集成 owner，按影响补测。
   最终交可审阅 Draft；合并、发布和关闭旧 PR 由主会话另行安排。

## C. 后续新增 session

| session／阶段 | 范围与写入边界 | 依赖 | 下一步 | 验收标准 | PR |
| --- | --- | --- | --- | --- | --- |
| `Codex-Catalog-Design`／M0 | 目录字段、修订、关联和权限契约；只读代码/设计文档 | 本轮设计、既有 data foundation | 与 M1 owner 对齐 schema/API/迁移拆分，确认 unresolved 与 run input scope 语义 | 不产生第二套身份/科学规则；字段、权限、查询和迁移可实施 | 设计随文档交付 |
| `Cursor-Catalog-M1` | 新目录领域/服务/repository/API/迁移及生成客户端 | M0 审阅；#212 共享后端接点稳定，推荐从其集成后 main 开始 | 先提交字段/接口契约；实现标准分类、实验修订和手工目录查询 | 升级不改旧修订，真实查询分页/索引、并发/权限和 provenance 边界通过 | 新 Draft |
| `Cursor-Catalog-UI` | 目录管理、组合筛选、管理员手工编辑页面 | M1 冻结契约；公共前端接点稳定 | 消费正式生成客户端，先完成无 AI 的使用路径 | 未知/冲突清楚；成员禁写新分类；管理员保存回读；真实 HTTP 浏览器闭环 | 新 Draft，必要时堆叠 M1 |
| `Cursor-Viewer-V2` | viewer 独立入口/构建/依赖和宿主薄接线 | #217 验收；轨道契约审阅 | 同仓同源 iframe/full-page，限定版本化消息 | 独立构建与样式隔离、真实 Range 绘图、worker/Peri 不参与读取 | 新 Draft |
| `Cursor-Peri-Boundary`／P1 | 专用 Peri 适配层及现有 authoring consumer 的窄接口；不新增通用 Agent 框架 | A1 C4 已完成，避免与 #212 重叠。下面的归属是设计要求，当前尚未实现 | 业务 API 保持 HelixWeave 的身份、项目范围、草稿、确认、回执和任务生命周期；Peri 启动参数、扩展协议、版本和运行身份留在适配层 | 更换或升级 runtime 不要求重写业务 API，但新实现单独验收。ACP 兼容只是接入条件；权限、工具、外发、取消和资源隔离仍逐项核验 | 新 Draft |
| `Cursor-Peri-Provider`／P2 | 模型出口/配置及隔离测试；不改目录权限 | P1 稳定的接缝、具体 provider/API 和 operator 凭据 | 分开资格 stub 与真实 API；限定外发后做最小联网验证 | 真实端点的投影/工具/撤权/超时证据；无本地模型要求；用量与来源可追溯 | 新 Draft |
| `Cursor-Catalog-AI`／M2 | 目录投影、提案/批量采用服务、审阅 UI 的明示子范围 | M1、P2、工具/外发设计复审 | 成员查询/提建议，管理员批量确认，独立工具资格 | 精确差异/来源/版本/期限，事务原子性、撤权、重取回执、无 SQL 能力 | 新 Draft，按消费者分批 |
| `Cursor-Catalog-CLI` | 终端客户端与同一公开 API 的契约测试 | M2 提案/预览/采用契约冻结 | 现有平台登录身份，显示差异并确认同一摘要 | 网页/CLI 结果与回执一致，不直写数据库、不绕过 CSRF/角色 | 新 Draft |
| `Codex-Run-Actions-Design` | create/preflight/start/cancel 的确认与权限设计 | A1、P2、既有运行生命周期 | 明确一次确认覆盖范围、幂等性、恢复和取消确认 | 不复制造任务/队列 API；设计通过后才开实现 session | 后续独立决定 |

M1 原则上按“分类与手工查询”和“可靠的结果接线”拆成两个消费者交付。
实验/assay 到运行的精确关联若需要新 snapshot binding，应单独提交契约，
不能用同一样本的偶然关联绕过。V2 不依赖 M1 完整实施，可复用现有 run-scoped 轨道目录。
P2 与 M1 可在各自前提满足后并行，M2 等二者。真实 provider 参数未确定时先完成设计，
不能把联网失败改为 stub 成功然后宣布交付。

## D. 文档、维护与交接

`Cursor-Docs-Delivery` 继续维护 Draft #213 的当前入口，只处理状态、链接和准确证据。
原根工作区单路径决策已在历史副本保留。不要在科学或 skills 分支重复写 roadmap。
#194/#195 先复核内容是否仍适用，不能整批合并过时环境/代理流程建议。

skills 随 #216 的真实科学接口变化更新；后续只按使用反馈迭代。
测试精简不自动开启第三批。部署/发布、更多科学 confirmation、生产容量、
Peri 从零构建、长期 fork、自由文本和第二 viewer 都是独立后续目标。

每次交接至少包含：

- owner、worktree、准确 base/head、工作区状态、PR；不得覆盖别人已有修改。
- 已完成消费者路径、未完成边界、公共 schema/投影/生成客户端身份及兼容变化。
- 精确测试证据：本地/CI、run/attempt、PR head 与实际 checkout、success/failure/skipped。
- 影响的科学/执行闭包与 Gate 裁决；模型 stub、真实 native、真实 provider 分别报告。
- 本任务资源清理；仅清理自己启动的进程、容器和临时文件。
- 下一个 owner 可以直接消费的契约与仍需审阅项，不用“全部完成”概括多个阶段。

文档已提交并 push 到 Draft #213，该 PR 仍未合并。后续状态更新继续使用这条分支。
实施 goal 仍按各自依赖单独开启；当前文档可审阅并不表示开发计划已经完成。

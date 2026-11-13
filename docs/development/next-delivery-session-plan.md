# 下一批交付：session 分工与验收交接

日期：2026-11-13。排期入口为 [roadmap](workflow-platform-agent-roadmap.md)，
准确 PR/SHA 为 [集成清单](current-delivery-integration.md)，新增契约见
[设计](../architecture/catalog-agent-viewer-boundaries.md)。
这份计划用于之后逐个发给 session，不表示本轮已经启动 coding、下载、Gate 或合并。

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
| `Cursor-A1-Backend` | 现有 authoring、受管目录/校验、后端测试、其 OpenAPI/client 生成 | `e00130e…` 最新边界定点复审已通过 | 交出目录及样本计划契约，支持前端/C4 集成中的具体问题 | 本轮 41 项检查通过；不加入元数据分类/真实模型/运行操作，不把 C3 可消费写成整个 A1 通过 | #212 |
| `Cursor-A1-Frontend` | 草稿页、样本行、回执 UI 及相关测试；不手改 generated | `e00130e2d5e55b3197624ef6bff9c25fd687cf87` 可供 C3 消费 | 核对祖先后纳入该契约；接目录→保存→刷新→真实 snapshot；最终新 head 跑 CI | 保留本地编辑/完整基线/操作代次；对照修订在计划行内；桌面移动端真实 HTTP；不得自动连接或启动运行 | #214 |
| `Cursor-Science-OnePath` | validator、科学 rules、相关身份/schema/docs/skills | `c44c4d6…` Picard 返修定点复审通过 | 支持单路径整体合入审阅与最终 schema 交接；不扩大科学范围 | 本轮 27 项定向检查通过；历史记录不改，默认不宣称生物学冻结 | #216 |
| `Cursor-WashU-W1` | 当前轨道服务/转换/页面及测试 | 最新 `bcdbd325…` 独立审阅 | 回答实际读取/绘图及身份拒绝问题；维持当前消费者范围 | 真组件读取并展示已知数据，有限资源与清理；不假称独立前端或支持所有接触矩阵 | #217 |
| `Cursor-Docs-Delivery`（现有维护 session） | roadmap、集成清单、架构设计及文档检查 | #215 已合入；本轮文档工作区和最新交接 | 整合本轮文档，与旧 #213 协调一次纳入；不开新功能或第三批测试精简 | 历史逐字保留，当前状态只一处权威入口；链接、身份、空白检查通过 | #213 或协调后的单一文档 PR |
| `Codex-Delivery-Review` | 只读复审、契约裁决和交接清单 | 各 owner 的准确提交与复现 | 区分阻塞缺陷、已满足契约和后续扩展 | 每条 finding 有代码/复现/影响/最小修复边界；不把新产品目标当修复循环 | 上述 PR |

共五个 Cursor owner，加 Codex 主会话为六个。后端本轮复审已完成，前端无需继续
等待 e00130e 的这两个修复，可开始正式接线；后续新字段仍需单独审阅。
冻结只表示可消费的接口版本，不代替最终实现审阅或 C4。遇到新安全/一致性缺陷继续修；
遇到新功能需求登记后续，不把它变成当前修复的无限追加项。

## B. 当前批次集成：准备好后启动

`Cursor-A1-Integration` 使用独立 worktree，消费已审阅的后端、前端和最终科学 schema。
不让 backend 或 frontend owner 在未告知另一方时更换接口或重写共同历史。

1. 核对 main、依赖 DAG 和来源提交，#214 已含的 #212 祖先只纳入一次。
2. 再生 OpenAPI/客户端/前端打包资产只能有一个 owner；报告合并冲突的行为处理。
3. 用隔离数据库、真实 API、真实浏览器、固定 native 与本地模型 stub 执行 C4：
   选择项目/样本/输入修订/参考，保存和回读，显式连接，读工具，提案/解释，采用/拒绝，
   真实 validator snapshot，重新登录发现回执；包含拒绝、冲突、丢响应及断连清理。
4. 没有真实 snapshot、native 或工具结果时不得用 fixture 成功文案代替。
   此层不引入真实 provider，不执行 managed run 或真实 FASTQ。
5. 最终组合的普通 CI 通过后，根据准确闭包运行需要的 Protected Bulk Gate。
   #212 已改闭包，不能拿 #210 的 Gate 代替。配置准入沿用已有授权，无授权的管理变更
   先给出具体缺口；不为了拿到绿灯降低保护。
6. #217 可独立验收再进入组合；导航和打包资产交给集成 owner，按影响补测。
   最终交可审阅 Draft；合并、发布和关闭旧 PR 由主会话另行安排。

## C. 后续新增 session

| session／阶段 | 范围与写入边界 | 依赖 | 下一步 | 验收标准 | PR |
| --- | --- | --- | --- | --- | --- |
| `Codex-Catalog-Design`／M0 | 目录字段、修订、关联和权限契约；只读代码/设计文档 | 本轮设计、既有 data foundation | 与 M1 owner 对齐 schema/API/迁移拆分，确认 unresolved 与 run input scope 语义 | 不产生第二套身份/科学规则；字段、权限、查询和迁移可实施 | 设计随文档交付 |
| `Cursor-Catalog-M1` | 新目录领域/服务/repository/API/迁移及生成客户端 | M0 审阅；#212 共享后端接点稳定，推荐从其集成后 main 开始 | 先提交字段/接口契约；实现标准分类、实验修订和手工目录查询 | 升级不改旧修订，真实查询分页/索引、并发/权限和 provenance 边界通过 | 新 Draft |
| `Cursor-Catalog-UI` | 目录管理、组合筛选、管理员手工编辑页面 | M1 冻结契约；公共前端接点稳定 | 消费正式生成客户端，先完成无 AI 的使用路径 | 未知/冲突清楚；成员禁写新分类；管理员保存回读；真实 HTTP 浏览器闭环 | 新 Draft，必要时堆叠 M1 |
| `Cursor-Viewer-V2` | viewer 独立入口/构建/依赖和宿主薄接线 | #217 验收；轨道契约审阅 | 同仓同源 iframe/full-page，限定版本化消息 | 独立构建与样式隔离、真实 Range 绘图、worker/Peri 不参与读取 | 新 Draft |
| `Cursor-Peri-Boundary`／P1 | authoring 到 adapter/repository 的窄接口、实际接线及测试 | A1 C4 已完成，避免与 #212 重叠 | 移走通用服务的 ENCODE 直接依赖，保留 owner/SQLite 权威 | 同一真实消费者走新接口，A1 行为/拒绝路径不变；无泛化框架 | 新 Draft |
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

`Cursor-Docs-Delivery` 协调本轮文档与 #213；只处理当前入口、链接和准确证据。
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

本轮文档 worktree 仅用于交付设计；没有提交、push 或新远端分支。
后续接手者需先整合这组文档到自己的文档分支，核对最新状态后再开相应实施 goal。

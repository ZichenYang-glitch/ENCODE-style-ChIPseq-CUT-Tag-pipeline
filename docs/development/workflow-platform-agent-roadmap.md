# HelixWeave Product Roadmap

## 2026-11-25 实验网页取回历史轨道

首轮试用包含「实验记录 → 结果第一版 → 在 WashU 查看」。用两份同名、同大小且信号不同的轨道验收，第二版存在时须能从网页重开第一版。操作员在终端处理长期盘上的文件。[历史版本设计与验收](experiment-archive-track-history.md)

## 2026-11-19 当前进度

以下按主干 `6e29d2d` 核对。合入状态和部署到具体机器的状态分别记录。

| 方向 | 主干交付 | 下一次验收 |
| --- | --- | --- |
| WashU | 独立同源页面内可以选择、加载和移除当前运行的轨道。[#222](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/222)、[#237](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/237)。早期草案 [#227](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/227) 的验收用例已纳入实现契约。 | 与目录网页合成同一份前端包后，从实验入口打开实际选中的轨道。 |
| Peri 与 adapter | 运行时会话已分离，离线建议传输和 Hi-TrAC mapq 的 adapter 审阅已合入。[#223](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/223)、[#232](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/232)、[#234](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/234) | 已保存草稿的实验校验入口待真实用户试用。[#242](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/242)。供应商外发及 Agent 写操作先审设计。 |
| 实验目录 | M0、目录基础、验证快照绑定及历史列表两种视图已合入。[#224](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/224)、[#225](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/225)、[#236](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/236)、[#241](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/241)。目录 CLI 可只读预览分类导入。[#230](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/230) | 目录网页、分类建议预览与采用、按结果反查分别验收。[#235](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/235)、[#229](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/229)、[#240](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/240)、[#243](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/243) |
| 存储与科学结果 | 外部结果登记读取契约、单份结果恢复演练、三种工作流的科学验收证据台账已合入。[#228](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/228)、[#233](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/233)、[#231](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/231) | 合成环境验证迁盘后旧链接读取、外部结果按版本下载及组合前端。[#238](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/238)、[#239](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/239)。实验室归档 SOP 见 [#244](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/244)。 |

**外部生信分析结果归档。** 实验室用其他流程或工具得到的结果也纳入长期归档。按课题、样本或细胞、实验整理；文件保存在长期存储目录，目录数据库记录来源、版本和存储位置。验收时，应能从实验找到指定版本并取回文件，无需为该数据类型实现查看器。Hi-C 矩阵、单细胞分析结果都是归档示例；[#248](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/248) 只记录可选的查看场景，有具体需求再单独立项。

**下一阶段设计。** Agent 创建与取消运行见 [#245](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/245)，本地试用到服务器发布的差距见 [#246](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/246)，Peri 首次获准外发见 [#247](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/247)。11_17 的试用单由 IRENE-34 准备，实际样本和同事的试用反馈是另一步。涉及供应商、预算、允许外发的字段和正式部署的选择留给项目负责人。

架构解耦的后续范围见下文「P1b 后续维护」。各里程碑以其主干实现和验收证据分别收口。

## 快照：2026-11-17T08:13:19Z

2026-11-16T08:22:10Z 及更早快照保留在后文，不倒写。
本次查询 origin/main 为 `ff3f7d981592a844750dba502f980a356becc874`，
即 [#220](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/220)
的普通 merge。父提交是 roadmap 合入 `b21d941…` 与 #220 head `a299db9…`。
本文件将由新的文档 Draft 交付，不写回本提交 SHA。根工作区未提交的 roadmap 未搬入。

已合入且不再是当前阻塞：

- A1 在 `1b1db4…`。Gate 与 C4 的实际运行 SHA 仍以 08:22:10Z 入口为准，不改写成在后续 SHA 上新跑。
- Viewer 在 `2a6c22a…`，批准 head `04e5075…`。
- roadmap 文档在 `b21d941…`。
- #220 合入后的 push attempt 1 全部 success：
  [CI 44039567042](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/44039567042)、
  [Lint 44039567053](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/44039567053)、
  [Lock Check 44039567123](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/44039567123)。
  coverage job `133900520618` success，完整覆盖率门槛通过。
  四个条件执行 job 为 skipped，不是执行通过。没有派发新的 Bulk Gate。

下一轮分工。下面未列出的实现都还没开始。ACP 兼容不是自动获得工具、外发或运行资格。
新增目录整理的管理员规则不改变已有 A1 草稿采用权限。

| 会话 | 范围 | 状态 |
| --- | --- | --- |
| `HW-Agent-P1-Boundary` | 先做 P1a runtime 边界 | 尚未实施 |
| P1b | ENCODE authoring 投影与窄 repository 接口 | 后续单独收口，不并进 P1a |
| P2 | 真实外部模型 API。依赖稳定接点和明确 provider 配置 | 设计与实施尚未开始。不要求本地模型 |
| `HW-Viewer-V2-Isolation` | 独立构建、同源嵌入与全页 | 尚未实施 |
| Frontend | 配合 Viewer | 不新开重叠实现 |
| Science | 无新增任务 | 不启动新确认批次 |
| 数据目录 M0/M1、AI 整理、Agent 运行操作 | 各自契约和验收 | 仍待开始 |

非阻塞维护，不因此自动开新测试批次：#220 的全进程 FD 断言隔离建议，
以及 relay 第三连接断言补强。两者都不改生产行为，也不改写已有 C4 或 Gate 身份。

### P1b 后续维护：实验配置服务与工作流适配器解耦

`services/authoring.py` 直接调用 ENCODE 的 `project_encode_fields()` 和
`normalize_encode_workflow_sample_plans()`。在 P1b 中将字段投影和样本组织规则交给对应的
workflow adapter，由公共服务通过适配器接口调用。草稿保存、版本检查和提交仍由公共服务负责。

验收时，已有 ENCODE 草稿应能正常回读、校验和提交。选取另一个工作流的真实配置需求
（例如 HiTrAC 的专属字段），验证改动能集中在该适配器及其测试中，公共服务无需增加工作流专用分支。

## 快照：2026-11-16T08:22:10Z

2026-11-15T12:18:06Z 及更早快照保留在后文，不倒写。
本次查询 main 为 `1b1db421ba3172bc39cf31ce0d03976511379fda`，
即 [#219](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/219)
的普通 merge。父提交是 `101983f…` 与批准 head `e112223…`。
两者文件树同为 `1f9ac53f6b34625fcbb6c09965427a3d6867e1ef`。
本文件在 Draft #213，仍未合并，不写回本提交 SHA。根工作区既有修改未动。

**受限配置建议与校验 A1 已交付并进入 main。**
C4、launcher 生命周期和最终 Bulk Gate 不再是当前阻塞。
它们的证据分别属于下面的准确 SHA，不能写成在 merge SHA 上重新跑过。

| 证据 | 实际运行身份 |
| --- | --- |
| Bulk Gate | `e112223…`，[CI 43792685991](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43792685991) attempt 2 success。job `133061548271` 通过 4 项执行测试和 5 项产品浏览器测试。Artifact `14195706467` 的 ZIP SHA256 为 `54ffa1c4aa24454ffc0dc63e8d6a0728004345c977c2633b0e05b5985adb88a6` |
| C4 | `82d091278a5cf46f5c33b40f4cbb26d3d12fe2a8`。Vite 页面、真实平台、固定 Peri、模型 API stub，5 passed |
| 合入后 push | [CI 43816767856](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43816767856)、[Lint 43816767900](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43816767900)、[Lock Check 43816767895](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43816767895) 均为 attempt 1 success，head 是 merge SHA |

相同文件树支持把 `e112223…` 的 Gate 和 `82d0912…` 的 C4 复用到这次 merge。
push CI 的 real-execution、platform-real-execution、container-smoke、
bulk-rnaseq-real-execution 为 skipped，不是执行通过。
本次查询核对了 run、attempt、head 和 job 结论，没有重新下载 Gate artifact。

GitHub 已把 #212、#214、#218 标为 merged，其已审阅内容包含在这个 main 中。
#216 的内容此前已在 `101983f…`；本次查询它也被标为 merged，记录的 merge 提交仍是那个较早 SHA。
不需要再次合并。

当前进行中：

- #217 仍是 open Draft，head 仍为 `c4c4909…`。它独立推进，与新 main 的组合仍需核验。
- #213 本身仍是文档 Draft，未合并。

后续路线。设计尚未实施；ACP 兼容只是接入条件，不是自动通过接入。
现有 A1 草稿采用权限不因下面的目录或 AI 设计而改变。

| 后续 | 状态 |
| --- | --- |
| 真实模型 API | 设计与实施都还没开始。不要求本地模型 |
| Agent runtime 适配 | 设计要求：平台拥有授权、项目范围、草稿、确认、回执和生命周期；Peri 协议与版本留在适配层。尚未实施 |
| Agent 提交与取消运行 | 独立操作契约、确认和幂等仍是设计，尚未开始 |
| 数据目录 | 按课题、细胞类型、实验类型分类与查询仍是设计，尚未实施 |
| AI 整理 | 设计：成员查询和建议，管理员在网页或终端确认。不改变已有 A1 草稿采用权限 |
| managed execution、生产容量、科学确认 | 另列，不属于这次 A1 交付 |

非阻塞注记：`test_relay_refuses_a_third_connection_and_closes` 仍只证明某个 handler 启动。
可用两个 handler 屏障加第三连接结束和计数断言补强，不需要改生产 relay。
若以后另交测试提交，保留 `82d0912…` 的 C4 身份，不整批改写。

## 快照：2026-11-15T12:18:06Z

2026-11-15T07:21:49Z（#213 `23681d9…`）及更早快照保留在后文，不倒写。
本次查询 main 仍为 `101983f…`。#216 代码已在 main，PR 对象仍 open，不因此再次合并。
本文件在 Draft #213，仍未合并，不写回本提交 SHA。

集成 Draft [#219](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/219)
的已审阅 head 是 `6d570822620ccb873c26a1dabad1a8baebe5593b`。
当前不再把 `d2a0cc…` 加未提交接线、没有集成 Draft 当作现状。
短命创建线程已由常驻 launcher 修复。新的待返修是创建超时后晚到进程可能无人接管；
不能靠加长超时掩盖。

C4 本地记录有 5 passed，完整验收未完成：没有采用后再校验；丢响应测试实际阻止请求发送；
未执行 Reject；所谓过期保存实际是 revision 冲突。原先 2 passed、1 failed、2 skipped 仍是更早历史证据。
pull_request [CI 43693526422](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43693526422)
attempt 1 success。同一 head 的 workflow_dispatch
[43693536449](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43693536449)
attempt 1 failure：Bulk 在 runner 执行前被分支保护拒绝，coverage 为 85.37%，低于 87.28%。
不能写成只剩环境准入。

#217 `c4c4909…` 的取消恢复增量已通过独立定点复审。`db3d8a5…` 的 P2 不再列为未修复。
该准确提交的生产浏览器证据仍待收口。#214 `a566943…` 与 #218 `ca40bce…` 已纳入 #219。
Frontend 待命，Science 没有新增任务。169 项真实执行证据仍属于 `0fef7a2…`。
C4、整个 A1、最终组合 Gate、真实 provider 和 Agent 运行操作仍分别记账。

`23681d9…` 的 CI 里 PID 123 的 SIGKILL 是测试 mock 记录的 `os.kill` 调用，
不是已证明向宿主 PID 123 实际发信号。

| 工作线 | 下一步 |
| --- | --- |
| #219 | 返修创建超时后的进程归属，并补齐尚未执行的 C4 场景。覆盖率按未覆盖行为补齐，不降低门槛 |
| Viewer | 收口 `c4c4909…` 的生产浏览器证据。不把取消恢复复审写成 Viewer 全部通过 |
| 前端 / 科学 | Frontend 待命。Science 无新增任务，不重跑 `0fef7a2…` 的 169 项 |

## 快照：2026-11-15T07:21:49Z

2026-11-14T19:18:06Z 及更早快照保留在后文和集成清单，不把后来的结果写进那些段落。
本次查询 main 仍为 `101983f…`。#216 代码已在 main，PR 对象仍 open，不因此再次合并。
本文件已在 Draft #213，仍未合并，也不写回本提交 SHA。

C4 已尝试，未完成。手工保存、刷新回读和真实 snapshot 成功。
native 在短命请求线程中创建，线程正常退出后触发 `bwrap --die-with-parent`，native 收到 SIGKILL。
这是创建线程生命周期缺陷，不是 OOM 或资源不足。原 JUnit 为 2 passed、1 failed、2 skipped。
集成 HEAD 是 `d2a0cc683aefd3a05fe507dce0fba60e79ecc320` 加未提交接线；交接时没有集成 Draft。
下一步只修 owner 线程生命周期，然后回到该 worktree 完成采用、拒绝、解释和重新登录回执。

#217 `db3d8a5…` 的生产外发 P1 已通过复审。新 P2：取消预检后再选择不能恢复，等待窄修复。
这不是 Viewer 全部通过。#214 `a566943…` 继续可消费，没有新的 C3 任务。
#212 `ff4bdd5…` 的文案与署名增量、#218 `ca40bce…` 的单行文档增量已通过复审，生产和测试代码未变。
169 项真实执行证据仍属于 `0fef7a2…`。
C4、整个 A1、最终组合 Gate、真实 provider 和 Agent 运行操作仍分别记账。
某个 session 的 goal complete 不表示这些全部完成。

| 工作线 | 下一步 |
| --- | --- |
| C4 | Backend 在现有集成 worktree 修复 native 创建线程生命周期，再完成被阻断的提案及后续步骤。不提高内存、不删除 `--die-with-parent` |
| Viewer | 只修取消预检后的失败状态残留，保留已通过的生产 CSP。不扩大成 V2 |
| 前端 / 科学 | #214 保持可消费。#218 不重跑科学层，不改写 `0fef7a2…` 的测试身份 |

## 快照：2026-11-14T19:18:06Z

本文件记录当前产品边界、优先级和交付顺序。准确 PR/SHA/CI 快照见
[集成清单](current-delivery-integration.md)，实施分工见
[session 计划](next-delivery-session-plan.md)，新设计见
[数据目录、Agent 与 Viewer 边界](../architecture/catalog-agent-viewer-boundaries.md)。
本次只更新当前状态，没有修改产品、运行迁移、新的科学数据批次、Gate 或发布。

2026-11-14T17:40:08Z 的入口保留为集成清单中的同名快照，不在这里改写。
2026-11-14T19:18:06Z 查询，main 仍为
`101983f9480bb72306da69875e79dc6ad65bca5a`。#216 的代码已在该 merge 中，
PR 对象仍 open、未标记 merged；这个异常不表示需要再次合并。
本轮只采纳各 owner 交回的最终 head，不追逐其间的中间提交。
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

#214 `a566943…` 的保存返修已通过复审和普通 CI，可进入 C4 准备。这不是整个 A1 批准。
`399c3b3…` 上的返修要求留在该 SHA。
#218 `0fef7a2…` 的科学真实执行已通过：169 passed、零失败、零错误、零跳过。
workflow `43628729299` 的整次结论是 success；bulk job 为 skipped，不是执行通过。
`101983f…` 上原来的 dispatch 失败仍只属于那次 run。禁止再次合并 #216。
#212 现为 `a08aa62…`。相对 `e00130e…` 的窄增量是绑定已接受的 `5cc82caf…` 摘要，
该增量复审通过，普通 CI 也已 success。`e00130e…` 的边界复审不自动扩大到别的提交。
#217 `aa0ef9f…` 的 gzip 与空白返修已通过，但新增 P1：生产读取失败会尝试外部 CORS 代理。
这不是 Viewer 全部通过。Viewer 返修独立推进，不阻塞 A1 准备，也不能当作已验收 Viewer 纳入交付。
C4、整个 A1、真实 provider、Agent 运行操作和最终组合 Gate 仍分别记账，本轮都未执行。

原 binary `60144bbe…` 的缺失和历史 A0 验收记录保留，不改写成失败，也不转给新摘要。
`5cc82caf…` 已在固定受限范围内被接受为替代运行身份。该接受不是整个 A1、C4、
真实 provider 或生产容量结论，也没有恢复旧 binary 身份。

## 后续交付计划

每项记录范围、依赖、下一步、验收标准和关联 PR。细节以设计及 session 计划为准。

| 工作线／范围 | 依赖 | 下一步 | 验收标准 | 关联 PR |
| --- | --- | --- | --- | --- |
| A1：受限配置建议与校验 | #212 `a08aa62…`、#214 `a566943…`、#218 `0fef7a2…` 均已有本轮最终复审 | Backend 在独立 worktree 主导 C4 准备，Frontend 配合。不加入新分类系统 | C4 本身仍待执行：真实 HTTP、snapshot、回执和已接受的 `5cc82caf…` runtime。整个 A1 与最终 Gate 另计 | #212、#214、#218；集成 PR 待创建 |
| 科学单路径 | #218 `0fef7a2…` 的整次 workflow `43628729299` success；真实执行 169 passed | 保留该 head。非阻塞的 sample-sheet 文案可另改，不因此重跑科学层。禁止重复合并 #216 | 代码与测试身份不改写。不是生物学确认，也不是 Gate | #218 |
| WashU 首版 | `aa0ef9f…` 的 gzip/空白返修通过；生产 CORS 代理回退仍是 P1 | 在现有生产页面限制受保护轨道的外部代理，并用生产约束证明失败路径没有外部请求 | 正例仍可读；失败不发外部请求、不显示部分成功。不是 V2，也不是 Viewer 全部通过 | #217 |
| 数据目录 M1 | 本轮设计；A1 共享 API/迁移接点稳定 | 字段契约复核后开发课题/样本/实验分类与查询 | 手工管理先可用；组合筛选/分页准确；旧修订不改写；不夸大结果归属 | 独立 PR 待创建 |
| Viewer 隔离 V2 | WashU 首版通过；版本化轨道契约 | 同仓独立构建与样式/依赖隔离，保持同源认证 | 独立构建可嵌入；不直接查 SQL/调执行；worker/Peri 停止后仍可读已有结果 | 独立 PR 待创建 |
| Peri 收口 P1／真实模型 P2 | `5cc82caf…` 只在固定受限范围被接受；旧 binary 记录保留。业务适配层尚未实现 | C4 可使用该已接受摘要。真实 provider 和 Agent 运行操作仍不开始 | 更换 runtime 不改业务 API，但新实现单独验收。不新增通用 Agent 框架 | #212 后续独立 PR |
| AI 元数据整理 M2 | M1、P2、目录投影及工具资格 | 先读与建议，再接管理员批量采用；网页/终端共用契约 | 无 SQL 工具；有来源、差异、版本、期限和回执；撤权/冲突拒绝写入 | 独立 PR 待创建 |
| Agent 运行操作 | A1、真实 provider、业务确认设计 | 设计 create/preflight/start/cancel 与既有生命周期的映射 | 不另造任务 API；确认范围、重试、取消请求/确认及审计有明确契约 | 尚未开始 |
| 文档、skills、测试与 CI 运维 | 实际反馈及已确认缺陷 | 当前入口保持一份；按具体收益处理维护 | 历史不倒改；skip 不当通过；不为省测试弱化契约 | 本状态在 Draft #213，未合并；#194/#195 待审 |

## 顺序与并行规则

1. 现有五个 Cursor session 为 Backend、Frontend、Science、Viewer 和维护。
   维护继续在 Draft #213 更新当前入口，不另开文档 PR。加 Codex 主会话共六个。
   Codex 复审准确 head；CI 成功不会自动冻结接口或批准合并。
2. #216 的代码已经在 main，禁止再次合并。C4 准备使用 #212 `a08aa62…`、
   #214 `a566943…` 和 #218 `0fef7a2…`。Backend 主导独立 worktree，Frontend 配合。
   Viewer 的 CORS 返修独立推进。集成前检查 #212 与 #214 的祖先，避免重复补丁。
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

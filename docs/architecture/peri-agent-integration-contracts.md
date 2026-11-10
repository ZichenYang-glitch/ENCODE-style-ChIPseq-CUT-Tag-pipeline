# Peri ADR 附录：接口、工具和验收交接

> ACP最终补测（2026-11-09）：原策略单测1/0/0通过，源码和运行binary身份不变；**固定版本、固定补丁和已声明边界内的A0验收通过，交回主会话审阅**。见[追加验收材料](../../scripts/peri_a0_closeout/evidence/acp-final/RESULTS.md)。A1未获准、未决定长期fork。以下保留此前各轮当时结论，不将旧“暂不通过”当作最新状态。

> 补测前收尾结论（2026-11-09，ACP 单测尚未执行）：固定受限候选 native 矩阵23/23通过，目录库回归通过；新增ACP Rust单测编译触发4GiB OOM，未执行，**当时完整A0暂不通过，不建议进入A1**。见[当时证据](../../scripts/peri_a0_closeout/evidence/RESULTS.md)。这不是补测后的当前状态。

配套[ADR](peri-agent-integration-decision.md)。状态：设计修订，2026-11-09；仅A0获实施授权；
初次研究及 A0 实验基线 HelixWeave `622778f64e59f935d822d92a91bade78807dabce`、
Peri `c7d154c9985cbe08ba25e964cb7468d15014e6d0`。
以下“拟议工具/投影”不是已实现 endpoint。下列验收表是契约目标，不能当作已通过；A0实际执行结果由原型证据单独报告。

2026-11-09 接续核验：远端 main 为 `b2f15868c6f74ab8764a73ae6fe4394f67304475`，
#210 已合入并验收。仅更新平台取消依赖，不迁移原型或改写旧实验身份。写这段时
ACP 单测尚未补上，当时把固定接法记为不通过；该句不是补测后的当前状态。当前
状态只见 [A0 归档](../development/peri-a0-archive.md)。原版未测硬门仍以
[逐项证据](../../scripts/peri_a0/evidence/RESULTS.md) 的当时范围解释。

本次“接法修订与复验”新增[候选矩阵](../../scripts/peri_a0_revision/evidence/RESULTS.md)：
采用完整catalog策略实验补丁，保留原版失败。初始工具面/一次审批拒绝已有证据，
完整ready仍被辅助prediction出口阻断；双实例、恢复、允许工具成功等门仍未测。
准入须覆盖主推理和所有辅助模型调用，不能从主请求工具目录推导全部出口已准入。
辅助能力必须有消费者、投影、会话绑定和中断契约后另行决定，不能直接加入白名单。
relay期限/连接预算及临时cgroup已补；控制器探针与native运行中故障验收分开记账。
本次不改变下述A1/B2授权与前置门，也不因局部通过降低A0验收要求。

## 1. 现有接口到代码的映射

HTTP 路径以下均带 `/api/v1` 前缀；与[已提交 OpenAPI](../../frontend/openapi.json)
逐项核对。模型不会获得这些路由的通用 HTTP 调用器。

| 能力与真实路径 | 实现/证据（相对仓库根） | 可复用与缺口 |
| --- | --- | --- |
| `POST /auth/login`、`POST /auth/logout`、`GET /auth/session` | [dependencies.py](../../src/encode_pipeline/api/dependencies.py)；[authentication_service.py](../../src/encode_pipeline/services/authentication_service.py)；[auth tests](../../test/api/test_auth_routes.py) | cookie/CSRF、即时账号/session撤销；没有 Agent bearer token或project ACL |
| 所有业务 router 的 member gate | [routes/__init__.py](../../src/encode_pipeline/api/routes/__init__.py)；[认证 ADR](lan-authentication-decision.md) | 用户共享项目；Agent另加服务端scope检查，不给管理员额外工具 |
| 项目、样本、输入引用 | [ValidationRequest](../../src/encode_pipeline/api/models.py)；[data_registry](../../src/encode_pipeline/platform/data_registry.py)；[绑定测试](../../test/services/test_project_sample_snapshot_bindings.py) | 有内部registry/binding；无project/sample/input列表路由，不能发明 `/projects`；MANAGED_REVISION_V1 handoff当前拒绝，不能绕过资格门 |
| `GET /workflows/`、`GET /workflows/{workflow_id}` | [workflows routes](../../src/encode_pipeline/api/routes/workflows.py)；[workflow_info](../../src/encode_pipeline/services/workflow_info.py) | descriptor/capabilities/authoring及execution availability |
| `GET /workflows/{workflow_id}/schema`、`.../reference-profiles` | 同上；[reference_profiles service](../../src/encode_pipeline/services/reference_profiles.py) | adapter拥有schema/version/限制，enabled reference revisions；模型面需更窄字段投影 |
| `POST /workflows/{workflow_id}/validate` | [validated_inputs](../../src/encode_pipeline/services/validated_inputs.py)；[snapshots](../../src/encode_pipeline/platform/snapshots.py)；[tests](../../test/services/test_validated_inputs.py) | ValidationRequest config/samples/options+可选bindings；成功可snapshot=null；需要完整准入才有可提交snapshot |
| `POST /workflows/{workflow_id}/runs` | [runs routes](../../src/encode_pipeline/api/routes/runs.py)；`ValidatedRunCreationService`；[SQL repository](../../src/encode_pipeline/persistence/repositories.py) | snapshot_id+tags；201新建/200重放；`VALIDATED_SNAPSHOT_REPLAY_CONFLICT`；没有现成按snapshot GET |
| `POST /runs/{run_id}/preflight` | [preflight route](../../src/encode_pipeline/api/routes/preflight.py)；[service](../../src/encode_pipeline/services/preflight.py) | created→validating；BackgroundTasks规划到planned；重复409；会物化workspace/dry-run，不是A校验 |
| `POST /runs/{run_id}/start` | [RunSubmissionService](../../src/encode_pipeline/services/run_submission.py)；[tests](../../test/services/test_run_submission.py) | planned才可排队；build/admission复核；稳定job、reservation/dispatch收敛；禁止Agent直接Redis |
| `GET /runs`、`GET /runs/{run_id}`、`.../events`、`.../logs` | [RunService](../../src/encode_pipeline/services/runs.py)；[models](../../src/encode_pipeline/api/models.py) | SQLite真值、分页非SSE；history无project filter，public run无project/snapshot；必须内部归属resolver |
| `POST /runs/{run_id}/cancel` | [RunCancellationService](../../src/encode_pipeline/services/run_cancellation.py)；[jobs](../../src/encode_pipeline/workers/jobs.py)；[timeouts](../../src/encode_pipeline/workers/timeouts.py) | 接受意图≠worker ack；[race tests](../../test/services/test_run_cancellation.py)；public run暂未直接给requested/ack时间 |
| `GET /runs/{run_id}/artifacts`、`.../{artifact_id}`、`.../{artifact_id}/download` | [artifact routes](../../src/encode_pipeline/api/routes/artifacts.py)；[ArtifactDownloadService](../../src/encode_pipeline/services/artifact_downloads.py)；[tests](../../test/services/test_artifact_downloads.py) | run/artifact/generation/revision核验、descriptor安全、审计；保持下载在可信UI，不将字节送模型 |
| `GET /runs/{run_id}/qc-metrics` | [QC routes](../../src/encode_pipeline/api/routes/qc_metrics.py) | generation cursor；换代冲突从第一页重读，不能拼接新旧结果 |
| `POST /workflows/{workflow_id}/agent/chat` | [AgentService](../../src/encode_pipeline/services/agent.py)；[tools](../../src/encode_pipeline/services/agent_tools.py)；[defaults](../../src/encode_pipeline/services/defaults.py) | 实际固定explain_issues；session_id仅回显；registry validation用ValidationService而非durable snapshot；不是通用agent loop。route的mock-only文案不足以代表factory实际provider能力 |
| 审计、错误安全 | [security_audit](../../src/encode_pipeline/platform/security_audit.py)；[run audit](../../src/encode_pipeline/services/run_security_audit.py)；[Agent audit](../../src/encode_pipeline/services/agent_audit.py)；[redaction](../../src/encode_pipeline/services/agent_redaction.py) | 安全审计持久且闭合，Agent audit有界内存；缺业务确认/operation/tool-version关联。不得原样透传Peri诊断 |
| 前端草稿/提交 | [ValidatedSubmission](../../frontend/src/features/input-workbench/ValidatedSubmission.tsx)；[RunProgressPanel tests](../../frontend/src/features/run-progress/RunProgressPanel.test.tsx)；[AgentSidebar](../../frontend/src/features/agent-sidebar/AgentSidebar.tsx) | semanticRevision使旧snapshot失效；create后导航并preflight，start独立；Agent消息仅React state、session_id=null；复用既有workbench，不重造配置表单；project/sample选择控件及绑定接线仍需补齐 |

审阅了这些测试的断言；本轮没有执行它们，不把存在测试写成当前环境已通过。

## 2. 最小工具契约（拟议 `hw.*@1`）

共用 envelope：`contract_version, call_id, observation_time, ok, data, issues`。
`issues`只含白名单 code/severity/safe_message/retry_class；不带traceback/private context。
对写操作增加 `operation_id, disposition`：`rejected | accepted | outcome_unknown`。
accepted只描述该操作的已知回执；`run.status`永远来自平台。结果分页返回
`next_cursor/generation`，大小/页数/超时由服务端限制，不允许任意循环扩展数据集。

调用上下文由 API 绑定 `user_session + conversation + project_scope + generation`，
不由模型传入principal/token。每次执行入口统一校验工具版本/参数schema、当前身份、
真实对象归属和外发字段。未知工具/字段/跨scope引用拒绝，不交给反射执行。

| 工具 / 分级 | 输入 | 输出 | 权限、确认、幂等语义 |
| --- | --- | --- | --- |
| `list_workflows` / 只读 | cursor/limit | 安全descriptor、能力与availability | 当前member；全局目录，不带工作区配置 |
| `get_workflow_schema` / 只读 | workflow_id | 安全schema投影、schema/adapter version、必要本地填写提示 | schema文本/default也筛选；模型不获得私有输入路径；不承诺投影可取代真实validator |
| `list_reference_profiles` / 只读 | workflow_id | 可选revision的安全引用/适配信息 | 仅enabled compatible；展示名/说明亦过滤 |
| `get_run` / 只读 | run_id | canonical status、更新时间、安全取消意图/ack投影 | 从真实binding检查scope；无归属不能猜测；请求/ack分别报告 |
| `get_results` / 只读 | run_id, kind=artifact或qc, cursor, generation, limit | 白名单metadata/QC、安全UI链接 | scope+原结果契约；禁止自动download。generation变更显式重启分页 |
| `propose_draft` / 草稿 | workflow_id, base_draft_revision, 结构化patch | draft_id/revision、diff、安全待填项 | CAS revision；只能schema允许字段。输入/路径由本地选择引用提供，不能模型造文件路径。重复call只回同草稿回执 |
| `validate_draft` / 校验 | draft_id/revision | safe issues + snapshot安全身份或null | 用户已采用该revision；复用ValidatedInputService；它会写snapshot，**不是纯只读**；绑定同revision结果缓存，失败/丢失可重新校验但不得自动批准新snapshot |
| `submit_confirmed` / 写 | operation_id（不接受config/snapshot替换） | run_id及create/preflight/start已知进度、平台状态 | 平台批准有效且scope/身份匹配；读取服务器固定内容，由services窄coordinator推进三段；幂等依据原snapshot/稳定job；多call不新建任务 |
| `cancel_confirmed` / 写 | operation_id | target_run_id、request_accepted、observed_status、ack evidence | 独立取消批准绑定run/执行身份；重复请求复用原取消契约；不可将202显示为已停止 |
| `get_operation` / 只读 | operation_id | 安全确认摘要、原snapshot/run关联、实时平台结果、需重新确认原因 | actor/scope访问检查；过期可查历史，不能触发副作用；补现有丢回执查询缺口 |

不暴露通用 `http_request`、SQL、shell、filesystem、上传、任意MCP配置、管理员恢复、
原始日志或任意artifact正文工具。工具清单没有 `approve`：确认是UI直接到平台的
用户操作，Peri无权自批。`submit_confirmed`的Peri工具许可即使被允许，缺业务确认
仍拒绝；业务批准存在也不能绕过关闭的工具能力或撤销后的身份。

## 3. 权限装配的源码结论与阻断原型 A0

固定Peri源码追踪（链接均在ADR固定来源）：

1. `main.rs` ACP入口创建host Bypass。
2. `host/workspace.rs:374` 为每个正常SessionEnvironment创建
   `SharedPermissionMode::new(host.permission_mode.load())`，不是clone host Arc。
3. `host/requests.rs:97-105` 将非生命周期请求分派到真实session的environment.cfg；
   `config_options.rs:37,67` 写这份mode。缺session/closing拒绝；无environment有host fallback。
4. `prompt_dispatch.rs:77-82` → `prompt.rs` → `host/stage_builder.rs` →
   `peri-agent/src/session/exec/stage_builder/agent.rs` → `middlewares/assembly.rs`
   把该会话Arc传给PermissionMiddleware，实际决策读取该mode。
5. 因而正常native会话的A改mode不沿该Arc传播到B；已有
   `requests_workspace_cases_test.rs::worktree_scheduled_approval_uses_session_permission_and_rejects_closed_owner`
   只证明env/host隔离片段，不能替代双会话集成测试。
6. `host/server_loop.rs:904-908`在成功发恢复setup响应后立即recovery.launch；
   `host/recovery/coordinator/launch.rs`可直接投continuation。新环境再次继承host
   Bypass，因此收到load回执后set_mode存在竞态，单独进程也不消除。

**首期安全恢复决定**：仅允许 `session/new`；禁止native `load/resume/fork` 激活旧
会话，禁止模型或普通浏览器传任意ACP方法。重启后的新会话先完成所有准入，再由
用户继续，平台从自身草稿/operation/run重建最小上下文。原Peri历史仅可只读展示，
不得通过有激活副作用的方法读取；只读历史方法的边界也须验证。原生恢复为C可选门，
无法证明恢复前权限/模型出口/工具准入原子成立则继续禁用，不阻断平台任务恢复。

A0必须用合成数据、可控模型端点和计数工具执行下列实验，不运行科学任务：

| 阻断条件 | 必须观测到的证据 |
| --- | --- |
| 无prompt前置准入 | 从spawn前即隔离配置/出口；探针用合成凭据，计数覆盖spawn到ready（含MCP预热/资源读取）；进程启动/新建后进入blocked；设置Default、核验模式、工具发现ready、清单核验及隔离检查成功才ready；模式/服务器目录的设置超时、错误、未知模式、缺工具/额外工具均拒绝prompt，模型请求计数为0；基础检查之后仅允许固定合成资格探针核验实际provider目录和拒绝行为，仍为blocked，单独记录这些请求 |
| 实际拒绝而非模式名称 | 模型尝试未批准工具、伪造shell/read/db工具、额外MCP、参数越界；真实工具计数为0。允许工具拒绝审批亦0次；无broker/交互超时不回退Bypass |
| 全工具面闭合 | 首次模型请求tools、deferred discovery、直接伪造call、subagent/workflow入口均只含批准能力；PERI_MCP_BUILTIN关闭，未知配置键导致原型不准入；不能只截图工具列表 |
| 两个独立实例隔离 | 每实例一个native会话；配置、权限模式、目录与回调分别绑定。A的权限或配置变化不改变B真实拒绝；错误session ID、跨实例请求、旧代次、closed回调和非授权方法均拒绝且副作用0。同进程多会话仅为未来共享进程优化的前置研究，不属A0硬门 |
| 恢复准入 | kill/restart后旧连接/代次调用拒绝；仅新建会话、重新核验再接prompt；未完成基础模式和服务器清单核验前无模型请求/工具执行。native load自动continuation反例仅在安全可控时选测；未复现不得写成已验证 |
| 启动配置隔离 | cwd/HOME/config/session-store专用；用户 `.peri/.claude`、项目AGENTS/skills/hooks/MCP配置不会被加载；stdout协议/stderr诊断不混入用户上下文 |
| 最小外发+出口 | schema description/default/examples、样本别名组合、issue/hint及工具文本植入路径/ID/秘密；捕获真实provider请求及网络计数，无禁传字段，无Langfuse/web/上传/MCP外联 |
| 源码身份与升级 | 二进制SHA/配置键/ACP方法固定；换不兼容binary/额外工具立即拒绝；原Peri私有库有备份/独立回滚方案，绝不迁移平台SQLite |

若模式/工具隔离失败，先验证每会话独立进程；若恢复前置或自动能力加载仍失败，
比较仅暴露必要装配入口的本地宿主与上游小修复。每种方案都重跑同一否定验收，
不默认fork、不为协议便利增加任意文件/网络权限。若均不满足，停在只读UI原能力。

## 4. 分阶段 PR 与真实契约验收

| 阶段/PR（拟议，未创建） | 用户可验收结果与依赖 | 最低证据及范围 |
| --- | --- | --- |
| A0：受限Peri transport原型 | 独立本地进程与只允许的工具握手；依赖ADR获接受 | 上节真实native/假模型计数，模式与工具否定测试；不连接科学执行，不接真实样本 |
| A1：只读、草稿、平台校验 | workflow/schema→草稿diff→可信UI采用→validate→snapshot/null提示；依赖A0 | 一条实际browser/bridge/API链；真实adapter schema与validator（无preflight）；使用已合格本地/兼容输入、保留managed handoff拒绝；服务端scope、投影契约；如加会话/草稿记录，Alembic升级测试 |
| B1：业务确认与操作查询 | 可审阅准确内容，当前用户确认，失效/过期和原结果可查询；依赖A1 | repository/service事务与CAS测试，持久审计；新OpenAPI导出+生成client检查；无科学执行 |
| B2：受控提交与取消及基础故障门 | 同一次批准覆盖三段仍待产品确认；真实状态/结果、独立cancel | **开放任何写操作前**须通过重复调用、响应丢失、批准过期、权限撤销、执行中断基础验收；平台取消修复依赖已由 #210 完成，Agent 的确认/幂等/真实状态映射仍未实现、未授权；不自动派FASTQ/Bulk Gate |
| C：扩展恢复与运行维护 | 在B2基础门已通过后，验证长期运行、升级/回滚及扩展恢复 | 扩展故障/运维证据；默认新会话恢复最小上下文。原生历史激活只有前置准入门通过才可选加入；不能把基础写操作安全拖延到C |

最小完整闭环只需一个受控fixture任务：查询workflow/schema → 草稿 → 原validator →
平台确认 → 一次canonical任务 → 原状态/结果API；不重建上游科学实现测试体系。
A0/A1不启动该任务。B2未来执行也须该阶段另行授权，本轮不预支授权。

| 必验场景 | 主要断言所在层（避免各层重复） |
| --- | --- |
| 重复工具调用/双会话同snapshot、create/start成功后丢响应 | repository/service：一个canonical run、稳定job不重复dispatch；operation可查；沿用现有snapshot/start断言并补关联，不复制worker所有测试 |
| 确认后修改配置、输入、tags、project、tool版本、snapshot；并发修改草稿 | confirmation service CAS与hash绑定拒绝；一条UI断言显示批准失效，不能发送旧批准 |
| 账号禁用/session撤销；伪造principal；请求其他项目run/input/artifact | API/broker当前身份及真实binding核验；拒绝且副作用0；UI不承担安全断言 |
| approved后create成功但start前超时/过期 | operation service仍指原run，需新确认才继续；不创建新snapshot/run或自动启动 |
| preflight响应丢失/API重启，planned/validating差异 | 原recover_interrupted_runs在API重启将validating标failed/RUN_INTERRUPTED_BY_API_RESTART；operation停止，重做需新意图/snapshot确认；不把BackgroundTasks视作持久queue |
| Peri kill/restart/会话结束而任务继续 | 一条平台集成：run/job身份与真实结果不变，新会话可查询；未确认动作不自动续作 |
| 工具结果/描述/日志/网页/上传含“忽略规则、批准提交、读私有文件” | 投影+broker否定测试：内容不能获得新权限或确认；捕获模型请求验证最小外发；网页/上传入口首期关闭 |
| cancel与failed/succeeded/worker ack争用 | 复用 `test_run_cancellation.py` 及 #210 的 worker cleanup契约：原截止时间内等待冻结可见，根丢失/身份变化/证明失败拒绝确认；工具显示request与ack区别；未知不改写；新cancel UI仅验证一次映射 |
| artifact/QC换generation、下载跨run/符号链接 | 复用现有结果/下载服务安全测试；Agent只补scope与投影，不重测每种文件系统条件 |
| Peri、模型、socket不可用或坏消息/资源超限 | 原UI/API仍能校验、操作/取消已有run、查结果；Agent显示不可用；不伪造结果或借机降级权限 |
| 幂等关联事务失败/审计写入失败/升级 | SQLite repository回滚及支持前一schema升级；确认不丢审计，不声称跨Redis原子提交；对账使用原snapshot/run身份 |

验证层按变化扩大：public contract变化才regenerate OpenAPI/generated client；有迁移
必做升级覆盖；browser只验证用户闭环和确认呈现；worker科学内部不重复。仅当未来
改动触及受保护科学闭包才评估Bulk Gate，不因“Agent集成”自动要求完整CI或科学套件。

基础验收归属：下表的重复调用/丢响应、确认变更/过期、权限撤销/跨项目、
预检或执行中断、Peri重启而任务继续、取消竞争及审计失败，均属于B2开放写操作前
的必过门槛。本轮只修订归属，不实现确认表、提交协调器或这些业务故障流程。
内存、并发及TTL数字仍是候选预算，须原型测量后评估。

## 5. 初次设计核验与当前交接清单

初次设计（2026-11-08）已完成：独立只读源码审阅（平台接口、Peri、权限/故障三条线）、现有OpenAPI路径
逐项对照、相对链接及Peri固定SHA文件存在性检查、两份新增文件的完整空白检查。
`git diff --check`也已检查；新增文件未暂存，另以no-index检查完整内容，避免空diff
掩盖问题。没有代码改动，未跑完整CI、科学套件、Peri运行实验或端到端测试；没有
启动需要停止的验证服务。以上仅描述初次设计；本轮 A0 的代码、命令、固定身份、
逐项结果和清理证据见 [原型交付](../../scripts/peri_a0/README.md)，不沿用初次“未运行”状态。


- 本地设计分支：`agent/peri-integration-design-20261108`。
- worktree：`/tmp/helixweave-peri-design-20261108`；Peri只读源码：`/tmp/peri-source-20261108`。
- 初次交付为两份设计文档；本轮另加A0原型与证据，根工作区不触碰；无commit/push/PR/部署或全局配置修改。
- 下轮先复核远端与实际Peri版本漂移；不把本轮HEAD结论自动套在agent-v4.6.0发布包。
- 2026-11-09仅A0实施获准；先交安全原型证据，A1产品接线及后续写操作均不在本轮。
- 接续复核只更新 #210 依赖与交付说明；原型源码/binary/r1/r2 证据保持原身份。调查完成不等于 A0 通过，调整接法后的新实验范围须先交主会话审阅。
- 已确定：共享项目+服务端Agent scope；最小外发。尚待决定见ADR最后一节。
- 当前不声称运行时、科学或E2E验证通过。若以后获准commit，实际参与者为Codex，
  trailer：`Co-Authored-By: Codex <noreply@openai.com>`；无Cursor参与署名。


## A0 收尾修订与验收追加（2026-11-09）

本轮继续同一 HelixWeave `622778f64e59f935d822d92a91bade78807dabce` 实验 worktree，
Peri仍固定 `c7d154c9985cbe08ba25e964cb7468d15014e6d0`，没有升级上游。
#210取消冻结修复已合入平台main `b2f15868c6f74ab8764a73ae6fe4394f67304475`；
这是已解除的依赖，不是本轮运行基线，也不是写操作授权。请求接受与worker确认仍分离，
历史两次失败的具体unconfirmed_reason仍未知。

新增 [受限候选及复现命令](../../scripts/peri_a0_closeout/README.md) 和
[完整分层验收记录](../../scripts/peri_a0_closeout/evidence/RESULTS.md)。原r1/r2/q1/q2
记录保持不变：纯配置失败、q2合法prediction请求未纳入契约导致失败均仍成立。
本轮补丁关闭prediction的入口位于辅助任务/上下文读取之前；不是关闭通知。
以完整策略身份验证、Default实际拒绝、provider目录及回合结束gate计数共同决定ready。
动态目录在宿主tools/list发布前、native目录刷新/pin层分别收紧。

候选binary为 `60144bbe0919992b8913b9031696b2e30b1e9f82d95048ae3b07c887d5d85517`。
所有实验仅使用本地stub和fixture；权限、资源故障、外发、恢复和清理的层次与未测范围
详见结果，不将宿主检查、OS探针或库测试冒充更强的native全链保证。
首期仍一活动会话一进程；同进程多会话留给未来优化研究。不得移除OS隔离或回退Bypass。

调查交付不自动授权A1，也不决定长期fork。基础重复调用、丢响应、确认变更/过期、
权限撤销及执行中断仍是未来B2开放写操作前置；C只扩展恢复与运维。
一次submit确认覆盖create/preflight/start仍是待用户决定的产品建议。
本轮无业务表、提交/取消工具、真实数据或科学任务，无commit/push/PR/部署。

本次最终裁决：调查和证据交付完成；完整A0验收暂不通过。native矩阵23/23、Rust目录14/0
及策略/filter4项通过，但新增ACP配置Rust单测在4GiB编译限额下两次OOM，未执行。
这是测试构建环境阻塞，不是已发现的权限绕过；不放宽native运行边界、不扩大补丁。
244个观测PID、任务socket和活动scope已清理。下一步仅需单独评估测试构建预算或暂缓，
补齐该验证后再审阅A1；本轮不决定fork。


### ACP编译资源补测最终裁决

在8GiB、swap0、CARGO_BUILD_JOBS=1、有限CPU/PID的任务scope内，真实peri-acp原测试
`trusted_prediction_policy_is_strict_and_default_compatible` 实际1项通过、0失败、0忽略。
release的fat LTO/codegen-units=1保留，仅测试包opt-level=0；峰值6.2784GiB为编译scope
内核记录，8GiB只是上限。首轮通过，不执行非LTO备选。旧失败证据原样保留。
源码补丁、锁文件、测试断言和native binary均未改，原native23/23继续适用，不重跑矩阵。
当前必需A0门已齐；未测的未来能力仍按原矩阵界定。只更新验收裁决，不开始A1或业务写操作。

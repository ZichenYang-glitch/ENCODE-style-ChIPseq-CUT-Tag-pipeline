# Agent 运行操作边界

本文件固定 create、preflight、start 和 cancel 相对现有运行生命周期的边界。它复用
[平台总览](platform-overview.md) 的只读 Agent 约束、[Peri 契约附录](peri-agent-integration-contracts.md)
里已经对照过的 HTTP 入口，以及 [目录 / Agent / Viewer 边界](catalog-agent-viewer-boundaries.md)
对运行操作另立契约的约定。不新增任务 API，不注册写工具，不开放受管输入执行，也不把
一份设计当成运行已经成功。

对照的是当前服务与路由，不是拟议中的 `submit_confirmed`。附录里“一次批准覆盖
create、preflight、start”仍是当时未决的建议；这里不采用它。启动必须是人对具体 run
的另一次确认。

## 现有入口

这些路由都在 `/api/v1` 下，并经过成员会话和 CSRF。create、start、cancel 另外把当前
principal 写入运行安全审计。它们不要求管理员，也没有按项目再收一层 ACL。模型文本、
工具结果和“用户已同意”都不能代替这次会话。

| 步骤 | 操作 | 服务 | 当前结果 |
| --- | --- | --- | --- |
| snapshot | `POST /workflows/{workflow_id}/validate` | `ValidatedInputService` | 执行可用且校验通过时保存不可变 snapshot。执行不可用时成功体里的 snapshot 为 null，没有可提交身份。 |
| 创建 | `POST /workflows/{workflow_id}/runs`，`createRun` | `ValidatedRunCreationService.create_run` | 请求体只有 `snapshot_id` 和 `tags`。首次消耗为 201、`created`；同一 snapshot 配同一 tags 再提交为 200 和原来的 run。 |
| 预检 | `POST /runs/{run_id}/preflight`，`triggerPreflight` | `RunService.transition_run`，随后 `LocalPreflightService` | 只接受 `created`。成功先变成 `validating` 并返回 202，干跑在 HTTP 响应之后的后台任务里做完才到 `planned`。离开 `created` 后再调用是 409 `PREFLIGHT_ALREADY_TRIGGERED`。 |
| 启动 | `POST /runs/{run_id}/start`，`startRun` | `RunSubmissionService.start_run` | 显式提交。`created` 和 `validating` 为 409。`planned` 才预约执行身份并入队；`queued` 或已带 dispatch 标记的后续状态返回该 run，不再入队。 |
| 取消 | `POST /runs/{run_id}/cancel`，`cancelRun` | `RunCancellationService.cancel_run` | 只作用于路径里的这个 `run_id`。尚未 running 时直接进入 `cancelled` 并返回 200。running 时先记下停止意图，202 只表示已请求停止。已终态的 run 原样返回，不再发新的停止。 |
| 读回 | `GET /runs/{run_id}`，以及 `/events` | `RunService.get_run` | SQLite 里的 status 是终态来源。公共 `RunRecordResponse` 含 status、时间和 `cancellation_reason`，不含原始 inputs、snapshot id，也不含取消请求时间或 worker 确认时间。 |

Agent 聊天仍走 `POST /workflows/{workflow_id}/agent/chat`。`ReadOnlyToolRegistry`
只登记 `list_workflows`、`get_workflow_schema`、`validate_workflow_input` 和
`explain_issues`，并拒绝 `run_`、`submit_`、`execute_`、`delete_`、`update_`、
`apply_`、`write_`、`modify_` 开头的工具。`validate_workflow_input` 调用的是
`ValidationService`，不写 snapshot，也不能拿去创建 run。

工作台现有按钮与上述入口对齐。人点击 “Create run from validated inputs” 时，请求只带
这份 snapshot，tags 用接口缺省的空映射。成功后页面打开该 run，并带上一次 preflight。
若这条 run 已经离开 `created`，服务器回 409，不会再排一次预检，也不会因此启动。
“Start run” 和 “Cancel run” 都在该 run 的进度页上。

## 人确认什么

确认发生在平台页面，使用当前成员会话。一次确认只覆盖下面的一件事。

| 人的动作 | 确认的对象 | 不包含 |
| --- | --- | --- |
| 用当前草稿校验并看到非空 snapshot | 服务器刚保存的那份 snapshot | 不创建 run，更不启动 |
| Create run from validated inputs | 该 `snapshot_id` 加上缺省空 tags；成功后页面为这条 run 触发一次 preflight | 不启动执行；run 已离开 `created` 时预检得到 409，不另排一次 |
| Start run | 这一条已经 `planned` 的 `run_id`，以及预检时冻结的构建身份 | 不换 snapshot，不另开一条 run |
| Cancel run | 路径中的这一条 `run_id` | 不取消同 snapshot 的概念本身，也不把 202 当成已经停止 |

snapshot 过期、构建身份变了、tags 与已消耗的 run 不一致，或执行资格不再可用时，
创建或启动拒绝，人需要重新校验后再确认。旧确认不能套到新 snapshot 上。

## 幂等键绑定什么

这些路由没有 `Idempotency-Key` 头。键就是已经存在的对象身份，不能另造一把客户端
键去开第二条 run 或第二个 job。

| 操作 | 键绑定 | 断线后怎么重试 |
| --- | --- | --- |
| 创建 | `snapshot_id`，加上同一组 tags，以及该 snapshot 已冻结的项目/样本、输入、参考和实验绑定 | 同一 snapshot、同一 tags：200，返回原来的 run。tags 或绑定对不上：409 `VALIDATED_SNAPSHOT_REPLAY_CONFLICT`，仍不新建。没有按 snapshot 查询 run 的 GET，丢失的创建响应靠这次重放找回来。 |
| 预检 | 这个 `run_id` 且 status 仍是 `created` | 第一次成功后 status 已离开 `created`，再调得到 409，不会再排一次后台预检。响应丢失时先 GET：仍是 `created` 才可以再触发一次。 |
| 启动 | `run_id` 决定 job。`build_execution_job_id` 是 `run-execution-` 加该 run id 的 SHA-256；同一 run 只保留第一条执行分配 | `queued`，或 running 及终态上已经有 dispatch 标记时，再次 start 返回这条 run，不第二次入队。队列没有确认时 run 保持可重试，不把失败写成 `queued`。 |
| 取消 | 这个 `run_id`；running 时还要带上已领取分配的 job、backend 和 queue | 终态上再取消返回当前记录。running 上重复取消不另写一条取消意图，也不新派一个执行 job；它可能对同一分配再请求一次停止。是否已停止仍看 GET 的 status。 |

重试使用同一成员会话。会话已经无效时请求在写入前被拒绝，不能为了补发而换成模型侧的另一套身份。

## 怎样读真实终态

终态只有 GET 返回的三种 status：`succeeded`、`failed`、`cancelled`。`created`、
`validating`、`planned`、`queued`、`running` 都不是终态。

- 首次创建的 201 表示新 run，status 为 `created`。同一 snapshot 的 200 返回这条 run 的当前记录；status 可能已经离开 `created`。200 不表示重新创建，也不表示重新入队。
- 预检的 202 只说明已经进入 `validating`。`planned` 要等后续 GET。
- 启动的 202 只说明这次调用被接受。入队与否看响应体或下一次 GET 的 status。
- 取消在 running 上的 202 只说明停止已请求。status 仍可能是 `running`。公共响应没有请求时间和 worker 确认时间，不能从 202 或单独的 `cancellation_reason` 推断已经取消。
- 预检期间 API 重启时，`recover_interrupted_runs` 把仍为 `validating` 的 run 标成 `failed`，issue 为 `RUN_INTERRUPTED_BY_API_RESTART`。这就是终态。该 run 不能再预检；恢复提示里的 “submit a new preflight” 与当前路由不符，因为路由只接受 `created`。要继续，必须重新校验并确认一条新 run。
- `queued` 或 `running` 但没有对应的 dispatch 或 claim 时，同一次恢复可以把它们标成 `failed`（`RUN_ORPHANED_AFTER_API_RESTART`）。已经有确认标记的派发不靠这次恢复改写。展示以恢复之后的 GET 为准。

## 例：一次获准的创建和启动

成员会话有效，目标工作流的执行资格为 available，校验返回非空 snapshot
`vsnap_0123456789abcdef0123456789abcdef`。人看到这个 id，用空 tags 点击创建。

1. `createRun` 首次返回 201。响应给出 `run_id`，status 为 `created`。人确认的是这份 snapshot 和这组 tags。
2. 进度页对这个 `run_id` 调用 `triggerPreflight`，返回 202，status 为 `validating`。GET 后来变成 `planned` 才表示预检完成。此时还没有执行 job。
3. 人在这一条 run 上点击 Start run。确认的是这个 `run_id`。job id 由该 run id 算出，形如 `run-execution-` 加它的 SHA-256。调用返回 202 后，用 GET 读到 `queued` 或更后的真实 status。
4. 终态再次 GET。只有 `succeeded`、`failed` 或 `cancelled` 可以当作这次运行的结果。

创建响应若丢失，用同一个 snapshot id 和同一组空 tags 再调用 `createRun`。结果是 200 和同一个 `run_id`，没有第二条 run，也没有因为这次重放就入队。

## 例：撤权，以及断线后的启动重试

上例的 run 已经是 `planned`，人点击了启动，但浏览器没有收到响应。

- 会话仍有效，而且 GET 已经是 `queued` 或更后并带有这条 run 的 dispatch：再调用 `startRun` 得到同一条记录，job id 不变，不第二次入队。页面继续显示这次 GET 的 status，而不是“启动失败所以再派一次”。
- 若 GET 仍是 `planned`，且上一次是队列未能确认的 503：可以用同一次启动确认重试。重试仍绑定同一个 job id。在 GET 变成 `queued` 或更后之前，不把它展示成已经启动。
- 重试前会话被撤销，或账号被禁用：`resolve_session` 不再给出 principal，下一次 start、preflight 或 cancel 为 401，不写新状态。已经创建的 run 保留。撤权不会补做启动，也不会收回一个已经派发的 job。之后要取消，仍须一个有效成员会话，并且路径里仍是这个 `run_id`。在 GET 变成 `cancelled` 之前，展示保持 `queued` 或 `running`。
- 启动前执行资格撤回，或构建身份与预检时不一致：start 返回 503 或 409，不入队。GET 仍是 `planned` 或其他拒绝时读到的 status。

若断线发生在预检的 `validating`，并且 API 已重启，GET 的 status 是 `failed`，issue 为 `RUN_INTERRUPTED_BY_API_RESTART`。不再对这条 run 预检或启动。人若要再跑，先重新校验，再确认新的 snapshot。

## 当前不能称为可运行的路径

- 工作流执行资格不是 available：校验可以成功，但 snapshot 为 null。没有 snapshot 就不能创建，也不能把这次校验说成可运行。
- 选择了受管修订而执行资格不可用：校验失败，码为 `INPUT_USE_CAPABILITY_UNAVAILABLE`，没有 snapshot。
- 输入用途的 provenance 为 `managed_revision_v1`：`ManagedInputVerificationService.verify_run` 在执行前返回 `MANAGED_INPUT_EXECUTION_UNAVAILABLE`，不进入实际工作流。即便 snapshot 或 run 已经写上，这条路径仍未开放。本轮不改这个交接。
- Agent 的 `validate_workflow_input`、聊天里的说明，以及附录里尚未实现的 `submit_confirmed` / `cancel_confirmed`：都不是运行入口。
- 过期、陈旧、缺失或绑定损坏的 snapshot：创建返回 409 或 404，不产生新 run。
- `validating` 和 `planned` 不是已启动。取消的 202 不是已取消。

只读的实验反查不在本文件开放执行。成员读一条已有 run 的终态，仍用上面的 GET。

## 本轮不做

不增加 operation 表、写工具、确认协调器、OpenAPI 操作或迁移。不修改预检、启动、取消或受管输入的服务行为。不派发真实或合成科学任务。后续若要实现 Agent 侧的获准提交，必须先有单独授权，并沿用本文件的四步顺序、现有幂等身份和 GET 终态；不能用一次工具调用走完创建、预检和启动。

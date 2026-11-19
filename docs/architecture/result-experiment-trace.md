# 从已发布结果反查冻结的实验修订

M0 契约不把 `experiment_binding` 放进运行记录或发布摘要。本读取沿用该边界：发布详情和 `getRun` 保持原样，只增加从一条已发布结果回到运行消费时复制的修订。

## 已有读取能回答的问题

- `getArtifactPublication` 回答这一代发布的身份，以及运行级样本范围。它不指向实验修订，也不把某个样本解释成单条轨道的贡献。
- `catalogListExperimentRevisionRunResults` 从实验修订正向列出发布结果。无绑定的修订返回空结果和 `experiment_snapshot_binding_absent`。
- `getRun` 回答运行状态、阶段，以及与绑定一致的参考身份。参考证据不一致时返回 `RUN_REFERENCE_EVIDENCE_INVALID`。
- `listRunQcMetrics` 回答已保存的 QC。没有结果代状态时，原接口报错。
- `listRunEvents` 回答已记录的事件，其中包括 `preflight_failed`、`dry_run_failed` 和带 `PREFLIGHT_FAILED` 的状态事件。触发预检仍是单独的 POST。

## 补上的反向查询

`GET /api/v1/catalog/projects/{project_id}/artifact-publications/{run_id}/{artifact_id}/experiment-trace?generation=...`

操作 ID 是 `catalogGetArtifactPublicationExperimentTrace`。它只读取 `artifact_publications`、`run_experiment_bindings` 和 `experiment_revisions`，不新增表，也不回填旧运行。

可见范围与目录的其他项目读取相同。发布不存在、项目不对、Legacy 项目、或代数不对时返回 `catalog_not_found`。信封始终带上这一条已发布结果，关联类型保持 `run_input_scope`。

- 运行没有实验绑定：`trace_status=unresolved`，`reason=experiment_snapshot_binding_absent`。不返回修订 ID 或摘要。
- 复制的修订 ID 能找到同一项目的修订，且修订摘要与运行上复制的摘要一致：`trace_status=experiment_snapshot`，返回那一对复制值。实验之后出现新修订，不改变旧结果的答案。已归档实验仍可按修订 ID 读取。
- 绑定在，但修订缺失、项目不一致，或摘要不一致：`trace_status=unresolved`，`reason=experiment_revision_evidence_inconsistent`。不回显不可信的 ID 或摘要。

正向列表的未解析信封仍然是空结果。反向查询不同，因为成员打开的就是这一条发布。

## 结果页与 Agent 的同一次只读

`GET /api/v1/agent/projects/{project_id}/artifact-publications/{run_id}/{artifact_id}/context?generation=...`

操作 ID 是 `readPublishedResultContext`。结果页和 Agent 使用调用者当前的成员会话。能打开发布详情的结果才能读到上下文；请求的项目与这条发布不一致时仍是 `catalog_not_found`。Legacy 发布沿用详情的可见范围，返回 `trace_status=unresolved` 和 `experiment_snapshot_binding_absent`，不把 Legacy 项目放进目录追溯。目录接口对 Legacy 项目和错误项目保持原来的 404。它不调用触发预检、创建运行或启动运行。

缺值用证据字段表示，不补猜：

- `reference_evidence` 为 `absent` 或 `present`。参考绑定缺失时是 `absent`。身份不一致时返回 `RUN_REFERENCE_EVIDENCE_INVALID`，与 `getRun` 相同。
- QC 只在结果状态把 `qc_artifact_generation` 绑到这一代产物时才可读。请求的不是当前产物代时为 `unverifiable`，不返回其他代的指标或 QC 代。当前代没有指标行、也没有 QC 摘要时为 `absent`。当前代仍有指标行，但 QC 代号、清单摘要或关联产物代缺失时，与 `listRunQcMetrics` 一样返回 `RUN_QC_METRIC_DATA_INVALID`，不把这种不一致说成从未记录。摘要绑到当前代时为 `present`，并沿用 QC 列表的分页。
- `preflight_failure_evidence` 只反映已经记录、并且挂在这条已发布结果上的失败。事件上下文不外传；带路径的消息改成固定说明。事件页读不完时返回 `RUN_EVENT_READ_INCOMPLETE`，而不是声称没有失败。
- 预检在任何发布之前失败时，没有产物代可供这个入口读取。Agent 只凭 run ID 读取 `getRun` 和 `listRunEvents`。运行状态说明失败，事件的 `context.reason_code` 和 `issue.context.reason_code` 保留具体原因，例如 `REFERENCE_PROFILE_BINDING_INVALID`。这两次 GET 不改变运行状态，也不增加事件。

项目权限和历史绑定数据保持不变。历史列表的显示修复留在 IRENE-42。

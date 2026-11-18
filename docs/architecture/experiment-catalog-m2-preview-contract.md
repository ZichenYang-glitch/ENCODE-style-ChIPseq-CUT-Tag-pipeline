# 实验数据目录 M2：离线建议预览、子集选择与修订冲突契约

状态：IRENE-27 设计契约与可执行验收基线（2026-11-19，已按复核意见修订）。本文固定 M2
离线建议预览的用户动作、角色边界、合成目录夹具与可复现验收用例，配套可运行的内存契约夹具与测试；
**不新增** Alembic 迁移、正式批量写入 HTTP 路由、OpenAPI/生成客户端改动或目录网页改动。

上游依据：

- [实验数据目录 M0 契约](experiment-catalog-m0-contract.md)（M1 Slice A 为唯一读写契约基线；分类任务不依赖 M1 Slice B 的 `snapshot_experiment_bindings` / `run_experiment_bindings`）
- [数据目录、受控 AI 整理与可视化边界设计](catalog-agent-viewer-boundaries.md) 第 4 节（`M2：Peri 数据整理`）
- [下一批交付计划](../development/next-delivery-session-plan.md)（`M2` 与 `IRENE-27` 验收范围）

## 1. 边界与职责拆分

M2 将**提案生成（Proposal）**与**用户确认采用（Adoption）**严格分开：

| 层次 | 拥有方与职责 | 明确不做 |
| --- | --- | --- |
| 离线提案 / M2 预览契约（本契约与可运行 fixture） | 固定提案结构、建议引用（`citations`）、角色投影、管理员子集选择、不可变 `preview_digest`、精确基线与来源核验、一次性采用状态、M1 A 范围幂等与回执形状 | 不写目录数据库迁移；不实现正式批量写入 HTTP 接口；不修改目录网页；不依赖 M1 B；不调用真实模型/外部 provider |
| 正式迁移、公开 HTTP 路由与生成客户端 | 仍由目录负责人（`HW-Catalog`）在后续采用工单（如 IRENE-32）统一整合，避免并发改写 Alembic 头和 `frontend/src/api/generated/` | 不在本契约中提前占用 Alembic revision ID 或修改 `frontend/openapi.json` |
| 在线 Peri 整理与模型出口 | 仍受 P2 与边界设计第 4 节约束 | 不从模型直接提交修订；不把未脱敏路径、原始日志、序列、供体标识或自由文本外发 |

默认边界数值沿用[边界设计](catalog-agent-viewer-boundaries.md)第 4 节：

- 单次变更集最多 `100` 条对象修改，单次预览/采用请求体不超过 `256 KiB`（`262144` 字节）。
- 提案采用窗口为 `30` 分钟（`expires_at`）；未采用提案过期后保留 `7` 天只读回执窗口（`retention_until`）。
- 已采用回执和对应修订作为平台审计证据保留，不随提案 TTL 清除。
- 写请求必须携带符合 M1 A 规则的 `Idempotency-Key`（长度 `1..128`，正则 `^[A-Za-z0-9][A-Za-z0-9_.-]{0,127}$`）。

## 2. 用户行为与状态流转

### 2.1 参与者与角色投影

- **普通成员（`member`）**：
  - 在活动项目内可以发起或读取有界离线建议，**整条建议按成员可见范围投影**：
    - `ambiguous` 槽位省略 `candidates`（仅保留 `{"status": "ambiguous", "values": []}`），不含词条 `aliases`。
    - `citations.source_refs` 同样按角色过滤，**省略任何指向歧义候选词条修订的引用及可识别候选的 `evidence_token`**，成员只看到目标对象修订及已确认（`known`）词条的引用。
  - **不能**将提案推进到预览采用或采用写入；成员调用预览采用或采用入口固定返回 `403 ADMINISTRATOR_REQUIRED`，零写入。
- **管理员（`administrator`）**：
  - 在自己的会话中查看待处理或已采用提案的逐项字段前后差异、完整**建议引用**（含 `ambiguous` 候选来源与 token，供管理员核对候选依据）以及**当前服务端头修订**。
  - 在不修改提案原文的前提下**勾选任意非空子集**（`selected_item_ids`），先生成不可变的子集预览摘要（`preview_digest`），确认被排除项不会暗中带入。
  - 携带 `preview_digest`、`selected_item_ids`、每个目标对象的 `expected_revision_number` 及 `Idempotency-Key` 提交采用。
  - 若请求响应丢失或重新登录，可使用同一 `Idempotency-Key` 与相同请求体安全重放并取回原回执，或重开已采用提案 / 按 `proposal_id`、`project_id` 重新查询原采用结果。

### 2.2 提案生命周期、一次性确认与精确校验

```mermaid
flowchart LR
    Catalog[M1 A 分类表头修订] --> Proposal[离线建议 Proposal pending]
    Proposal --> AdminSelect[管理员勾选子集 selected_item_ids]
    AdminSelect --> Preview[生成不可变预览摘要 preview_digest]
    Preview --> Check{同事务校验: 角色/幂等键/采用状态/期限/规则词表/基线ID+编号/引用来源}
    Check -->|全部通过| Commit[原子追加修订 + 更新头指针 + 标记提案 adopted + 写入采用回执]
    Check -->|提案已采用或规则/词表/引用失效| Invalidated[409 catalog_proposal_invalidated 整批回滚]
    Check -->|基线修订ID或编号与现头不等| Conflict[409 catalog_revision_conflict 整批回滚]
    Check -->|同操作者同键同操作同目标同体| Replay[返回原回执 replayed=true]
    Check -->|同操作者同键但操作/目标/体任一不同| IdemConflict[409 catalog_idempotency_conflict]
```

1. **一次确认产生一份采用结果**：提案状态为 `pending` / `adopted` / `expired`。采用提交在同一事务中将选中子集修订、头指针、`proposal_state = adopted` 与 `adoption_receipt` 一起提交。提案进入 `adopted` 后：
   - 重新打开提案会看到 `proposal_state = "adopted"`、`adoptable_by_caller = false` 以及原 `adoption_receipt`。
   - 任何管理员尝试对该提案生成新预览，或拿着此前预览的**不同子集**（即使目标与首次采用完全不相交）提交旧确认，均拒绝写入（`409 catalog_proposal_invalidated`），项目中不会出现第二份采用回执。
   - 原操作者使用**原 `Idempotency-Key` 与原请求体**重试仍命中幂等重放，返回原回执并标明 `replayed = true`。
2. **核验服务端提案基线与引用证据**：预览与采用均以服务端持久化提案及当前目录状态为权威，请求方传入的 `expected_revisions` 或 `preview_digest` 不能替代校验：
   - 校验目标对象属于 `proposal.project_id` 且项目处于活动状态（非归档、非 Legacy）。
   - 校验目标当前头修订的**修订 ID 与修订号同时匹配**服务端提案记录的 `base_revision_id` 和 `base_revision_number`（“编号相同但引用另一修订 ID”必须拒绝并返回 `409 catalog_revision_conflict`）。
   - 校验 `citations.source_refs` 引用的实验修订、样本修订、词条修订或投影行在当前目录中真实存在、归属同一项目且未归档（引用失效或不存在时拒绝写入，返回 `409 catalog_proposal_invalidated` 或 `409 catalog_term_archived`）。
   - 校验提案及每条引用的 `rule_version`、`vocabulary_version` 与 `projection_version` 仍符合当前目录生效版本；规则或词表换版后旧提案拒绝写入（`409 catalog_proposal_invalidated`）。
3. **沿用 M1 A 的幂等键范围**：幂等回执主键为 `(actor_user_id, idempotency_key)`（不把 `proposal_id` 放进查找键）。命中已有键时再比较 `operation`、`target_key`（本场景为 `proposal_id`）与 `request_sha256`：完全一致则重放原回执；同一管理员把同一 `Idempotency-Key` 复用到另一份提案、另一操作或不同请求体时，返回 `409 catalog_idempotency_conflict`，第二份提案零写入。

## 3. 与 M1 Slice A 契约的对齐

本契约复用 M1 A（`438a943`）的以下规范，不另造第二套身份或槽位格式：

- **身份前缀**：Project `prj_`、Experiment `exp_`、ExperimentRevision `expr_`、Sample `smp_`、SampleRevision `smpr_`、Term `term_`、TermRevision `trmv_`（32 位小写十六进制）。本契约新增提案身份 `cprop_`、预览回执身份 `cadopt_`（同样为前缀加 32 位小写十六进制，供未来目录负责人落地表与路由时使用或调整）。
- **栏位语义**：
  - 实验分类字段：`assay`、`target`、`batch`（单值槽位），以及每个样本钉上的 `replicate_role`（单值）与可选 `replicate_group`。
  - 样本材料字段：`species`、`tissue`（单值）；`cell_type`、`cell_line`、`treatment`、`tag`（多值 `{status, values}`）。
  - `ambiguous` 栏位及其对应候选 `source_refs` 仅对管理员可见；成员视图省略 `candidates`、候选引用及 `aliases`。
  - 非 `catalog` 自由属性（如既有 `attributes.legacy_note`）在样本分类修订时原样复制到新 `SampleRevision`，绝不丢失。
- **复用错误码**：`catalog_revision_conflict`（409）、`catalog_idempotency_conflict`（409）、`catalog_payload_invalid`（422）、`catalog_not_found`（404）、`catalog_term_archived`（409）、`catalog_experiment_archived`（409）、`catalog_term_type_mismatch`（409）、`ADMINISTRATOR_REQUIRED`（403）、`API_REQUEST_INVALID`（400）。
- **M2 新增错误码（供后续采用接口整合）**：`catalog_proposal_expired`（410）、`catalog_proposal_invalidated`（409）、`catalog_preview_digest_mismatch`（409）。

## 4. 合成分类表与管理员预览用例

静态 JSON 样例见 [M2 离线预览契约样例](experiment-catalog-m2-preview-examples.json)，可运行 Python 契约与夹具实现在 `src/encode_pipeline/platform/catalog_curation_preview.py`，验收测试在 `test/platform/test_catalog_curation_preview.py`。

### 4.1 合成目录夹具（Project `prj_11111111111111111111111111111111`）

初始均处于 `revision_number = 1`：

| 目标对象 | 初始头修订 | 当前关键分类状态 | 离线提案项 | 建议变更（仅允许字段） | 建议引用（`citations`） |
| --- | --- | --- | --- | --- | --- |
| `exp_44444444444444444444444444444444`（`Treg H3K27me3 CUT&Tag`） | `expr_4444...1`（rev 1） | `assay=known(CUT&Tag)`、`target=known(H3K27me3)`、`batch=unknown` | `item_1` | `batch: unknown -> known(batch-2026-11)` | `rule_version=m2-catalog-rules/1`、`reason_codes=["BATCH_TERM_MATCH"]`、钉住 `trmv_5656...` 与 `expr_4444...1` |
| `exp_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa`（`Treg H3K4me3 CUT&Tag`） | `expr_aaaa...1`（rev 1） | `assay=known(CUT&Tag)`、`target=unknown`、`batch=known(batch-2026-11)` | `item_2` | `target: unknown -> known(H3K4me3)` | `rule_version=m2-catalog-rules/1`、`reason_codes=["TARGET_TERM_MATCH"]`、钉住 `trmv_4545...` 与 `expr_aaaa...1` |
| `smp_22222222222222222222222222222222`（`Treg sample 1`） | `smpr_3333...1`（rev 1，含非目录属性 `legacy_note`） | `cell_type=unknown`、`tag=unknown` | `item_3` | `cell_type: unknown -> ambiguous([regulatory T cell, keratinocyte])`；`tag: unknown -> known([m2-curated])` | 公共引用钉住 `smpr_3333...1` 与 `trmv_1010...`；仅管理员可见的候选引用钉住 `trmv_bbbb...`（`alias-candidate-treg`）与 `trmv_dddd...`（`alias-candidate-kc`） |

### 4.2 验收用例（`M2-P1` 至 `M2-P7`）

| 用例编号 | 操作场景与并发/校验条件 | 预期可观察结果 |
| --- | --- | --- |
| `M2-P1` | 查看离线提案与建议引用；成员尝试预览采用或提交采用 | 管理员预览可见完整 `citations`（含歧义候选 `term_revision_id` 与 token）及 `candidates`；视图不含宿主路径、供体标识、原始日志或自由文本。**成员视图在 `item_3` 上同时隐藏 `candidates`、`aliases`、候选词条引用及其识别 token**（仅见目标 `smpr_3333...1` 与已确认 `tag` 引用）；成员提交预览或采用返回 `403 ADMINISTRATOR_REQUIRED` 且零写入 |
| `M2-P2` | 管理员 A 勾选子集 `["item_1", "item_3"]`（排除 `item_2`），生成 `preview_digest` 并提交采用 | 预览响应明确列出 `selected_item_ids=["item_1", "item_3"]`、`excluded_item_ids=["item_2"]` 与不可变 `preview_digest`。若提交篡改后的子集但沿用原摘要，返回 `409 catalog_preview_digest_mismatch`。合法采用后仅 `exp_4444...` 与 `smp_2222...` 升至 `revision_number=2`（`legacy_note` 原样保留），`exp_aaaa...` 仍停留在 `revision_number=1` |
| `M2-P3` | **并发冲突（手工修订领先）**：提案与管理员 A 预览均期望 `exp_4444...` 为 `1`、`smp_2222...` 为 `1`；采用提交前，管理员 B 以 `expected_revision_number=1` 先行修订 `exp_4444...` 至 `revision_number=2` | 管理员 A 随后提交采用时，整批返回 `409 catalog_revision_conflict` 与有界冲突明细（`target_id=exp_4444...`、`expected_revision_number=1`、`current_revision_number=2`）。`exp_4444...` 保持管理员 B 的 rev 2（不合并字段），`smp_2222...` 保持 rev 1（**不部分提交**） |
| `M2-P4` | **一次性确认与并发采用**：一名管理员预览第 1、3 项，另一管理员预览不相交的第 2 项（以及另一份重叠目标提案）；第一人先行提交成功 | 第一人成功后重开该提案得到 `proposal_state="adopted"`、`adoptable_by_caller=false` 与原 `adoption_receipt`；第二人拿着第 2 项的旧预览提交或重预览均被拒绝（`409 catalog_proposal_invalidated`），`exp_aaaa...` 停留在 rev 1 且项目内仅存在一份回执；另一提案若再针对已升至 rev 2 的 `exp_4444...` 提交预期 rev 1，返回 `409 catalog_revision_conflict` 且零写入；第一人原键原体仍可重放 `replayed=true` |
| `M2-P5` | **丢响应后幂等与 M1 A 幂等键范围**：管理员 A 提交 `Idempotency-Key=idem-admin-a-adopt-1` 成功（rev 1 -> 2），客户端丢失响应后重放；同一管理员把同一键复用到另一提案 `cprop_2222...` 或不同体 | 同操作者、同键、同操作、同目标、同体重试直接返回原回执并标明 `replayed=true`（即使提案已是 `adopted` 且头修订已为 rev 2）；按 `proposal_id`/`project_id` 可重查该回执；同键异体或同管理员把同一键复用到第二份提案均返回 `409 catalog_idempotency_conflict`，第二份提案零写入 |
| `M2-P6` | **精确基线与来源/版本失效检查**：(a) 提案 `item_1` 的 `base_revision_number=1` 与现头相等，但 `base_revision_id` 指向另一修订 ID；(b) 提案 `citations.source_refs` 指向不存在/跨项目/已归档修订；(c) 当前生效的 `rule_version` 或 `vocabulary_version` 已升级 | 即使请求方传入匹配编号与有效 `preview_digest`，(a) 仍整批返回 `409 catalog_revision_conflict`（注明 `expected_revision_id != current_revision_id`）；(b)/(c) 整批返回 `409 catalog_proposal_invalidated` 或 `409 catalog_term_archived`；同批第 2 项同样**零写入** |
| `M2-P7` | 提案超过 30 分钟采用窗口；管理员在提交前被降权为 `member`；请求携带宿主路径或超过 100 项限制 | 过期提案拒绝采用并返回 `410 catalog_proposal_expired`（7 天内仍可只读查看）；被降权管理员提交返回 `403 ADMINISTRATOR_REQUIRED` 且零写入；宿主路径及超限变更集返回 `422 catalog_payload_invalid` |

## 5. 未纳入本阶段的项目（待后续集成）

- Alembic 数据库迁移、持久化 SQLAlchemy 表（`catalog_proposals`、`catalog_adoption_receipts` 等）及并发事务实现（由目录负责人统筹，接 IRENE-32）。
- 公开 FastAPI 路由与 `frontend/src/api/generated/` 客户端重新生成。
- 目录网页（IRENE-23）与目录 CLI 的实际消费接线。
- 在线 Peri 工具调用与真实外部模型出口（P2）。

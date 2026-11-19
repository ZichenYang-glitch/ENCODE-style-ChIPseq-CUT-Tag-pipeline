# 目录绑定入口

普通校验页提交已有 `validateWorkflow` 的可选 `experiment_revision_id`。Peri 草稿只提供跳转，把已保存的配置、项目和样本修订带到这一页。跳转、预选实验和旧的校验结果都不是新的实验绑定。

## 日常入口

会留下可执行 snapshot 的入口是普通校验页：`validateWorkflow`，然后 `createRun`。配置草稿页写明校验不启动运行。`authoringValidateDraft` 的请求仍只有 `expected_revision`。M1 契约也不让这条草稿校验接收实验修订，也不在 `createRun` 消费 snapshot 时补实验。

因此这里不扩展草稿校验 API，也不改共享目录 API、迁移或生成客户端。

## 从已保存草稿带走的内容

只读已保存草稿，不读未保存的编辑。

| 字段 | 处理 |
| --- | --- |
| `config`、`options` | JSON 对象才复制到普通校验页 |
| `project_id` | 活动项目。Legacy Project 不复制 |
| `sample_revision_ids` | 非空、无重复、形状合法才复制，并保留原顺序 |
| `samples` | 字符串单元格组成的表才复制。路径字符串不复制 |
| `reference_profile_revision_id` | 形状合法才复制 |
| `input_selections` | 形状合法才复制 |

`workflow_samples`、草稿 id、草稿修订、旧 snapshot 和草稿校验结果都不进入 `validateWorkflow`。样本计划留在原草稿。

## 实验修订

普通校验页用活动项目的 `catalogListExperiments` 列出未归档实验修订。默认不选。用户可以改提交的样本修订顺序，或把其中一个改成另一修订。

只有请求里的实验修订属于该活动项目，且与提交的样本修订集合都非空、无重复并精确相等时，服务端才冻结。顺序可以不同。不相等时整次校验失败，不写绑定。

界面只根据这次响应的 `experiment_binding.status` 显示：

- `unresolved`：未绑定
- `experiment_snapshot`：这次 snapshot 冻结了返回的实验修订

省略 `experiment_revision_id` 的请求保持可用，响应为未绑定。

## 不改变的请求

- `authoringValidateDraft` 仍只提交 `expected_revision`
- `createRun` 仍只提交 `snapshot_id`
- 原草稿和历史运行不补绑定

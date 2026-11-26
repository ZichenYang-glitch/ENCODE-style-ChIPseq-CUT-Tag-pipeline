# 工作流字段投影（authoring-field-projection/1）

本版把建议字段从公共服务里的 ENCODE 专用函数，改成由工作流 adapter 声明。
摘要算法保持 `hw-a1/1`。公共 HTTP 模型、登录授权、草稿修订、采用和校验入口都不改。
Peri 传输和目录 API 不在本版范围。

路线依据是 [#221](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/221)
的 `3554f5f`：字段投影交给 adapter，草稿保存、版本检查和提交仍由公共服务负责。
样本计划不在字段投影这一版的接口里。后续契约见下文「样本计划」。

## 接口

adapter 用同一个函数提供建议目录：

```text
project_authoring_fields(schema_payload) -> AuthoringFieldProjection | None
```

函数可以是 adapter 类上的方法。类文件若已在已验证的执行闭包里，则不要改那个文件，
改为由同一包的 `authoring` 模块导出同名函数。平台先调用类方法；没有该方法时，
再加载该类所在包的 `authoring.project_authoring_fields`。同一个包里的子类，
包括定义在 `hitrac_preprocess.results` 的结果 adapter，使用这份导出。

`schema_payload` 是平台刚刚纳入 `schema_digest` 的那份 schema 文档。
返回值只含建议可以读写的标量字段，以及这些建议可以引用的问题标签。
字段约束是 adapter 策略和这份 schema 的整数交集。默认值、描述和示例不进入投影。
返回 `None`，或函数抛错，表示这个工作流没有建议目录。平台因此给出
`CAPABILITY_UNAVAILABLE`，手工草稿和 adapter 校验仍可用。

`AuthoringFieldProjection` 含有：

| 部分 | 归属 |
| --- | --- |
| `fields` | adapter：token、标签、所在文档、键、整数上下界 |
| `issues` | adapter：建议可以解释的问题代码和固定句子 |
| 授权、草稿修订、采用、校验 | 平台 `AuthoringService` |
| `schema_digest`、`projection_digest` | 平台，算法仍是 `hw-a1/1` |

没有这个方法的 adapter，例如 `bulk-rnaseq`，保持没有建议目录。

## ENCODE 旧草稿

ENCODE 仍只投影 `config.threads` 和 `options.cores`，问题标签的内容和顺序不变。
因此已保存草稿上的投影摘要不变：

`b3ea5bb01c2457fc788508187d967445c1eca0f1c95f26d0a6560b6a497ac71b`

打开连接时，平台仍用这份摘要和当前 schema 摘要比对。摘要一致就可以打开、保存和校验。
摘要若被改写，已有连接会按现有规则变成 `CAPABILITY_CHANGED`，建议不能再被采用。

## 对调用方

- 建议目录仍从会话响应的 `authoring_capability.fields` 读取。响应模型没有新字段。
- ENCODE 会话看到的字段和问题代码与本版之前相同。
- `hitrac-preprocess` 会话以前没有建议目录。本版起目录里只有 `options.mapq`。
- 前端按返回的字段列表渲染，不需要重新生成客户端。
- 工具读取仍是 `get_authoring_schema` 和 `get_draft_summary`。信封形状不变。
- 越界的提问偏好返回 `INVALID_REQUEST`。越界的建议正文返回 `MODEL_OUTPUT_REJECTED`，草稿不改。
- 采用后的值仍由平台按投影边界再检查一次，然后写入对应文档的键。

## 第一个非 ENCODE 字段

Hi-TrAC 的 `authoring.py` 投影 `options.mapq`。策略边界是 0 到 255，再与当前 `option_schema` 取更紧的整数范围。
未配置运行时时，默认 registry 选择 `HiTracPreprocessAdapter`；运行时获准后选择
`HiTracPreprocessResultsAdapter`。两者读取同一份投影。
`hitrac_preprocess/adapter.py` 和 `results.py` 都留在 Bulk 执行实现清单里，本版不改它们的字节。
校验仍调用现有 `validate_hitrac_inputs`。合法值会带上默认 `threads`，但 `threads` 不是建议字段。

## 产品路径

成员会话不变。路径前缀是 `/api/v1/agent/authoring-sessions`。

1. `POST /api/v1/agent/authoring-sessions`，`workflow_id` 为 `hitrac-preprocess`。
2. `POST /{sid}/connection` 打开建议连接。
3. `POST /{sid}/turns`，意图 `suggest`，字段 `options.mapq`。
4. `POST /{sid}/proposals/{proposal_id}/adopt` 采用合法值。
5. `PUT /{sid}/draft` 修改已采用的值。
6. `POST /{sid}/validate` 重新校验。

本版验收使用平台服务和受控建议正文，不调用模型，也不把模型桩当成在线验证。

默认安装的 Hi-TrAC 没有获准运行时。产品校验在参考绑定里会先跑 adapter 校验：
越界 `mapq` 返回 `HITRAC_OPTION_INVALID`；范围内的 `mapq` 通过这一步后，因为没有运行时，
绑定失败并表现为 `REFERENCE_PROFILE_BINDING_INVALID`。这不是成功快照。
成功快照仍要等获准运行时和服务器拥有的参考配置，本版不提供。

## 本版仍不支持的建议字段

| 工作流 | 未进入建议目录 |
| --- | --- |
| Hi-TrAC | `options.threads`、样本列、空的 config、参考配置选择 |
| ENCODE | 除 `config.threads` 与 `options.cores` 以外的配置、样本计划和参考选择 |
| bulk RNA-seq | 全部字段；`authoring` 模块没有导出 `project_authoring_fields` |

样本计划已由下文的 adapter 契约处理，公共服务不再为 ENCODE 单独分支。

## 样本计划（authoring-sample-plan/1）

registry 选中的 workflow adapter 决定样本计划。公共服务仍负责会话、项目内修订、本地路径解析、幂等、校验缓存和 snapshot。

adapter 实现 `WorkflowSamplePlanAdapter` 时提供四件事：

| 方法 | 作用 |
| --- | --- |
| `supports_workflow_sample_plan` | 是否要求显式样本计划 |
| `normalize_workflow_sample_plans` | assay、layout、role、peak mode 和配对规则 |
| `describe_workflow_sample_inputs` | 按输入使用顺序给出 `WorkflowSampleInputBinding` |
| `materialize_workflow_sample_rows` | 用平台已解析的样本名和路径生成样本行 |

绑定只含样本修订、输入修订和 occurrence。平台据此做项目内核对，并组成已有的 `InputFileRevisionSelection`。occurrence 与 `declare_input_uses` 的读段坐标相同。服务不比较 ENCODE 工作流 ID，也不读取 `reads`。

仍留在共享层的是公开草稿文档：`normalize_draft_content` 限制 `workflow_samples` 的修订身份、读段编号只能是 1 或 2，以及它和 `sample_revision_ids` 的顺序一致。这是所有工作流保存前的文档边界，用来保持已有草稿和幂等，不是 ENCODE 的实验规则。路径是否属于本项目、是否归档、能否打开，也留在平台，因为文件不在工作流字段里。

支持该契约时，缺少计划仍是 `WORKFLOW_SAMPLE_PLAN_REQUIRED`，保存时形状不合法仍是 `INVALID_REQUEST`，计划与修订不符、越项目或已归档输入仍是 `WORKFLOW_SAMPLE_PLAN_INVALID`，路径核验失败仍是 `WORKFLOW_SAMPLE_INPUT_UNAVAILABLE`。输入修订按绑定顺序进入快照。同一修订再次校验仍复用未过期快照。

没有实现该契约的 adapter，包括 Hi-TrAC 和 bulk RNA-seq，不要求样本计划。它们的草稿即使带上 ENCODE 不接受的 `workflow_samples`，仍按各自原有校验处理。方法写在 registry 返回的 adapter 上，不从包里的 `authoring` 模块兜底，因此不必改动已封闭的 Hi-TrAC 文件。

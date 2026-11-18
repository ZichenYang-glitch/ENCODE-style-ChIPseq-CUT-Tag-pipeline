# 外部与下游登记结果：身份、来源和读取契约 v1

状态：IRENE-29 设计契约与可机读例子。基线是 main
`d85c3b6a50a9de0f108691bd60efc678d851a2be`。本文不表示相应表、API 或客户端已经存在；
本交付不写迁移、公开 API、生成客户端或文件解析器，也不登记真实数据。

本文扩展[实验目录 M0](experiment-catalog-m0-contract.md)，沿用
[数据基础 ADR](data-foundation-decision.md) 的 Project、StoragePool、安全相对路径、
不可变修订和 SQLite 权威边界。现有 `Run`、`ArtifactPublication` 和
`downloadRunArtifact` 只覆盖平台运行产物，不能用来冒充外部分析或未由平台执行的下游结果。

## 1. 用户动作与验收

管理员在一个活动 Project 内登记已经位于已绑定 StoragePool 中的一份常规文件。管理员提供
项目内 stable key、显示名、格式、来源、方法、参考和可选父结果；平台只观察并记录文件，
不移动、复制、改名或写入该文件，也不创建 Run。成员随后查询当前结果、查看历史修订，
并按精确结果修订读取字节。

完整合成例子见 [registered-result-read-v1.json](examples/registered-result-read-v1.json)。
后续实现必须用临时 SQLite 和合成常规文件证明以下行为；本 PR 的测试只校验例子，不是运行时验收。

| 编号 | 用户动作 | 必须观察到的结果 |
| --- | --- | --- |
| R1 | 登记一份外部 Hi-C 或单细胞结果 | 得到独立 `rreg_` / `rrev_` 身份；没有 `run_id`、artifact generation、publication 或生命周期事件 |
| R2 | 登记有两个父结果的下游结果 | 两个父边分别保存角色和精确 `rrev_`；不保存“当前父结果”或仅稳定 ID 的弱引用 |
| R3 | 只修改显示名 | 同一 `rreg_` 追加修订；字节摘要、来源、方法、参考和父边不变；旧修订仍显示旧名 |
| R4 | 管理员确认同一逻辑结果已有上游新版 | 追加新修订并移动当前头；已有下游仍指向旧父修订，同时显示“父结果已有新版” |
| R5 | 外部来源、方法或参考没有证据 | 对应字段显式为 `unknown`，界面显示“未知”；没有空字符串、伪造值或猜测的参考基因组 |
| R6 | 成员读取精确修订 | 每次重新鉴权并复核大小和 SHA-256；响应不含物理路径、StoragePool 配置或来源凭据 |
| R7 | 已登记字节缺失或改变 | 读取失败关闭，不返回 200、不自动改挂当前头或另一份同摘要文件 |
| R8 | 执行本设计的最小登记 | 数据库只追加目录元数据；原文件内容、路径和名称不变，不启动 worker 或科学计算 |
| R9 | 有效成员用一个旧修订链接下载，随后撤权并原样重用该链接 | 撤权前为 200，正文匹配旧修订，并有一条带操作者和精确修订目标的成功审计；撤权后同一请求为 401，不产生新的成功下载审计 |
| R10 | 匿名访问下载链接，或有效成员把路径改为另一 Project | 匿名与已撤权会话均为 401；有效成员的错项目请求为 404 `catalog_not_found`，两种失败都不读取源字节或写成功审计 |
| R11 | IRENE-39 合成样例在响应私有副本校验完成后，让外部工具原位改写源文件 | 已建立的 200 只发送私有副本，完整正文 SHA-256 仍等于请求修订；若源变化使副本在校验完成前不符，则在响应前返回 409 |

## 2. 身份与修订语义

新增两个 workflow-neutral 对象，ID 使用小写 `uuid4().hex`：

| 对象 | 前缀 | 形态 | 含义 |
| --- | --- | --- | --- |
| RegisteredResult | `rreg_` | 前缀 + 32 位十六进制 | 一个 Project 内由管理员维护的逻辑结果系列；与 Run 无关 |
| RegisteredResultRevision | `rrev_` | 前缀 + 32 位十六进制 | 一次不可变的描述、来源、父边及字节观察；所有读取坐标 |

`registered_results.stable_key` 在 Project 内唯一，规则与 Sample/InputFile stable key 相同。
显示名不是身份。相同 SHA-256 不去重或合并两个登记；不同来源、角色或科学含义可以共享字节。

- 改名、来源补充、方法/参考修订、父边修订都追加 RegisteredResultRevision。
- 路径、大小或内容摘要改变也必须追加修订。不得原地更新旧修订。
- 上游换版只有在管理员确认它属于同一逻辑结果系列时才沿用 `rreg_`；另一数据集或另一
  科学对象必须新建 RegisteredResult。平台不从 accession、文件名或显示名自动判断。
- `created_at` 是平台提交修订的时间，不冒充上游生成时间或现有产物的 `published_at`。
- 归档只阻止新修订和新绑定，不删除历史读取证据；归档本身不伪造成一个新修订。

修订摘要方案为 `sha256-framed-registered-result-revision-v1`。长度前缀算法沿用现有目录
修订，但使用本方案的独立字符串。摘要覆盖 Project、结果及修订 ID、修订号、完整规范化
payload、私有存储坐标、大小、内容 SHA-256，以及有序父边和父修订摘要。父边关联表是这份
payload 的查询镜像，不能独立修改。

## 3. 规范化 payload 与未知字段

每个修订都包含 `schema_version: 1`、`display_name`、有界格式 token、`origin`、`method`
和 `reference`。模型 `extra=forbid`；未知字段不能用 `null`、空字符串或省略必需键表示。

### 3.1 外部结果

```json
{
  "origin": {
    "kind": "external",
    "source": {
      "status": "known",
      "namespace": "NCBI GEO",
      "accession": "GSE000000",
      "version": "original-2026-01"
    }
  }
}
```

`known` 来源必须同时具有 namespace、accession 和上游 version；缺证据时使用
`{"status": "unknown"}`，且不保留这三个正向字段。v1 不接收任意 URL、带凭据的 locator、
本地绝对路径或自由文本来源。accession 只是所声明的外部来源证据，不会自动变成平台身份。

### 3.2 下游结果与多父边

```json
{
  "origin": {
    "kind": "downstream",
    "parents": [
      {"result_revision_id": "rrev_11111111111111111111111111111111", "role": "primary"},
      {"result_revision_id": "rrev_55555555555555555555555555555555", "role": "control"}
    ]
  }
}
```

下游结果至少有一个父边。父边是有序数组，`role` 是有界 safe-key token，由登记人显式给出；
平台不从 Rep1、control、文件名或目录名推断角色。每个父修订最多出现一次，父修订必须已经
存在、属于同一 Project，并核对其不可变摘要。新修订只能引用已经存在的修订，因此不能形成
环；引用同一结果系列的旧修订是允许的。

若无法指出任何父结果，应按实际情况登记为 `external + source=unknown`，不能创建空的
`downstream` 并声称来源完整。v1 不把 InputFile、Run 或 ArtifactPublication 混进父边类型；
这些关系若将来有需要，必须另加显式版本化 union，而不是添加可选 `run_id`。

### 3.3 方法与参考

`method` 只能是：

- `{"status": "known", "name": "Juicer", "version": "1.6"}`；或
- `{"status": "unknown"}`。

`reference` 只能是以下一种：

- `{"status": "known", "kind": "reference_profile_revision", "reference_profile_revision_id": "refpr_…", "public_identity_sha256": "…"}`；
- `{"status": "known", "kind": "external", "namespace": "NCBI Assembly", "accession": "GCA_000001405.15", "version": "GRCh38"}`；
- `{"status": "unknown"}`；或
- `{"status": "not_applicable"}`。

ReferenceProfile 引用必须钉住 `refpr_` 修订及其公开身份摘要，不能读取当前 enabled 修订代替。
外部参考必须有明确 namespace、versioned accession 和 version；单独的 `hg38`、`mm10`、物种、
项目标签或 Viewer 显示名不构成 `known`。v1 对方法/参考采取保守的整字段 unknown，不把知道
软件名但不知道版本包装成完整已知方法。

## 4. 公开读取投影

默认列表只返回每个 RegisteredResult 的当前头。历史列表按 `revision_number` 升序；精确详情
始终解释请求中的 `rrev_`，不混入当前头的名称、来源、方法、参考、父边或字节。

公开修订 DTO 的有界形状如下，完整值见例子文件：

```json
{
  "project_id": "prj_11111111111111111111111111111111",
  "result_id": "rreg_aaaaaaaaaaaaaaaaaaaaaaaaaaaaaaaa",
  "result_revision_id": "rrev_11111111111111111111111111111111",
  "revision_number": 1,
  "revision_status": "current",
  "display_name": "GM12878 Hi-C contact map",
  "format": "hic",
  "origin": {"kind": "external", "source": {"status": "known"}},
  "method": {"status": "unknown"},
  "reference": {"status": "unknown"},
  "bytes": {"size_bytes": 1048576, "content_sha256": "64 lowercase hex characters"},
  "created_at": "2026-11-18T00:00:00Z"
}
```

`revision_status` 只是读取时由当前头计算的 `current` / `superseded` 注释，不是运行代际状态。
DTO 不包含 `run_id`、artifact generation/revision、StoragePool ID、配置键、相对或绝对路径、
文件 URI、私有 payload 或异常文本。上述示例中的缩略 source 和摘要字符串只用于说明形状，
机读例子使用完整合法值。

下游详情的每个父边返回钉住的 `result_id`、`result_revision_id`、该旧修订自己的
`display_name` 和 `role`，并附加：

- `revision_status`：父修订相对其当前头的状态；
- `current_result_revision_id`：父结果当前头，仅供提示和显式导航。

因此父结果改名不会重写下游历史详情里的旧显示名。父结果换版后，界面显示“父结果已有新版；
本结果仍基于修订 N”，不能自动改边、重新计算方法或把旧字节标成当前父版。

## 5. 字节读取与权限

建议由目录负责人在共享 API 中冻结下列 `/api/v1` 路由和 operation ID；本 PR 不实现它们：

| 方法与路径 | operation ID | 行为 |
| --- | --- | --- |
| `GET /catalog/projects/{project_id}/registered-results` | `catalogListRegisteredResults` | 当前头列表；默认 20、最大 50，沿用目录过滤游标规则 |
| `GET /catalog/projects/{project_id}/registered-results/{result_id}/revisions` | `catalogListRegisteredResultRevisions` | 同一结果的历史修订 |
| `GET /catalog/projects/{project_id}/registered-result-revisions/{result_revision_id}` | `catalogGetRegisteredResultRevision` | 精确修订详情及父边的新版提示 |
| `GET /catalog/projects/{project_id}/registered-result-revisions/{result_revision_id}/content` | `catalogDownloadRegisteredResultRevision` | 精确修订的整文件字节流 |

列表、历史、详情和 `/content` 都使用 `require_principal`；登记、修订和归档使用
`require_administrator`、现有 CSRF 及事务内角色复核。尤其是 `/content` 必须在每次请求重新
校验当前会话、账号启用状态，以及路径 Project 与精确修订的所属关系；不能复用列表或详情时
留下的授权结论。匿名、会话已撤销或账号已停用均为 401，即使请求只是重用一个过去有效的旧
链接。有效成员请求错项目或不存在修订均为 404 `catalog_not_found`，并在解析任何存储坐标前
结束。成员可读共享项目，但目录可见性不扩大其他平台产物的权限。

字节读取必须带精确 `rrev_`，没有“下载最新版”的隐式坐标。授权通过后，服务从该修订取得
StoragePool 和安全相对路径，通过 operator-private 配置定位根目录，以 no-follow descriptor
打开常规文件。服务随后在本次响应私有的临时目录中新建权限为 `0600` 的常规文件，从源
descriptor 复制字节并同时计算大小和完整 SHA-256。只有私有副本与登记证据完全相符后，才
rewind 该副本并准备 200；响应正文只从这个已验证副本发送，不再读取源 descriptor。
副本在响应结束或任一失败后清理，其路径不进入 SQLite、公开 DTO、审计或日志。它是单次
响应缓冲，不是新的长期字节位置、发布或备份。

因此，外部工具即使在私有副本校验完成后原位改写源文件，也不能改变已经建立的响应；任何
成功响应的完整正文都与请求修订的 SHA-256 一致。源文件缺失、类型不安全，或复制结束时大小/
摘要不同，都在发出 200 前返回 409 `registered_result_bytes_unavailable`，不暴露路径或异常。
v1 与现有普通结果下载一样使用整文件响应和 `Accept-Ranges: none`；Hi-C/单细胞格式解析、
Viewer 轨道和 Range 均需独立消费者契约，不能由此端点宣称已经支持。

成功下载沿用[认证 ADR](lan-authentication-decision.md) 的封闭 SecurityAuditEvent，新增
`registered_result.download` action 和 `registered_result` resource kind。事件记录当前
`actor_user_id`、UTC `occurred_at`、`outcome=succeeded`；resource 只能由 bounded factory 从
`(project_id, result_id, result_revision_id)` 这三个逻辑 ID 派生为不透明目标，从而固定到精确
修订，不能使用文件路径、StoragePool 坐标、来源 URL 或任意 context map。只有副本校验和审计
写入都成功才发出 200；审计写入失败时清理副本并失败关闭。401/404/409 不产生成功下载审计。
父结果的新版提示不影响本次读取的授权、审计目标或字节坐标。

## 6. 最小实施交接

以下只是目录负责人整合时的逻辑持久化方案，不是本 PR 的 migration：

| 逻辑表 | 必要内容 |
| --- | --- |
| `registered_results` | `result_id` 主键、`project_id`、项目内唯一 `stable_key`、`created_at`、可空 `archived_at` |
| `registered_result_revisions` | `rrev_`、所属结果、连续修订号、不可变 canonical payload、摘要方案/摘要、`created_at`；私有 StoragePool、安全相对路径、大小和内容 SHA-256 |
| `registered_result_revision_heads` | 每个结果一条当前头，与修订追加在同一事务更新 |
| `registered_result_parents` | 子修订、`ordinal`、父 `rrev_`、父修订摘要、角色；外键 `ON DELETE RESTRICT` |

一次登记或修订在一个 SQLite 事务内复核管理员、活动 Project、Project–StoragePool 绑定、
同项目父修订和期望头修订，然后追加修订、父边镜像和头指针。文件观察沿用 InputFileAccess 的
no-follow、常规文件和前后变化检查，但结果身份不复用 InputFile。期望头冲突返回
`catalog_revision_conflict`；幂等和回执沿用 M0，不另造规则。

实现顺序固定为：先由目录负责人评审并冻结字段/表/路由及上述审计扩展，再写迁移与领域校验，
然后接公开读取和最小管理员工具，最后用 R1–R11 的临时 SQLite + 合成文件验收。生成客户端
只能由该实现 PR 通过仓库命令生成。本工单不直接修改这些共享面。

## 7. 非目标

本契约不创建或补造 Run，不把文件移动、复制或上载到另一受管/长期位置，不设计跨存储迁移或
恢复；第 5 节的响应私有临时副本不属于这些操作。本契约不实现目录网页，不推断样本贡献者，
不自动追随上游最新版，不解析 h5ad/hic，不提供科学有效性背书，也不承诺跨机访问、对象存储、
生产容量或长期保留。恢复演练由 IRENE-26 独立处理。

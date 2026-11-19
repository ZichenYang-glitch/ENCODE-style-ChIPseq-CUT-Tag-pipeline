# 登记结果实现 v1

IRENE-39 把 [读取契约](registered-result-read-contract.md) 落到临时 SQLite。
登记只追加目录元数据。源文件的路径、文件名和字节保持不动，也不创建 Run。

## 写路由

| operationId | 方法与路径 |
| --- | --- |
| `catalogRegisterRegisteredResult` | `POST /api/v1/catalog/projects/{project_id}/registered-results` |
| `catalogReviseRegisteredResult` | `POST /api/v1/catalog/projects/{project_id}/registered-results/{result_id}/revisions` |
| `catalogArchiveRegisteredResult` | `POST /api/v1/catalog/projects/{project_id}/registered-results/{result_id}/archive` |

读路由保持契约中的四个 operationId。写请求要求管理员、CSRF 和 `Idempotency-Key`。
相同键和相同规范化请求返回原来的修订。相同键配不同请求返回 `409 catalog_idempotency_conflict`。
`expected_revision_number` 与当前头不一致时返回 `409 catalog_revision_conflict`，并且不写回执。

客户端不提交摘要。服务用已有的不跟随打开方式计算 `size_bytes` 和 `content_sha256`。
只改显示名时重新观察同一池内相对路径；大小或摘要变化则 `409`，不会静默变成新字节。
调用方显式给出新的 `pool_relative_path` 时，才按新文件登记。

缺来源、方法或参考保持 `{"status":"unknown"}`。参考 pin 必须等于该 `refpr_` 修订已存储的身份，而不是当前启用修订。

## 表与迁移

`20261119_39` 的父修订是主干当前唯一头 `20261119_18`（实验 snapshot 绑定，已经合入，不改写）。它创建：

- `registered_results`
- `registered_result_revisions`
- `registered_result_revision_heads`
- `registered_result_parents`
- `registered_result_write_receipts`

同时把 `registered_result.download` 和资源种类 `registered_result` 加入审计约束。
成功下载的资源身份只由 `(project_id, result_id, result_revision_id)` 经现有封闭工厂派生。

同一条链是 `20261118_17`（实验目录）之后的 `20261119_18`（snapshot 绑定）。IRENE-32 在自己的分支上保留 `20261119_32`，父修订仍是 `20261118_17`。本修订不占用 18 或 32。

## 下载

每次 `GET .../content` 重新检查会话、账号启用状态，以及路径中的 Project 是否就是该修订的 Project。
匿名、撤销或停用账号是 `401`。Project 不匹配、修订不存在、Project 已归档是 `404 catalog_not_found`，此时不解析存储坐标。
字节缺失或摘要不符是 `409 registered_result_bytes_unavailable`。这些失败都不写成功审计。

`200` 只在两件事都完成后发出：响应私有副本（目录 `0700`、文件 `0600`）在复制过程中完成摘要核对，并且该副本的只读句柄已经打开；然后才写入成功审计。响应持有这个句柄，不再重新打开文件。句柄打不开是 `409`，成功审计增量为 0。这个故障只注入在已校验副本的打开函数上，核对钩子必须先执行。审计写入失败时关闭句柄、删除副本，并且不返回 `200`。公开的 `200` 只有 `application/octet-stream`。生成客户端的 `catalogDownloadRegisteredResultRevision` 使用已有的 `blobFetcher`，成功类型是 `Blob`；错误响应仍按 JSON 契约解码。响应头还有 `Accept-Ranges: none`、`Cache-Control: no-store`，文件名只使用 `registered-result.{format}`。副本路径不进入 SQLite、公开 DTO、审计或日志。副本在响应结束后删除，不是新的长期位置。

确定性 Python 分片只安装 Python 依赖，默认验收不启动前端工具。真实生成客户端下载的显式入口是，在 `frontend` 用锁文件执行 `npm ci` 之后：

```bash
HELIXWEAVE_REGISTERED_RESULT_CLIENT=1 python -I -S scripts/checkout_bootstrap.py \
  --repository-root . pytest test/api/test_registered_results_http.py
```

它调用锁文件安装出的 vitest 1.6.1 可执行文件，不通过 `npm exec` 另行解析版本。

因此，核对完成后源文件被原地改写，已经密封的响应体不变，其 SHA-256 仍等于所请求修订。核对完成前源文件变了，则返回 `409`，没有 `200`。

## 第一份外部文件

`test/fixtures/irene39/fixture39.tsv` 是合成单细胞簇计数。登记格式为 `tsv`，来源、方法和参考都是 unknown。验收解析下载体，确认 `synthetic-A=3`、`synthetic-B=2`，并核对服务器计算的字节摘要。这个夹具不是 h5ad，也没有查看器。

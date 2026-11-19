# 目录 CLI：预览实验 assay 分类（IRENE-28）

用户先从同一项目的目录中记录实验 ID 和当前实验修订 ID，再填入想采用的 assay 词条修订 ID。拿到 TSV 后，运行：

```bash
helixweave catalog preview --api-origin http://127.0.0.1:8000 \
  --project-id prj_00000000000000000000000000000001 classifications.tsv
```

预先在 `HELIXWEAVE_CATALOG_SESSION_TOKEN` 中安全提供该服务的有效会话 cookie **值**（不是密码、完整 Cookie 请求头或 CSRF token）。不把值写进命令、TSV、日志或 issue。客户端直接连接 HTTPS 或数字 loopback HTTP 地址，选择对应部署模式的会话 cookie 名，不使用环境代理或跟随重定向；未提供会话时请求仍由服务器鉴权。

## 输入与输出

采用 UTF-8 TSV，表头与列顺序固定。以下是格式示例，ID 均为占位值；实际验证使用独立测试目录实时返回的 ID。

```tsv
experiment_id	expected_experiment_revision_id	assay_term_revision_id
exp_11111111111111111111111111111111	expr_22222222222222222222222222222222	trmv_33333333333333333333333333333333
```

一次只读一个活动项目。本命令沿用 [M0 契约](../architecture/experiment-catalog-m0-contract.md) 与 Slice A `438a943`：按 `experiment_id` 精确匹配，按 `expected_experiment_revision_id` 比较目录列表里的实验修订；不按名称、词条别名或行序匹配。只调用项目实验列表和实验修订列表，显式使用 `revision_view=linked`，沿用服务端默认页长并跟随 `next_cursor`；不请求仅管理员可读的修订明细。

成功读取后输出一个 JSON 对象，`rows` 保留 TSV 顺序和源行号。对于已匹配实验，报告读到的当前实验修订和当前 assay 词条修订（没有已钉住的 assay 则为 `null`）。每行状态如下：

| 状态 | 含义 |
| --- | --- |
| `matched` | 待比较的实验修订就是列表中的当前修订；可比较当前与建议的 assay ID，不表示建议可以提交 |
| `duplicate` | 同一实验 ID 在表里出现多次；报告所有出现，不选第一行或最后一行 |
| `unmatched` | 该 ID 不在当前项目的可见活动实验列表中；不会跨项目查询 |
| `stale_revision` | 预期修订出现在同一实验的可见历史中，但不是列表中读到的当前修订 |
| `unknown_revision` | 实验匹配，预期修订却不在该实验的可见历史中；不把别的实验的修订认作旧版本 |

行级结果是预览数据，退出码为 `0`。输入格式错误退出 `2`，HTTP、传输或响应契约失败退出 `1`；失败只在 stderr 输出安全诊断，不输出部分预览。跨页读取不构成快照，不能代替任何写入时的并发检查。

## 验收

- 已安装的 `helixweave catalog preview` 用有效会话、单行 TSV 和真实服务端默认页长运行，退出 `0`，报告 `matched` 及正确的当前实验修订。
- 小 TSV 覆盖当前、重复、未匹配、已知旧修订和外部修订；保持顺序。真实 HTTP 的后续页能够命中实验和旧修订。
- 匿名收到 `401`；成员和管理员可得到相同的列表预览，且成员不能读取旧修订明细。预览仅发 GET 请求，报告不含会话值、服务端错误正文、候选或别名。
- 不改变目录服务、迁移、生成客户端或写规则。只比较 assay 的修订 ID，不校验词条可采用性，也不支持样本分类或其他实验字段。显式提交、写入冲突处理和回执是单独的后续契约。

从已按 [源码验证指南](source-provenance.md) 安装本 checkout 的环境运行：

```bash
python3 -I -S scripts/checkout_bootstrap.py --repository-root . pytest \
  test/cli/test_catalog_preview.py test/cli/test_primary_cli.py \
  test/api/test_experiment_catalog_http.py -q
```

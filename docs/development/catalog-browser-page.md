# 实验目录网页

状态：本分支的页面说明。基线是 Slice A `438a943f158a1dba0943d5c3defa3eb4531d2d5b`
（PR #225，`agent/hw-catalog/5cd11aff828c`）。A 尚未合入主干，因此本分支堆在该
提交上，Draft PR 的 base 是 `agent/hw-catalog/5cd11aff828c`，不是 `main`。

本页只调用已经生成的目录客户端和现有的课题列表。不改目录 API、迁移、
`frontend/src/api/generated/`，也不改 WashU 面板。

## 页面入口

登录后的主导航增加 Catalog，路径是 `/catalog`。成员和管理员都能打开。
页面先选择课题，再按 assay 词条和细胞类型词条查找实验。同一字段多选是 OR，
不同字段是 AND。留空不是“查未知”。结果按目录契约每页 20 条，用返回的
cursor 翻页。点开一条实验后，用修订列表查看当前头修订和更旧的修订；画面上的
名称来自该修订钉住的词条，不改成词条的当前头名称。

## 角色

| 角色 | 页面 |
| --- | --- |
| `member` | 查找、翻页、打开旧修订、阅读实验追溯和样本运行范围 |
| `administrator` | 成员能做的全部，另外可以录入词条、用已有样本修订录入实验，以及修订当前头修订的显示名和 assay |

成员页面不渲染录入和修订表单。写请求仍由目录 API 拒绝非管理员。修订是追加新的
实验修订；正在查看的旧修订不会被改写。头修订上的 target、batch 和样本钉原样带入
新修订。歧义栏不在本页代选候选。

## 请求

读：

| 动作 | 生成客户端 |
| --- | --- |
| 课题 | `authoringListProjects` |
| 词条选项 | `catalogListTerms` |
| 实验查找 | `catalogListExperiments`，`limit=20`，可选 `assay_term_id`、`cell_type_term_id`、`cursor` |
| 旧修订 | `catalogListExperimentRevisions` |
| 实验追溯 | `catalogListExperimentRevisionRunResults` |
| 样本运行范围 | `catalogListSampleRevisionRunResults` |

管理员写，每次带新的 `Idempotency-Key`：

| 动作 | 生成客户端 |
| --- | --- |
| 录入词条 | `catalogCreateTerm` |
| 录入实验 | `catalogCreateExperiment` |
| 修订当前头修订 | `catalogReviseExperiment` |
| 选择样本修订 | `catalogListSamples` |

默认修订视图是服务端的 `linked`。本页不传 `revision_view`。

## 三条验收

1. 管理员可以录入实验分类，或修订当前头修订的显示名和 assay。成员看不到这两个表单。
   修订请求带当前头修订号；旧修订的显示名仍可打开。
2. 成员选定课题、assay 和细胞类型后查找实验，每页 20 条，cursor 翻到下一页，
   并能打开比头修订更旧的一版。
3. 实验追溯在 `trace_status=unresolved` 且
   `reason=experiment_snapshot_binding_absent` 时显示“未绑定”，不把空结果画成产物。
   样本运行结果若存在，标成 `run_input_scope` / “运行输入范围”，关联样本写成运行样本，
   不写成文件贡献者。

## 本页不做、接口仍在

- 不新建课题，也不新建逻辑样本。没有样本修订时，录入实验的表单说明这一点。
- Slice A 的实验追溯固定为未绑定。本页按响应显示，不补做 Slice B 的绑定。
- 不调用 `catalogGetExperimentRevision`、`catalogGetTerm`，因此不展示候选或别名。
- 不调用归档路由，也不调用 `catalogReviseSampleCatalog`。细胞类型筛选读的是已经钉在样本修订上的分类。

## 浏览器路径

`node frontend/scripts/catalog-page-browser.mjs` 用临时 SQLite 和真实 HTTP 走完三条验收：
成员查找、20 条分页、旧修订、未绑定；管理员录入词条、修订当前头修订、录入实验。
Python 解释器可用 `HELIXWEAVE_PYTHON` 指向能导入本仓库 `src` 的环境。脚本不进入默认 pytest。

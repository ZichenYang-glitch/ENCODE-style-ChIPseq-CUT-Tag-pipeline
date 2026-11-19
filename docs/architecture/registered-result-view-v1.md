# 已登记 Hi-C 与单细胞结果的查看候选 v1

状态：IRENE-40 的查看需求与验收设计。查看器基线是 main
`f52a2c20adae6c6e61375b34e47211c84aa035f7`。外部登记的精确修订是 IRENE-39
验收提交
`29d508647c18df94016f94d7da437e794b846f9e`
（[PR #239](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/239)）。
该提交不在本分支上。本轮不新增 adapter，不增加路由、迁移、页面或执行器。
实验室选定真实动作之前不实现其中任何一个候选。

本文是登记读取契约所要求的独立消费者契约。身份、精确修订和整文件下载沿用本分支的
[读取契约](registered-result-read-contract.md) 与
[登记例子](examples/registered-result-read-v1.json)。相对上述验收提交，这份契约只少了指向实现说明的两行；实现说明不在本分支，见
[29d5086 的实现说明](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/blob/29d508647c18df94016f94d7da437e794b846f9e/docs/architecture/registered-result-implementation-v1.md)。
登记例子与该提交中的同名文件一致。轨道分类和页面文句沿用
[WashU 只读结果浏览](../development/washu-browser-read-contract.md)。嵌入页与全页的边界沿用
[目录、Agent 与 Viewer 边界](catalog-agent-viewer-boundaries.md) 第 6 节。
合成例子见 [registered-result-view-v1.json](examples/registered-result-view-v1.json)。

在 `29d508647c18df94016f94d7da437e794b846f9e` 与本查看器基线之间，
`src/encode_pipeline/adapters/browser/catalog.py`、
`frontend/src/features/run-artifacts/GenomeBrowserPanel.tsx`
和 WashU 只读契约的文件内容相同。下面的读取器描述对这两个提交都成立。

## 1. 精确修订

查看坐标是某一次 `rrev_`。列表里的当前头不能代替请求中的修订。旧链接在改名或上游换版之后仍读取该修订自己的字节、显示名、方法和参考。登记内容路由保持整文件响应，`Accept-Ranges: none`。该路由不解析 hic 或 h5ad，也不因此成为查看器。

IRENE-39 的临时库没有把当时的 `rrev_` 写进仓库。该验收提交里能够复核的字节只有合成簇计数文件 fixture39.tsv：50 字节，SHA-256
`d6cb16f73da05e9d2e6f9a70e838dc9955dafdf560a955f6f62bb1dc2907817a`，
`synthetic-A` 为 3，`synthetic-B` 为 2。格式是 `tsv`。方法与参考都是 unknown，页面词是「未知」。该文件不是 h5ad，验收说明里写明它没有查看器。

读取契约例子里的 `external_v1` 使用 `rrev_11111111111111111111111111111111`，格式 `hic`。那里的长度和摘要只说明公开 DTO 的形状，没有一份被解析过的接触矩阵。`missing_source` 的格式是 `h5ad`，也没有对应的查看器。

两个候选互相独立。选定或日后实现其中一个，不自动开始另一个。

## 2. 候选动作

| 动作 | 选择 | 所需格式 | 现有可复用读取器 | 本轮结果 |
| --- | --- | --- | --- | --- |
| `view-hic-region-resolution` | 用已知区域、分辨率和归一化核验一个接触值 | `hic` | 运行产物上的 WashU `hic` 轨道 | unsupported |
| `view-single-cell-qc-cluster-counts` | 单细胞三项里选择 QC 簇细胞数 | 列为 `cluster_id`、`cell_count` 的 `tsv` | `catalogDownloadRegisteredResultRevision` 的整文件字节 | 字节已核对，没有表页面 |

### 2.1 Hi-C：区域、分辨率、归一化

请求带上项目、精确 `rrev_`、区域 `chrom:start-end`、分辨率（bp）和归一化名称。合成请求沿用 `external_v1` 的修订身份，区域 `chr1:0-2000000`，分辨率 `1000000`，归一化 `NONE`。这是设计请求。该修订没有观察到的矩阵，所以本轮不给出接触值。

`src/encode_pipeline/adapters/browser/catalog.py:86` 把文件名以 `.hic` 结尾、或 `output_type` 为 `hic` / `contact_matrix_hic` 的运行产物标成 WashU 类型 `hic`、原始接触矩阵，并标为可显示。`frontend/src/features/run-artifacts/GenomeBrowserPanel.tsx:519` 交给 `GenomeHub` 的轨道只有 `type`、`name` 和 `url`；`viewRegion` 是区域字符串。分辨率和归一化没有入口。同文件把 `kind` 为 `contact_matrix` 的轨道排除在区域记录摘要之外，因此区域摘要也不是接触值。浏览器会话的轨道来自某次 run 的 artifact，坐标里没有 `rrev_`。

预检请求 `Range: bytes=0-0`，并把 HTTP 200 与 206 都当作第一字节可读。这个状态只说明第一字节的 HTTP 结果，不证明矩阵、分辨率或归一化已被读取。登记修订的内容路由是整文件，响应头为 `Accept-Ranges: none`。现有 WashU hic 轨道因此不能对一份登记修订完成这个动作。当前结果是 `unsupported`，原因 `registered_revision_has_no_matrix_reader`。

日后实现且实验室提供了真实 `.hic` 时，成功表现为：

- 嵌入页与全页显示同一个精确修订上的同一个接触值。
- 页面标出的区域、分辨率和归一化就是这次读取实际使用的三个值。
- 接触值绑定该修订的 `content_sha256`。摘要一致之前不画数值。
- 参考为 unknown 时显示「未知」。WashU 染色体表只在已有绑定且程序集标签对得上内置脚手架时用于导航，并且不能代替参考身份。没有绑定就不猜测 hg38 或 mm10。
- 旧修订链接保持该修订的矩阵。当前头换版后，旧链接不改画新矩阵。

不支持时的表现：

| 情形 | 表现 |
| --- | --- |
| `.cool`、`.mcool` | 页面写 Unsupported format。原因码 `unsupported_contact_matrix`。明细沿用目录中的现句：Cooler matrices are not loaded. This view does not run HiGlass. |
| 文件中没有所请求的分辨率 | 页面状态 `requested_bin_absent`。不改用文件里的另一个 bin。 |
| 文件中没有所请求的归一化 | 页面状态 `requested_normalization_absent`。不改用 NONE、VC、KR 或其他归一化。 |
| 会话撤销或匿名 | HTTP 401，页面状态 `no_matrix`。 |
| 错项目或修订不存在 | HTTP 404，页面状态 `no_matrix`。 |
| 字节缺失或摘要不符 | HTTP 409，页面状态 `no_matrix`。不把失败画成接触值 0。 |
| 把 `tsv` 或 `h5ad` 送进这个动作 | 页面状态 `no_matrix`。不把计数解释成接触值。 |

### 2.2 单细胞：QC 簇细胞数

单细胞候选从 QC、embedding、表达中选定 QC。所需格式是上述 50 字节 TSV，列只有 `cluster_id` 和 `cell_count`。embedding 没有读取器。表达所对应的 `h5ad` 在登记契约中明确不解析，因此不选表达。

可复用的读取器是登记修订的整文件下载 `catalogDownloadRegisteredResultRevision`。IRENE-39 已在验收提交上观察到：改名之后，旧链接仍返回这 50 字节，摘要仍是 `d6cb16f73da05e9d2e6f9a70e838dc9955dafdf560a955f6f62bb1dc2907817a`，两个簇计数仍是 3 和 2。该操作返回字节，不把 TSV 画成表。当前结果是 `bytes_verified_table_not_shown`。

日后实现时，成功表现为：

- 所请求修订的正文摘要等于上述 SHA-256 之后，表为 `synthetic-A = 3`、`synthetic-B = 2`。
- 嵌入页与全页是相同的两行。摘要一致之前不画细胞数。
- 方法与参考显示「未知」。
- 旧修订链接保持这两行。改名不改这两行的计数；换了字节的新修订使用自己的摘要，不回写旧链接。

不支持时的表现：

| 情形 | 表现 |
| --- | --- |
| `h5ad`、embedding 坐标或基因表达表 | 页面状态 `no_reader`。不画空的 embedding，也不画全零表达。 |
| 列不是 `cluster_id` 与 `cell_count` | 页面状态 `unsupported_table`。不把其他列解释成线粒体比例、UMI 或基因数。 |
| 会话撤销或匿名 | HTTP 401，页面状态 `no_table`。 |
| 错项目或修订不存在 | HTTP 404，页面状态 `no_table`。 |
| 字节缺失或摘要不符 | HTTP 409，页面状态 `no_table`。不把失败画成细胞数 0。 |

## 3. 加载时序

本轮没有新页面，因此没有新的浏览器验收。日后被选定的那个动作沿用 Viewer 边界第 6 节：同源嵌入页与全页，查看器不查 SQL，不提交任务，不启动 worker 或 Peri。宿主只传版本化会话和精确 `rrev_`。

两个表面的顺序相同：

1. 先显示修订身份：显示名、格式、方法、参考、`rrev_`、字节长度和 SHA-256。
2. 正文摘要与该修订一致之后，才显示接触值或细胞数。
3. 401、404、409 结束等待状态。页面上不保留上一次成功的接触值或细胞数。

嵌入页不另选全页没有提交的分辨率或归一化。全页也不读取嵌入页会话之外的当前头。

## 4. 本轮验收与非目标

本轮验收只检查合成例子与本文一致：两个候选的格式、读取器、当前结果、成功形状和不支持表现都写在例子里；例子不出现 `run_id` 或存储路径；Hi-C 例子不编造接触值；单细胞例子的 TSV 摘要等于 IRENE-39 已核对的 50 字节。检查命令：

```bash
python -I -S scripts/checkout_bootstrap.py --repository-root . pytest \
  test/docs/test_registered_result_view.py test/docs/test_internal_links.py
```

非目标：不新增 `adapters/` 下的执行 adapter，不把 cool、mcool 接到 HiGlass，不解析 h5ad，不从簇计数推断 embedding 或表达，不把运行产物的 artifact generation 写成登记修订，也不宣称 WashU 已经按分辨率和归一化核验过一份登记矩阵。

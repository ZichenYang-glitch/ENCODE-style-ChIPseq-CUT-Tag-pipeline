# WashU 轨道浏览

本说明记录内嵌查看器的公开读取契约。运行产物的下载、科学执行和参考注册分别沿用各自的入口。

## 从本机开始或叠加已有结果

首页的「打开本地轨道」进入 `/viewer.html`。用户可以在同一个窗口里选择参考组装、打开文件、按基因或坐标定位。由现有运行进入时，选轨从该运行开始。两种入口都允许继续添加文件和有权访问的其他运行。帧内查看保留同源宿主握手。

| 用户选择的来源 | 读取与身份 | 出现问题后 |
| --- | --- | --- |
| 本机 WashU longrange | 一次选择 BGZF 数据和 TBI 索引，按基因组窗口读取。显示原文件名和用户指定的组装 | 缺少索引时提示补选；有多个可能配对时由用户指定 |
| 本机 bigWig、普通 BED 峰文件 | 浏览器会话直接读取选中的文件。文件内容与本机路径留在浏览器 | 给出不支持的格式或已知参考冲突，保留其他轨道 |
| 已授权运行中的产物 | 每个运行分别读取自己的会话和版本化地址。来源包括运行、产物、generation、revision 和参考绑定 | 某条轨道失效时标明原因，由用户决定是否重新选择 |

本地轨道的“用户指定 GRCh38”记录使用者的判断。平台运行的参考证据来自各自的绑定。平台公开的参考修订摘要还包含 profile、名称和工作流信息。相同组装下的运行可以使用不同参考修订，查看时保留双方的身份。染色体名称只用于导航，已知组装冲突会阻止叠加。文件和本机路径留在浏览器。刷新页面后须重新选择本机文件。

同一搜索框既接受基因名，也接受 `chr:start-end`。选择唯一命中或歧义列表中的某个位置后，所有已选轨道使用同一个窗口。查无匹配显示提示，注释不可用时仍可直接输入坐标。

## 组件

- npm：`wuepgg@62.0.2`（精确版本），样式 `wuepgg/style.css`
- 实际导出：`GenomeHub`、`GenomeViewer`、`useSessionData`
- 首版只用 `GenomeHub`。已发布的 `GenomeViewer` 在读取远程轨道时调用未定义的
  `formatDataByType`，不能完成绘制。`wuepgg3-track@9.0.3` 随依赖安装
- 同伴依赖：React 18 或 19。本仓库前端是 React 18.3.1
- 许可证：华盛顿大学非商业许可，全文在
  `frontend/src/features/run-artifacts/wuepgg-LICENSE.txt`
- 已核对的 62.0.2 包行为，而不是文档站示例：
  - `GenomeViewer` 只接收 `genomeName`、`tracks`、`viewRegion`、`windowWidth`
  - 文档中的 `customGenome` 和 `indexUrl` 没有接到这个组件
  - 基因组名必须是包内染色体表的键，例如 `hg38`、`mm10`
  - 轨道模型不接收单独的索引 URL；tabix 索引请求是数据 URL 再加 `.tbi`
- 文件字节用 XHR。同域请求默认携带 Cookie；组件没有把
  `withCredentials` 打开，因此不能靠跨域或查询参数令牌
- 轨道地址必须是同域绝对 URL。相对地址失败时，组件会把路径改写到
  `https://epigenome.wustl.edu/cors/`。首版不使用该代理
  - 内置 `geneAnnotation` 会访问外部注释服务。首版不添加这种轨道
  - `GenomeViewer` 没有提供取消进行中 XHR 的句柄

不部署 MongoDB，也不把 SQLite 换成别的库。

## 公开 HTTP

这些路由挂在现有 artifact 路由器上，因此沿用成员会话和 CSRF。
下载路由 `GET /runs/{run_id}/artifacts/{artifact_id}/download` 保持
`Accept-Ranges: none`。

| 方法 | 路径 | 作用 |
| --- | --- | --- |
| GET | `/api/v1/runs/{run_id}/browser-session` | 参考身份、可显示轨道、明确不支持的格式 |
| POST | `/api/v1/runs/{run_id}/browser-tracks/{artifact_id}/derivative` | 为区间或 PET 生成可追溯的 bgzip/tabix 衍生物 |
| GET, HEAD | `/api/v1/runs/{run_id}/browser-tracks/{artifact_id}/bytes/{generation}/{revision}` | 原产物或衍生物字节 |
| GET, HEAD | 同一路径加 `.tbi` | 衍生物的 tabix 索引 |
| GET | `/api/v1/runs/{run_id}/browser-tracks/{artifact_id}/region-summary/{generation}/{revision}` | 当前窗口的记录数，或说明该轨道不是计数轨道 |

浏览器只提交 run、artifact、generation、revision。响应不包含宿主路径、
环境值或异常文本。generation 或 revision 与当前索引不一致时返回 409，
旧地址不能继续读到字节。

字节响应：

- 无 `Range`：200，`Accept-Ranges: bytes`，`Content-Length` 为文件长度
- 单一 `bytes=start-end`、`bytes=start-` 或 `bytes=-suffix`：206 和
  `Content-Range`
- 越界、空范围或无法解析：416，`Content-Range: bytes */length`，空正文
- 多范围：416，不返回 multipart
- HEAD 与对应 GET 的状态和长度头一致，不返回正文
- 描述符在发送完成、失败或连接中断时关闭

身份或授权失败仍使用 JSON issue。Range 语法失败使用 416，因为 WashU
按 HTTP 字节范围读取，不解析平台 JSON。

大文件若使用 descriptor revision，单次 Range 只读取请求的字节，不从头
扫描文件。不超过 16 MiB 且使用 content revision 的产物，现有下载准备
过程会为了核对内容修订而读完整文件；这是既有修订规则，不是浏览器把整文件
交给页面。衍生物索引自己的 Range 只读取请求区间。

## 科学语义

平台保存身份、授权、generation/revision 和字节范围。adapter 目录
`adapters/browser/` 按 workflow id 和公开 `output_type` 解释轨道，
不在通用 API 里写 Hi-TrAC 路径或 ENCODE 样本列。

| 来源 | WashU 类型 | 含义 |
| --- | --- | --- |
| 已有 bigWig / bigBed 信号产物 | `bigwig` | 原文件。数值是信号，0 是真实值 |
| Hi-TrAC `hitrac_bedpe_all`、`hitrac_bedpe_no_bg` | `longrange` | 展示衍生物。源文件是 PET，不是显著 loop，也不是接触矩阵 |
| 其他 BEDPE | `longrange` | 成对区间。不把源分数或 MAPQ 解释成互作强度 |
| 已有 peak/BED 区间 | `bed` | 展示衍生物。保留来源分数列时，标签写明它是峰文件自己的分数 |
| 已有 `.hic` | `hic` | 原文件接触矩阵 |
| `.cool`、`.mcool` | 不加载 | 首版没有 HiGlass 服务 |

PET 衍生物使用 BED 0-based half-open 坐标。每条 PET 记 1 对。两个锚点
各写一行，第四列分数固定为显示标志 `1`。这个 `1` 不是 MAPQ、PET 计数或
显著性；两行分数不能相加。相同锚点只写一行。`pair_count` 写在衍生物
来源摘要里。

衍生物保留原产物。摘要包含转换器 `washu-display-derivative-1`、源
artifact revision、坐标约定和分数含义。转换在解压和解析过程中停止：压缩或明文源不超过
32 MiB，解压后不超过 256 MiB，单行不超过 64 KiB，源记录不超过 200 万，展示行不超过
200 万。PET 的两条锚点都计入展示行，不改变 `pair_count`。任一上限被打破时转换失败，
不写入部分衍生物，也不覆盖原始科学产物；页面把该轨道标为过大，不把抽样或空白画成完整结果。
256 MiB 足以容纳未超过展示行上限的常规 BED/BEDPE 文本，并在高压缩比输入继续膨胀前停下。
索引身份用文件戳核对。现有回归和浏览器负例都是向 `track.gz.tbi` 追加字节后再读，期望 `409 derivative_changed`；这覆盖追加导致的戳变化，不覆盖同长度、同 inode 的静默改写。
WashU 自身对一个区域最多绘制 10 万行；区域摘要若超过该数，页面说明绘制不完整。

独立验收不使用 CI 的 browser-e2e。在仓库根目录执行：

```bash
export PATH="/home/yangzichen/workflow/ENCODE-style-ChIPseq-CUT-Tag-pipeline/.local/node22/bin:$PATH"
export PLAYWRIGHT_BROWSERS_PATH="/home/yangzichen/workflow/ENCODE-style-ChIPseq-CUT-Tag-pipeline/.local/ms-playwright"
node frontend/scripts/washu-browser-accept.mjs
```

该入口只调用 `frontend/scripts/washu-browser-proof.mjs`：正例要求预检与后续数据和索引读取分开、已知窗口 `chr1:500-2700`、PET/峰/bigWig 的已知字节，以及标尺上的 1000 与 2000，并以 0 退出且打印 `WASHU_PROOF ok`。负例在一字节预检成功后追加索引字节，真实路由返回 `derivative_changed`。预期拒绝以 0 退出，打印 `WASHU_PROOF rejected`，并且不打印 `WASHU_PROOF ok`。其它证明错误以非零退出。汇总只在正例满足上述成功输出、负例满足预期拒绝时通过；负例非零、缺少拒绝文本，或负例结果里出现 `WASHU_PROOF ok`，都会让汇总失败。

没有绑定参考修订时，不猜测 hg38/mm10。有绑定时，页面展示 profile、
revision、identity SHA-256 和 binding digest。WashU 染色体表只在程序集
标签能对应到内置脚手架时用作导航，并单独标明它不是参考文件一致性的证明。
运行注释轨道只展示已有注释类产物。按名称定位使用查看器自己的离线索引，与是否登记了注释轨道分别处理。

## 在 WashU 内选择轨道文件（IRENE-21）

研究者从当前运行打开 WashU，在查看器内找到并勾选要看的文件。实现按设计提交 `fada3ed129dfcd014a3c7825c77c315667faf978`。

### 选择行为

直接打开时显示当前运行的候选文件，初始选择为空。用户可以添加、移除轨道，也可以取消最后一项。若宿主显式传入预选，只在进入该运行时用它初始化一次。同一代次内，翻页、刷新列表或一条轨道准备完成，都保留用户已有的选择。切换运行时重新初始化。

选择绑定运行、产物 ID、代次和修订。候选列表可读取元数据。只有用户选中的文件才会发起轨道预检和数据读取，需要准备索引时由用户触发。格式暂不支持、参考不兼容、版本失效或索引待准备的原因显示在对应文件旁，其他可读轨道仍可查看。

### 候选发现与版本变化

会话复用 `RunService.list_artifacts_page`，按稳定产物 ID 游标遍历有权访问的当前运行产物。第一页读到的代次用于之后的每一页。沿用本契约已有的轨道识别、同源授权和版本化读取。完整扫描同一代次后，才能宣告某类轨道缺失。

若两页之间代次改变，会话接口丢弃已扫描的部分结果，返回轨道读取接口已有的 HTTP 409 及 `generation_changed` 原因。查看器提示刷新和重选，由用户决定何时读取新版文件。

### 验收用例

| 入口及操作 | 预期结果 |
| --- | --- |
| 浏览器直接打开含 A、B、C 的运行，再勾选 A、B | 打开时轨道数据请求为零。勾选后只读取 A、B，C 仅出现在列表中 |
| 移除 A，刷新同代列表，最后移除 B | 刷新后仍只选中 B。取消最后一项后显示空选择，后续不再发起这些轨道的读取 |
| 显式预选 A，用户移除 A 后另一条轨道准备完成，再切到其他运行 | 准备完成不会重新勾选 A，切换运行也不会沿用旧运行的选择 |
| 在临时 SQLite 中登记 500 份普通报告，再登记一份互作轨道，经过真实 HTTP 接口打开 WashU | 能找到并选中第 501 个产物，选中后可以读取该轨道 |
| 在读取两页候选之间切换代次 | 会话接口返回 409。界面给出刷新入口，用户主动重选前不会读取新版轨道 |

打包后的页面由 `frontend/scripts/washu-selection-proof.mjs` 记录前两行和最后一行的浏览器网络请求。第 501 个产物和会话 409 的 HTTP 形状由 `test/api/test_routes_browser_tracks.py` 覆盖。显式预选由查看器组件测试注入，当前宿主打开时不带预选。

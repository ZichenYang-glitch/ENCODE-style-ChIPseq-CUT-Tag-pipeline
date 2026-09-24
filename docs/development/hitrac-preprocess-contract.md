# Hi-TrAC H1：固定预处理契约与 tiny 基线设计

状态：H1政策A及中间产物保留政策已获用户确认，H2已独立验收。
H3早期并发返修、H4结果接线与F1/F2限定返修均已获独立验收。
H5真实平台全链已交付，核心功能已经独立复验（未发现新增生产缺陷）；
限定测试/文档收尾并入最终候选。最终本机Gate在最终clean commit上执行，
状态以该轮交付报告与roadmap最新记录为准；远端protected CI未执行。
第6～9节保留各阶段的设计、停止记录及验证范围；第10节描述当前H4结果契约；
第11、12节分别保留F1/F2返修与H5验收当时的验证边界。
实施次序以[roadmap](workflow-platform-agent-roadmap.md#后续-hi-trac-预处理-adapter五个-pr)
为准。

## 1. 上游与接入范围

科学范围严格等于 cLoops2 **`tracPre2.py` 预处理及其已有 BEDPE/QC**。
独立 workflow ID 拟为 `hitrac-preprocess`；不增加 loop/domain/peak calling、
bigWig、差异分析或新算法。薄启动器只负责输入/身份准入、原脚本执行、失败识别和
结果发布，不重写 trimming、比对、去重或 QC。

本轮重新核对官方 tag `v0.0.5` 指向 commit
`de6cc732fa00b408551b9f4272933640c08447f1`；[原脚本][trac] 16,962 字节，SHA256
`c3c4ef4e6287fa4ade97a6f5980d20c81ea345b8e8e13c88ae3b7c4bd3a67aec`。
同名 PyPI `0.0.5` 不构成相同源码身份。H2 须锁定整个所需包/CLI 闭包和 Python、
pandas、joblib、Biopython、Bowtie2、samtools、bedtools、gzip 的实际构建产物与摘要。
[setup.py][setup] 未声明 Biopython，[utils.py:36][utils] 使用 distutils；旧环境描述
不是可直接使用的锁定 runtime，固定组合及安装字节锁见第7、8节。分发保留
[BSD-3-Clause 许可][license]；尚未完成全部工具的许可核验。

## 2. 输入、参数与参考

以下区分原行为与拟议平台校验；不是当前已经发布的 schema。

| 对象 | 固定源码行为 | 首版接入设计 |
| --- | --- | --- |
| 样本/多 lane | `tracPre2.py:108–129` 发现一级目录的 `_R1.fastq.gz` / `_R2.fastq.gz`；占位列表长度检查不能发现缺 mate | 多样本，每样本一对 gzip FASTQ；按已有 roadmap，用户先合并多 lane，首版不自动合并、不分 lane 去重。各样本独立去重和 QC，不混样本 |
| 配对完整性 | `:185–188` 以 zip 遍历，不比数量/ID；`:214–216` 改写双方 read name | 先完整核验 gzip/FASTQ、SEQ/QUAL 长度、两端记录数和逐对 ID。ID 建议只规范化首个空白前 token 及 `/1`、`/2`，存在 Illumina mate 标识时核对角色；不排序/补齐/丢弃配对，不新增按 read ID 去重 |
| 安全映射 | 样本名既进入 shell，又被 `_all`、`_unique` 拆分 | 私有 sample→R1/R2 映射，内部名 `s000001` 等；受控目录/参考别名不含空白或 shell 元字符，原数据只读。不得直接限制为“用户必须重命名原始文件”或开放任意 shell 参数 |
| linker | `:140–154,170–223` 固定 `CTGTCTCTTATACACATCT`，匹配其或 RC 的 9bp 前缀，取首个命中；两端剩余长度均≥10才保留 | 不开放另一 trimming 算法/linker。保留末尾恰好9bp seed 起点未被扫描的现有边界，H2 单独验证 |
| 并行 | `:79–99` 的 n=1、p=5；`:427–429` 转换 worker 数为 `int(min(样本数,n*p/2))`，samtools 固定 `-@2` | H2 资格入口 n 固定1、p 默认2并限定2～8，拒绝零 worker；原脚本默认5不变。实际工具辅助线程另计，参数不是硬 CPU/RSS 限制 |
| MAPQ | 默认10；`:299` `samtools view -b -F 4 -@ 2 -q N` 对 alignment 筛选，随后 name-sort、`bamToBed -bedpe` | 保留默认10；H2 私有入口允许0～255的整数阈值，不增加 proper-pair/blacklist/UMI 政策。单端未比对/低 MAPQ 导致的不完整 BAM 配对由 H2 实工具验证，不先假定每 raw pair 都产生 PET |
| 参考 | `:372` 只检查 `prefix*.bt2` 至少一项 | 初版准入完整 `.1/.2/.3/.4/.rev.1/.rev.2.bt2` 六文件、FASTA 摘要、contig 集、构建命令与 Bowtie2 身份；只读且执行前复核。仅 `.bt2l` 不被当前原入口接受，暂不支持；不在线补建 |

配对校验和参考闭包是 HelixWeave 的接入边界，不是上游已经保证的性质。
空原始输入及执行后零结果的处理见第5节。先沿用可验证的外部输入模式；
`services/managed_input_verification.py:68–83` 尚拒绝 managed uses 的执行，不能
借新 adapter 自动解除总闸。`platform/adapters.py:653–661` 的 WorkspacePlan 不是
现成 FASTQ symlink 接口；受控映射在后续启动器设计中落实。

## 3. PET、输出和 QC 单位

[tracPre2.py:259–287][trac] 先去背景：cis **浮点中点距离 `<1000` 且两端均无 linker**
才丢弃；trans 两端均无 linker 也丢弃。随后以原 BEDPE **前六列字符串 tuple 的 hash**
去重，保留首先遇到的一行；不含 strand/name/MAPQ，不等同 Picard/samtools 的重复键。

[qc.py:66–107,140–154][qc] 的 TotalPETs 是可解析 BEDPE 行数，含重复；QC 去重键则经
[ds.py:44–86][ds] 规范化端点且包含 strand。cis/trans 及距离分布基于该 QC 独特集合。
QC 距离先取整数中点：close ≤1000、middle 1000<d≤10000、distal >10000；不要把
背景 `<1000` 和 QC `≤1000` 混成一条规则。

计划正式产物仅为各样本原 `*_all.bedpe.gz`、`*_unique.bedpe.gz` 和一份原
`tracPre_summary.txt`。BEDPE 按原十列端点/name/score/strand 校验；区间采用
[bedtools BEDPE 约定](https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format)
的0-based start、半开区间。该格式说明不是工具版本锁；本轮 H2 实际构建和证据见第7节。
BAM 是 name-sorted，不承诺坐标排序或 BAI。按已确认政策，原脚本留下的 BAM、裁剪 FASTQ 保留
在私有 attempt，首版不公开下载也不自动删除；原脚本会删除 SAM、压缩前 BEDPE 及
中间 `*_bedpeQc.txt`，不能承诺这些仍留存。原始 FASTQ 永不由此接入清理。
日志含私有路径，只经现有受控脱敏入口消费。

**最终 summary 是15个指标列，加独立样本索引列，共16个物理 TSV 字段。**
`tracPre2.py:467` 曾先写中间表，`:493–511` 才产生最终表。空白索引表头不能被当作
第16个科学指标；不得仅以文件存在或“15个物理字段”验收。

定义 R/L=原始/trim后配对数；T=all行数；U=all按QC键独特数；C/K/M/D=其cis及
close/middle/distal数；N=noBg行数；V/Cn/Kn/Mn/Dn=noBg按QC键对应计数。
下表顺序为原最终列顺序，正常公式仅适用分母有效且相应cis>0的情况。

| 原字段 | 单位/分子分母 |
| --- | --- |
| `total raw sequences` | R pairs，不是两端 read records 之和 |
| `after linker removing sequences` | L pairs |
| `mapping ratio` | Bowtie2 overall alignment rate百分数/100；不等同T/R |
| `total mapped PETs (mapq>=10)` | T PET行数；标题10随配置N代入 |
| `total mapped PETs redundancy` | 1−U/T |
| `total mapped PETs intra-chromosomal ratio` | C/U |
| `total mapped PETs close ratio (distance<=1kb)` | K/C |
| `total mapped PETs middle ratio (1kb<distance<=10kb)` | M/C |
| `total mapped PETs distal ratio (10kb<distance)` | D/C |
| `unique noBg PETs` | N PET行数 |
| `Yield` | N/R |
| `unique noBg mapped PETs intra-chromosomal ratio` | Cn/V |
| `unique noBg mapped PETs close ratio (distance<=1kb)` | Kn/Cn |
| `unique noBg mapped PETs middle ratio (1kb<distance<=10kb)` | Mn/Cn |
| `unique noBg mapped PETs distal ratio (10kb<distance)` | Dn/Cn |

映射率实际接线为 Bowtie2 输出→FLAG_A日志→`:330`读取倒数第二行百分数→`:472`除100，
依赖固定工具日志形状，H2 必须核实。不以另算的比率覆盖原值。
无cis时qc.py将cis置1，无unique时又将unique置1；TotalPETs仍可为0。
因此空/no-cis报告的cisRatio可能非零，不能当作真实计数发布。

## 4. 最小封装及失败边界

拟在H2资格验证启动器中：核身份和完整输入→准备全新、空的attempt输出目录→
运行固定原脚本→核顶层exit、最终表头/样本集合、gzip/BEDPE结构、参考端点与计数/比例
一致性→一次性完成标记；后续H3/H4再接平台生命周期和结果发布。原始数据只读，
旧attempt不复用，取消/失败保留私有诊断，不以部分产物发布成功。

依据：`tracPre2.py:363` 在建目录前打开日志；`:365–389,516–517` 缺工具/输入/索引
可return并退出0；`:177–179,236–241,261–262` 遇旧产物会跳过；`:457` 扫描日志。
[utils.py:103–123][utils] 忽略子命令退出码；Bowtie2的stat取得于固定 `tracPre2.py:245`（不是247），bamToBed位于310；二者均未被检查。
启动器应在运行前检查空目录，不能用日志创建后的“目录非空”消息判断污染。

后置条件不是所有子命令成功的充分证明。H2 要逐阶段注入失败；若发现原脚本退出0且
完整但错误的输出仍可通过，提交最小退出码收集方案待裁决，不擅改原科学代码。
不把原进程exit0、一个summary、canary或工具存在当作科学成功。

`services/defaults.py:97–157` 的默认runner仍明确装配Bulk准入配置；注册新metadata
不会自动授予执行能力。H3公共execution保持不可用，H4结果契约完成后才能开放，
不在通用平台/UI硬编码Hi-TrAC。共享闭包的改动另按实际diff评估身份和Gate。

## 5. H2之前已确认的政策（2026-09-24）

- **零PET/无cis：用户已确认政策A（2026-09-24）。** 任一样本all或noBg为空，或
  需要cis分母的QC集合无cis，则整个attempt不发布成功，不允许部分样本成功发布。
  保留原始输出与诊断，明确标识触发样本、集合及原因；不得修改上游计数。
  这是首版接入成功语义限制，不得据此描述实验没有科学意义。政策B（保留合法PET、
  受影响QC标为不可用）留作后续能力，本轮不实施NA映射或扩展结果/发布契约。
  H2对空输入、全部trim掉、全部unmapped、noBg为空、全trans分别保存原结果，再验证
  获批的封装策略。
- **中间产物：用户已确认BAM/裁剪FASTQ私有保留、不自动删除。** 首版平台不提供其下载，
  公开仅已有BEDPE/summary。成功、失败及取消的attempt均保留已有中间产物和诊断；
  原始FASTQ始终只读。后续清理或扩大公开范围需另行授权。原脚本内部清理边界见第3节。

多样本/预合lane已按roadmap边界设计；两项政策均已确认。依赖锁、
资源上限和参考兼容性须H2实测，不把H1源码阅读当部署准入。

## 6. 真实 tiny 基线设计（H2执行，H1未生成/运行）

使用seed `20260924` 设计两条100kb人工参考chrA/chrB；每端120bp，从目标参考正链取R1、负链
反向互补取R2，合法SEQ/QUAL等长、唯一配对ID。避开窗口和连接边界中意外的linker seed，
检查设计唯一性；仍须真实Bowtie2验证唯一比对/MAPQ≥10。保存生成器、FASTA/FASTQ及
真实`bowtie2-build`六文件索引的SHA、argv、版本和provenance。

下表为0-based半开目标端点，strand为+/−；a1/a2坐标相同但read ID不同。

| pair | 端1 | 端2 | linker | noBg预期 |
| --- | --- | --- | --- | --- |
| a1 | chrA:1000–1120 | chrA:1300–1420 | R1的120bp后追加完整linker | 保留，cis300 |
| a2 | 同a1 | 同a1 | 同a1 | 重复仅留一条 |
| b | chrA:3000–3120 | chrA:3300–3420 | 无 | 背景丢弃 |
| c | chrA:5000–5120 | chrA:6000–6120 | 无 | 保留，cis1000 |
| d | chrA:10000–10120 | chrA:15000–15120 | 无 | 保留，cis5000 |
| e | chrA:20000–20120 | chrA:35000–35120 | 无 | 保留，cis15000 |
| f | chrA:40000–40120 | chrB:1000–1120 | R1追加linker | 保留，trans |
| g | chrA:45000–45120 | chrB:5000–5120 | 无 | 背景丢弃 |

若8对都如设计比对，all有8行、QC独特7（cis5/trans2），noBg为a/c/d/e/f共5行
（cis4/trans1）。按第3节顺序手推15项：
`8, 8, 1, 8, 0.125, 5/7, 3/5, 1/5, 1/5, 5, 5/8, 4/5, 1/2, 1/4, 1/4`。
H2比较解压的端点/strand多重集合、重复次数、linker位置编码及最终表；计数精确，
比率建议绝对容差1e−12仅容纳浮点序列化，不能吞科学偏差。MAPQ和mapping ratio另与
真实比对日志逐项核对；不以手推值覆盖真实失败。gzip时间戳、行序及被改写的用户ID
不能作为科学等价要求，a1/a2保留哪条须按实际顺序解释。

H2增补矩阵：

| 用例 | 预期与验证方式 |
| --- | --- |
| 双样本、预合lane | 两个样本各得同一独立结果；跨lane重复仍在整样本去重，不跨样本去重 |
| linker/长度 | 单端、双端、RC、多命中、末尾9bp；trim至9bp丢弃、10bp保留trim计数，不能强求10bp通过比对 |
| 距离边界 | 999/1000/1001、10000/10001及奇数长度，分别核背景浮点中点、QC整数中点 |
| 低MAPQ/未比对 | 加等同重复参考窗口产生歧义，并亲核MAPQ<10；加入不比对pair，预期trim数增而PET不增；单mate失败单列，不凭序列随机性断言结果 |
| 输入/参考负例 | 缺mate、数/ID不匹配、坏gzip/FASTQ、空输入、缺/混配索引、身份漂移：拟议准入拒绝，原输入不改 |
| 零/无cis | 全部trim/未比对、仅b、仅f、仅g；保留原exit与原QC，按获批政策验收，不为通过而修上游数字 |
| 失败/旧结果 | 缺工具、子工具非零、部分输出、中间summary、旧BAM/trim/QC；全新attempt拒绝旧结果，故障注入与真实科学正例分开报告 |

H2必须实际运行原脚本与固定真实工具，记录合法FASTQ、参考、命令、退出码及全部原输出；
H1没有该证据。H3–H5才覆盖平台执行/取消、QC与产物发布、真实下载和桌面/移动产品链。

## 7. 上一轮 H2 实测与停止条件（2026-09-24，历史记录）

以下描述 `/tmp/helix-hitrac-h2-bx4wf6lj/` 的旧候选，不是新入口的当前状态。
独立复核支持两项反例；获批的后续实现及新的正式测试见第8节，旧失败和断言保留。

原脚本及固定包源码未修改。任务科学环境锁定 Python 3.11.14、Bowtie2
2.5.4（bioconda `he96a11b_6`）、samtools 1.23.1、bedtools 2.31.1、gzip 1.14，
以及 pandas 2.2.3、joblib 1.4.2、Biopython 1.85 等完整124项 conda 构建。
此前两个2.5.5构建的实际比对带额外警告行，触发原脚本日志解析错误；失败记录保留，
未过滤日志或修改解析器。开发环境单独通过原 canonical bootstrap，未放宽源码校验。

真实 `bowtie2-build` 产生六件 `.bt2`；原脚本 `-n 1 -p 2 -mapq 10` 在两个独立样本
上均得到 all=8、QC unique=7、noBg=5，实际MAPQ均42；最终16物理列（15指标加索引）
与第6节条件预期一致。双样本按各自预合lane的整对FASTQ处理，没有自动合lane。
距离、linker/RC、trim9/10和低MAPQ等原工具边界另有记录；末尾9bp未剪除的读段出现
带插入的实际比对，不能把所有设计端点一律当成实际端点。

**失败检测反例已触发约定停止条件：** 在任务私有故障替身中，让真实samtools生成的
过滤后BAM少一对完整记录并返回73。原脚本仍退出0，首样本all/noBg变成7/4；BAM、
BEDPE、最终表头、原日志比率与QC计数相互自洽，候选后置检查错误接受。
固定诊断断言得到1 passed / 1 failed，未放宽断言。这是故障模型的实证，不是声称
正常samtools必然产生这种失败，也不能由此宣称已经有可靠资格启动器。

另一个无故障替身的真实用例中，两对输入各仅一个mate比对；过滤BAM剩两个不同read ID
的单端记录，bedtools却输出末记录的同坐标cis PET。重跑同一bamToBed只能证明结果可重复，
不能证明配对有效。拟增加每条BEDPE由同名BAM read1/read2支持的完成核验，矛盾时拒绝
整个attempt；不新增比对过滤、不删除异常PET或改写上游结果。

最小待裁决方案是：在Hi-TrAC私有运行时用透明子工具入口，原样委派固定可执行文件，
逐调用记录开始/结束及真实退出码；顶层同时核完整调用记录和后置产物。非零、缺少结束
记录或记录不完整都拒绝完成；不改原科学脚本，不接公共平台。方案及反例在
`/tmp/helix-hitrac-h2-bx4wf6lj/proposal.md`，本轮没有实施该方案或写成功完成标记。

空输入、全trim原脚本退出1；全unmapped、noBg空、全trans等可退出0，原输出均保留。
政策A只在诊断后验中验证，尚未接入正式资格入口。缺mate/错ID/坏gzip、混配索引、
完整分阶段故障矩阵、取消与正式CI回归仍未完成，不能把这次基线称为H2已验收。
原始输入、工具/参考/产物摘要、命令及全部失败保存在上述任务目录；未进入H3。

[trac]: https://github.com/YaqiangCao/cLoops2/blob/de6cc732fa00b408551b9f4272933640c08447f1/scripts/tracPre2.py
[qc]: https://github.com/YaqiangCao/cLoops2/blob/de6cc732fa00b408551b9f4272933640c08447f1/cLoops2/qc.py
[ds]: https://github.com/YaqiangCao/cLoops2/blob/de6cc732fa00b408551b9f4272933640c08447f1/cLoops2/ds.py
[utils]: https://github.com/YaqiangCao/cLoops2/blob/de6cc732fa00b408551b9f4272933640c08447f1/cLoops2/utils.py
[setup]: https://github.com/YaqiangCao/cLoops2/blob/de6cc732fa00b408551b9f4272933640c08447f1/setup.py
[license]: https://github.com/YaqiangCao/cLoops2/blob/de6cc732fa00b408551b9f4272933640c08447f1/LICENSE

## 8. H2 正式私有资格入口（2026-09-24，已独立验收）

用户已批准的技术方案是**透明子工具退出码收集＋独立配对来源核验**，不是H1的
政策B；零结果政策A、整批拒绝、私有中间产物保留及原始FASTQ只读继续适用。
正式实现限定于 `adapters/hitrac_preprocess/`、两份专用脚本和自有工具锁，不接
registry/API/UI、默认worker或Bulk runtime。隔离交付来自已验收维护提交的许可文件，
不是完整clean Git checkout；源工作区的既有内容没有导入为实现。

- [准备入口](../../scripts/prepare_hitrac_bindings.py)验证固定runtime与独立可信的真实
  `bowtie2-build`记录，生成待审核binding；不会以“文件存在”代替索引构建。
  [工具锁及用法](../../config/hitrac_preprocess/README.md)固定124个conda artifact、
  原cLoops2源码/CLI和安装字节；runtime-v3、原脚本和包未改。重新安装到其他prefix的
  完整可重现性仍需实际验证；相同版本号不能代替字节匹配。
- [资格入口](../../scripts/qualify_hitrac_preprocess.py)使用原canonical源码provenance。
  流式核gzip/四行FASTQ、数量/ID/mate角色、输入/参考身份；创建0700全新attempt，
  源输入不写，staging用安全token和0400副本。完整六件小`.bt2`和参考先绑定SHA，
  运行后再次核工具、参考、实现及staging字节。旧attempt、混配/漂移、危险路径拒绝。
- 原脚本逐参数调用锁定绝对工具别名。每次独立开始/结束记录，含私有argv、cwd、
  样本/阶段/工具身份及真实退出或信号；不注入Bowtie2输出，不重试科学命令。
  计划每样本6次（align/view/sort/remove_sam/bedpe/compress），全批2次QC和1次
  remove_qc；双样本共15次。探测独立标记，不抵扣科学调用。非零、未结束、缺失、
  重复、损坏或未知调用均拒绝。原rm仅准许对应SAM及两个中间QC文件；不增加删除，
  不删除BAM/裁剪FASTQ。gzip移除未压缩BEDPE是原工具行为。
- 原产物核验流式读gzip，用私有SQLite核计数、QC键、noBg集合和15指标，按实际MAPQ
  核16物理列表头；不重写上游表。政策A定位样本、集合及原因；原脚本已在空/全trim
  分支自行失败时，报告原失败及调用缺口，不制造零结果成功。
- 独立来源核验从原BAM用samtools读取SAM，逐记录存私有SQLite；不再调用bamToBed。
  每条all PET需同一样本、精确qname的两条不同record，分别以read1/read2角色支持
  完整contig/CIGAR参考端点/strand；允许两端整体互换。搜索所有匹配alignment，
  不静默取首条；角色歧义不能作支持，但无关孤立/歧义记录不使整样本失败。
  不增加proper-pair、MAPQ、duplicate、secondary/supplementary过滤。
  SQLite关联不依赖samtools自然排序等于字典序，也不把整BAM或大qname组常驻内存。
  重复BEDPE逐行核验但不重写/去重，匹配记录可支持多个输出行；这是存在性证明，
  不声称还原唯一转换历史、输入完整性或所有科学正确性。
- 取消/超时清理仅该attempt的进程组，保留文件与未结束记录；原工具孙进程、wrapper
  死亡及SIGKILL都有隔离实测。完成前较长身份扫描后与原子rename前再查取消。
  私有诊断写入先完成，全部条件满足后才原子生成`complete.json`。没有此标记不成功；
  不允许部分样本成功。此Linux私有CLI不承诺清理主动setsid逃逸的任意恶意进程。

正式测试分层：`test/adapters/test_hitrac_*.py`是普通CI自动收集的快速层；
`test/hitrac_qualification/test_qualification_real.py`用既有`real_execution`标记，
显式提供经过核验的runtime/reference/tiny坐标后运行，缺前提失败、不skip。
普通CI不依赖本机`/tmp`路径。现有科学CI路径不包含此新目录，本轮没有把专用科学
测试接入远端job或声称远端通过；专用资格命令和结果在本轮报告中。

本轮快速层159项通过，显式真实资格层24项通过（12项无替身科学、10项故障注入、
2项取消/超时），均无失败或skip；重复运行不累计。真实新参考、原脚本双样本all8/noBg5
及15指标、MAPQ10/17动态表头、linker/RC、
长度/距离/低MAPQ和真实singlemate反例均有正式测试。逐工具返回73（含原rm）、
少一对BAM后退出73、取消/超时另作为注入/生命周期证据。旧accepted=true反例及
冻结红灯保持原样，不把新入口通过写成旧候选已通过。
本轮证据 `/tmp/helix-hitrac-h2-impl-gzfw7n63/`；未启动平台服务、Bulk Gate或H3。


## 9. H3 adapter 配置与执行接线（2026-09-24，阶段记录）

H2 已独立复跑159项快速及24项资格测试通过，另有threads=3/8真实双样本CLI对照；
见 `/tmp/helix-hitrac-h2-independent-0d3nw0cc/report.md`。第7节为先前裁决点，
第8节的原始证据继续保留。H1政策A和中间产物保留政策均未改变。

### 作者输入与服务器权限

[authoring.py](../../src/encode_pipeline/adapters/hitrac_preprocess/authoring.py)与
[validation.py](../../src/encode_pipeline/adapters/hitrac_preprocess/validation.py)使用现有
adapter-owned schema/Issue：config仅空对象；samples为有序多行
`sample_id/fastq_1/fastq_2`；options仅threads（2～8，默认2）、mapq（0～255，默认10）。
每行一对外部gzip FASTQ，多lane预先合并；不开放managed输入、任意argv、运行时、
参考摘要或timeout字段。样本ID须唯一；JSON schema的uniqueItems只能判整行相同，
不能单独代替ID唯一性检查。校验Issue使用固定定位/提示，不反射非法字段或路径。
作者校验不读取FASTQ、不执行科学工具；实际配对、gzip及内容身份在服务端绑定时核验。

用户使用既有`reference_profile_revision_id`选择参考。私有参考配置的workflow键为
`hitrac-preprocess`，payload为`schema_version: hitrac-reference-profile-v1`及
`binding/sha256`（服务器已审核的H2参考绑定文件和摘要）。既有ReferenceProfile服务
核对修订与当前绑定。部署方通过`HELIXWEAVE_HITRAC_RUNTIME_BINDING`和
`HELIXWEAVE_HITRAC_RUNTIME_SHA256`提供经过审核的运行时；这是服务器坐标，
不是用户执行开关。绑定失败不授予可执行入口。

[adapter.py](../../src/encode_pipeline/adapters/hitrac_preprocess/adapter.py)无论未准入、
已准入或已绑定参考，公共availability始终为`not_configured`，沿用既有固定原因码。
metadata说明H4产物/QC发布尚未实现；现有服务隐藏不可用的执行能力、拒绝公共
executable snapshot和提交。默认registry可提供元数据和作者校验，不能公开跑科学任务。
没有artifact/QC capability。基础WorkflowAdapter协议强制的extract_artifacts方法仅返回
unsupported拒绝；未实现extractor、产物登记或发布。无环境变量或用户开关可绕过此门禁。

### 身份、计划与私有执行

[execution.py](../../src/encode_pipeline/adapters/hitrac_preprocess/execution.py)绑定规范化
输入/参数、原输入字节SHA、原行序与`s000001…`token映射、参考、固定runtime、
当前任务解释器及实现闭包。原始路径不传给上游shell；H2在新attempt内复制并核对输入。
此绑定将由已有validated snapshot服务消费；H3测试只在局部availability替身中验证
snapshot/重复创建/陈旧拒绝，不把这种注入写为公共可用。

WorkspacePlan只创建私有请求，不预建科学attempt；物化后command调用
[run_hitrac_preprocess.py](../../scripts/run_hitrac_preprocess.py)（正常安装的任务Python，
`-I -S -B`），经canonical provenance到原资格入口。执行前重核请求摘要、绑定、cwd和
全部输入；staging后核对已冻结输入及实现，完成前保留H2全部后验。旧attempt拒绝，
不删除或复用旧科学输出。完成标记只代表隔离资格完成，不能触发H3公开成功发布。
仅经过具体HiTrAC部署对象身份核验的解释器进入默认runner allowlist；不信任任意adapter
返回的executable，也不要求此adapter配置Docker。ENCODE/Bulk的原分支保留。

Hi-TrAC实现摘要包含包内全部Python、工具锁、入口、bootstrap及实际共享计划/运行/
worker控制文件。新增模块必然产生新身份。外层解释器在准入和入口重核，未声称运行
过程中替换解释器文件也一定在结束时检测。内容漂移拒绝不等于对同字节inode替换作承诺。

### 取消边界与证据等级

H2科学进程会新建session。原worker仅kill horse group存在嵌套进程存活及扫描后fork窗口，
同断言真实小进程红灯保存在本轮证据。`DurableWorker`对默认SIGKILL冻结匹配PID/starttime
的自有树，稳定重扫后终止并收割后代；horse等待结果留给RQ。wait边界等待清理完成，
失败不能形成取消确认；非SIGKILL原分支保留。不改SQLite状态、job identity或事件代码。
这是Linux受信科学调用树控制，不是任意恶意程序的安全沙箱。ProcessRunner自身已有
嵌套清理，本轮未改它。对应设计和故障边界见本轮`execution-review.md`。

快速层新增schema、reference/snapshot/plan/allowlist和公共关闭测试；真实资格层新增
`test/hitrac_qualification/test_adapter_real.py`与`test_adapter_control_real.py`，复用H2
显式坐标。前者实际经过adapter→planner→materializer→command→默认ProcessRunner；
后者区分真实fork horse控制与窄工具故障注入。没有Redis/RQ部署、API科学执行、H4公开
artifact/QC、下载或浏览器产品链。普通CI快速层不依赖本机科学工具；显式real_execution
层缺前提失败、不skip，尚无远端Hi-TrAC专用job验收。

本轮实施和实际日志：`/tmp/helix-hitrac-h3-0uj0yyu1/`。选择性隔离副本不冒称clean Git。
H3初次交付同步了Bulk111项manifest/qualification；独立复验随后确认共享runner的
Hi-TrAC准入helper仍在该闭包之外，原身份不能完整约束这项执行权威。返修将实际控制链
纳入正式闭包，并在受控部署入口固定其使用的tools.lock.json摘要；不把未被调用的科学
处理模块全部加入Bulk。Hi-TrAC自己的全部实现身份继续独立核对。

返修同时区分自有进程的存在、已消失/身份改变与无法确认：身份读取或子树发现失败
不能当作空树；完整冻结前根或已登记父进程消失/成为僵尸，拒绝清理确认。发现只沿
已确认自有节点的task/children进行，逐层冻结，不扫描/杀死任意新收养的子进程。
正常路径保留原RQ horse等待结果和非SIGKILL行为；故障路径的拒绝确认不等于自动
找回并清理所有失根后代。注入测试遗留进程由测试按已知PID/starttime回收，这不是
产品自动恢复能力；真实kill/wait也不能代替完整Redis/RQ持久取消链验收。

按用户本轮决定，完整Protected Bulk Gate在H5完成后的集成收尾统一执行；H3/H4仍
完成各自定向验证与必要身份同步，不因完整Gate尚未执行而单独阻塞交付。生成资格记录
不等于运行Gate。最终Gate只覆盖最终精确版本，不追认本H3或历史身份。完整ENCODE
身份因保护profile边界未捕获；既有维护commit的Gate及PR-6在线staging分别记录。
返修587项快速/相关回归及32项分层真实资格通过，零失败/skip；Bulk正式闭包
111→118，Hi-TrAC自身实现摘要同步变化。精确摘要、同断言红绿及原始输出见
`/tmp/helix-hitrac-h3-fix-jgvb1rzl/`。返修完成，待独立验收，不进入H4。

### H3 早期并发返修（2026-09-24，本轮验证通过，待独立验收）

上一轮顺序失败和执行身份修复获复验支持，早期查询的并发确认窗口另行补齐。
`DurableWorker.kill_horse` 的锁覆盖初始准备、group查询及重试、所有失败出口和
最终完成登记；原 `wait_for_horse` 的wait4保持锁外。匹配停止job、horse PID/starttime
的完成证明才允许停止确认，RQ已写stop标记但尚未进入kill也不能提前确认。
未完成冻结时失根仍拒绝确认，不承诺自动找回或杀死归属未知后代。

同字节9项并发红灯3失败/6通过，修复后与原worker回归共59通过，相关回归239通过，
分层真实资格9通过；证据 `/tmp/helix-hitrac-h3-h4-l8ug4hfs/stage-a/`。Bulk与Hi-TrAC
分别同步身份，完整Gate留H5后。用户授权冻结本阶段后直接继续H4，公共执行须待结果
契约实现并验证后才可开放；当前仍关闭。真实kill/wait不等于完整Redis/RQ持久取消链。

## 10. H4 结果、QC及通用界面（2026-09-24，待独立验收）

### 准入、原始来源与公开范围

[results.py](../../src/encode_pipeline/adapters/hitrac_preprocess/results.py)的Results子类
提供实际artifact/QC能力，部署组合仅在服务器runtime准入成功后提供执行能力。
未配置、绑定漂移、输入/参考/实现身份不匹配仍拒绝；H3基类仍保持关闭。
未加入用户绕过开关。真实结果测试经过当前adapter、原planner/materializer、
command及默认ProcessRunner；没有用availability替身开放产品。

提取重新核对本次私有计划、workspace/attempt、完整调用记录、退出状态、输入及
参考字节、当前实现和完成标记。随后复用H2后验及独立BAM配对来源核验。仅有
`complete.json`、summary或exit0不构成发布依据。内部token必须与冻结的有序样本
映射精确一致；缺少、重复、额外、错配样本或任一样本失败，整批拒绝发布。
H1政策A、原始输入只读以及所有已有中间产物和诊断私有保留均不变。

公开登记仅为每样本原`*_all.bedpe.gz`、原`*_unique.bedpe.gz`（上游noBg集合）
及原`tracPre_summary.txt`。分别使用`hitrac_bedpe_all`、`hitrac_bedpe_no_bg`、
`hitrac_summary`类型；只有summary作为QC来源。H4核验gzip、十列BEDPE、参考端点、
实际集合、计数、来源及输出SHA，再按原字节投影到现有`results/`边界。
不覆盖不同的既有投影，不删除原产物。BAM、裁剪FASTQ、请求、日志、调用记录和
其他私有文件不登记、列出或提供下载。正式下载继续使用平台原身份/路径检查，
没有新下载旁路。

### 原15指标与平台投影

原summary严格为15个指标加一个样本索引，共16物理TSV列；列名按实际MAPQ参数
生成，不硬编码10。原文件字节保留；内部token仅在平台metadata及QC坐标映射回
原用户样本名。计数/比例、分母、有限值及原输出一致性在投影前验证；非法或部分
结果不改写成零、NA或成功，也不新增科学好坏判级。

| 平台metric_key | 原统计含义 | 单位 |
| --- | --- | --- |
| raw_pairs / trimmed_pairs | total raw sequences / after linker removing sequences | count（read pairs） |
| mapping_ratio | 原mapping ratio | fraction |
| all_pets | MAPQ过滤后的all PET数 | count（PET） |
| all_redundancy | all集合的原redundancy | fraction |
| all_cis_ratio | all集合QC unique PET的cis比例 | fraction |
| all_close_ratio / all_middle_ratio / all_distal_ratio | all集合QC unique cis PET的原距离分组比例 | fraction |
| no_bg_pets | 原unique/noBg集合PET数 | count（PET） |
| yield | noBg PET / raw read pairs | ratio |
| no_bg_cis_ratio | noBg集合的cis比例 | fraction |
| no_bg_close_ratio / no_bg_middle_ratio / no_bg_distal_ratio | noBg cis PET的原距离分组比例 | fraction |

all集合QC中的unique统计不等于`unique/noBg`文件的计数。原算法和距离边界沿用
第3节定义，不重新计算新科学指标。平台既有Decimal契约最多12位小数：比例投影
沿用既有Bulk的ROUND_HALF_EVEN方式，最多舍入至12位；误差至多0.5e-12，计数不变。
原summary、完整科学核验值和BEDPE均不改变，下载保留完整原值。

### 原子发布与通用样本名

用户已批准内部原子bundle方案。Hi-TrAC在提取所有artifact、解析全部QC并完成
公共候选校验后，通过原RunService/repository提交同一bundle：artifact、QC、
二者generation、attempt、事件及append-only publication处于同一SQLite事务；
内存repository有对应原子语义。异常回滚不得留下部分可见的新结果；相同attempt
重复调用只观察同一批结果，陈旧或不同bundle拒绝。原ENCODE/Bulk普通分支保留。
这是实际Hi-TrAC消费者使用的内部协议，未增加公共字段或表结构。

用户另行批准只扩展QC `sample_id`：原ASCII字母、数字、点、下划线、连字符之外
仅增加ASCII空格U+0020，长度1～255。空串、纯空格、控制字符、其他空白、路径分隔符、
冒号及`.`/`..`拒绝。候选、来源metadata、持久化与公共响应同步；不strip、不折叠
空格、不替换为内部token。`experiment_id`、`assay`、数据注册表及其他ID契约不变。
既有metric ID算法不变，`ab`、`a b`、`a  b`为不同坐标，尾空格逐字节往返。

通用QC/产物页面使用[SampleIdentity](../../frontend/src/components/SampleIdentity.tsx)：
含空格ID以带引号、保留空白的等宽文本展示，并提供原样复制；不依赖workflow_id。
连续空格和末尾引号让尾空格可辨识。正式OpenAPI导出与原JSON逐字相同，因此未修改
生成客户端；字段的运行时字符校验改变在此明确记录。正式构建及资产打包已同步。

### 本轮验证与尚未覆盖

当前H4身份的真实双样本输出经正式提取、SQLite bundle、原认证HTTP列表/QC和
下载链验证；每样本all=8/noBg=5，原15指标和PET多重集合保持。MAPQ10/30独立结果
对照及MAPQ17发布链均验证动态表头。另覆盖marker/调用记录/输入/BAM/summary/gzip
漂移、真实singlemate和多样本空集合整批拒绝；退出73及外层取消/超时有分层证据。
样本ID专项57项覆盖两repository和HTTP；共享bundle回滚/幂等及原分支回归保留。

桌面1440×900、移动390×844页面直接消费本轮真实输出的数据库，每视口30项QC、
5项artifact；10次实际下载的SHA/大小与登记一致，连续/尾空格展示和剪贴板原值通过。
浏览器使用临时loopback测试应用，服务和浏览器已停止。这不是完整Redis/RQ科学产品
链；该链留H5，远端Hi-TrAC验证未执行。完整ENCODE摘要未绕过保护profile捕获。

各阶段独立diff、身份恢复材料、原始命令及所有失败记录位于
`/tmp/helix-hitrac-h3-h4-l8ug4hfs/`。H4 Bulk实际闭包119项，Hi-TrAC自身63项；正式
生成manifest/qualification只代表身份同步，不是Gate通过。完整Gate按用户决定留H5后，
仅覆盖最终精确版本，不追认中间/历史身份或PR-6在线staging链。

## 11. H3/H4 F1/F2 限定返修（2026-09-24，返修完成，待独立验收）

独立复验发现晚到停止标记与原子发布失败通知两处边界；本轮只修这两项，
未更改科学实现、H1政策、公共模型、表结构或发布算法。

[DurableWorker](../../src/encode_pipeline/workers/timeouts.py)让停止登记受清理锁约束，
并在原RQ monitor整个调用上下文中的停止标记消费处要求匹配job、horse PID和
starttime的完整清理证明。原wait4保持锁外；原停止回调分支以及随后原失败处理器
再次消费标记都受同一证明限制。诊断读取不等于停止确认。未知、失根和清理失败
继续拒绝确认；不承诺自动找回已失根后代，也不把真实进程测试当作Redis/SQLite
持久取消链验收。

[worker结果索引](../../src/encode_pipeline/workers/jobs.py)通过既有
AtomicResultPublishingAdapter内部协议识别原子发布分支，不按workflow_id特判。
只有本次精确artifact attempt提取成功且artifact/QC完整bundle匹配时才进入成功
通知点；准备失败、提交异常或不完整提交都不发送成功通知，也不以旧QC代替本次结果。
旧artifact/QC/generation及publication保持；科学运行SUCCEEDED语义不变。
非opt-in ENCODE/Bulk原通知政策、邮件发送失败处理及去重均不改变。

新增正式回归分别见
[test_late_stop_monitor.py](../../test/workers/test_late_stop_monitor.py)与
[test_atomic_bundle_notifications.py](../../test/workers/test_atomic_bundle_notifications.py)。
证据与身份恢复材料位于`/tmp/helix-h3-h4-fixes-iwa4yxgq/`。完整Gate仍留H5后，
只证明届时的最终精确身份；当前定向通过与qualification生成不等于Gate通过。
非阻塞的只读QC followup attempt观察未纳入本轮。

本次最终受影响集合476通过，F2通知专项18及原通知/worker35通过；当前新身份
真实资格9通过，均零skip。真实输出再次核验all=8/noBg=5、15指标和PET多重集合，
MAPQ10/30及发布链17保持；正式提取/SQLite/认证HTTP下载字节一致。该结果不包括
完整Redis/RQ持久链、浏览器或完整Gate，本轮未修改前端。

## 12. H5 全链验收（2026-09-24，本轮交付，待独立验收）

H3/H4 F1/F2 返修已获独立验收（`/tmp/helix-f1-f2-independent-Z3QODmwn/`）。
本轮在 `/tmp/helix-h5-XTvzkxCP/` 的隔离副本完成真实平台全链验收；新增正式
资格入口 `test/hitrac_qualification/test_platform_real.py`（真实 Redis/RQ 产品链）
与 `test_publication_real.py`（F2 发布失败定向集成），均为 `real_execution`
显式层，缺坐标明确失败、不静默 skip。部署/准入/升级与上游耦合记录在
[hitrac-preprocess-operations.md](hitrac-preprocess-operations.md)。

真实链：认证 HTTP API（登录+CSRF）→ 配置校验与 validated snapshot →
create/preflight/start → 任务专用真实 Redis 7.0.15 与 RQ 2.10 → 原
`encode_pipeline.workers.cli` worker 与默认 runtime 装配 → 固定 tracPre2 与
真实科学工具 → SQLite 原子 artifact/QC 发布 → 原 API 与实际下载。不使用进程内
ASGI、fakeredis、直接 runner 调用或预填数据库。双样本（含连续/尾部空格 ID）
各 all=8、noBg=5，15 指标与规范化 PET 多重集合符合 H2 基线；5 产物、30 QC、
5 publication 与 generation 一致；下载字节与登记摘要一致；无权/私有访问拒绝；
关闭重开后终态与结果可读；相同输入新 run 使用新 attempt，同 snapshot 重复创建
幂等回放。取消经真实 API/RQ stop/monitor/回调持久化，记录 job、horse
PID/starttime 与事件顺序并重开复核；超时按既有契约持久化
PROCESS_RUNNER_TIMEOUT，嵌套科学进程清理、无成功发布、无关 sentinel 不受影响。

桌面 1440×900 与移动 390×844 浏览器消费同一次真实 run（第三次打开同一数据库），
每视口 30 QC、5 产物下载 SHA 一致，空格样本名带引号展示与原样复制通过。
F2 发布失败的首次/旧 bundle 场景经原服务、SQLite repository 与真实 SMTP 客户端
到 loopback 捕获验证；该层不经过真实 worker，真实 worker 内 bundle 提交失败的
确定性注入接缝不存在，缺口在 operations 文档与本轮报告明确标注。

完整 Protected Bulk Gate 按用户决定留本轮之后的集成收尾，只覆盖最终精确身份；
本轮不追认历史/中间身份，G0、PR-6 在线 staging 全链与完整 ENCODE 摘要保护
边界继续保留。证据、命令与原始日志：`/tmp/helix-h5-XTvzkxCP/`。

## 13. H5 取消确认返修（2026-09-25，本轮修复，待独立验收）

H5 之后的完整 Bulk Gate-2 取消用例未达终态：run 永久停在 `running`
（`ended_at`/`cancellation_acknowledged_at` 均空）。根因不在 Redis/QoS，而在
[DurableWorker](../../src/encode_pipeline/workers/timeouts.py) 的冻结证明本身。

H5-era 冻结窗口内，任意一个已登记下行成员**正常退出**即触发致命分支。其父进程刚被
`SIGSTOP`，无法 `wait()` 回收，该成员以僵尸形态留在
`/proc/<ppid>/task/<tid>/children`；`_require_frozen_member` 把「不存在或已是僵尸」
一律判为失根致命 → `kill_horse` 抛错 → `_nested_cleanup_completed` 永不写入 →
`wait_for_horse` 抛错 → RQ 裸 `except` 吞掉并以 0 退出 → 停止回调不执行。Nextflow
提交突发（大量 `.command.run` 同时存在）下这是常态，不是异常；F1 契约只要求**根**的
死亡保持致命，并未要求对正常结束的下行成员中止——该行为是实现副作用。

根因由两条独立证据链确证：(1) 真实 Gate-2 现场 worker 私有 stderr 的两条 traceback
（停止线程 `command.py:141 → timeouts.py:434 → :183 → :80`；主线程
`base.py:636 → worker_classes.py:156 → :346 → :88 → :372`），`worker-exit.json`
为 `{"returncode":0}`；(2) 不依赖 Redis/容器/Nextflow 的独立微复现
`probe-cancel-freeze-abort-4`，触发点为僵尸成员、其父是树内成员、`is_root=false`。
RQ/Redis 层异常与回调硬超时已排除。

按用户裁决做**有界放宽 + 内核正向证明**（不改变契约本质、不延长窗口、不直接标
CANCELLED、不删除失根检查）：

- **根严格不变**：根必须存活且冻结；根的消失/僵尸仍致命。
- 已登记成员的正常退出（`_member_visibility == "gone"`）改为「已知归属的退出事件」：
  PID+starttime 记入 `exited`，从冻结证明中剔除。**身份改变**（starttime 不符）仍致命；
  身份读取被拒或格式异常仍致命，**不当作退出**。
- **登记时序**：`registered`（根 + 全部 tracked + 全部 exited）在第一个 kill 信号之前
  成型；`_record_owned_cleanup` 在写完成标记**之前**先落日志/报告（证据先于确认）。
- **正向证明①（同组完备性）**：`killpg(horse_pid)` 之后，所有 `pgrp == horse_pid` 的
  存活进程必须为空。进程组归属由内核维护、不随重挂丢失，覆盖已退出成员留下的**同组**
  孤儿。
- **正向证明②（孤儿归属检查，诊断）**：`_uncovered_processes` 报告既不在登记集合、
  又落在 horse 进程组、或作为 worker 新收养子进程出现的存活进程；只记诊断事件，
  不静默，也不作为确认依据。

**残余风险（显性记录）**：在本次清理作用域收养树之前就已退出、且已 `setsid` 的成员，
其独立 session 子树在内核中已无进程表项可枚举，既不在 horse 组内、也不一定在 worker
收养表内，因此不在证明覆盖范围内。取消路径对该形态保留诊断钩子（`uncovered` 记录 +
worker 警告日志），但这是「本轮证据中未出现」，不是「证明不可能出现」。关闭该残余需要
worker 全程 subreaper + 全作业期孤儿归属（更大的契约改造），不在本轮范围。

回归见 [test_relaxed_owned_cleanup.py](../../test/workers/test_relaxed_owned_cleanup.py)：
动态覆盖「成员在冻结窗口内正常退出仍确认清理」，并含否定用例「根丢失仍致命」「身份
改变仍致命」，以及登记时序（第一个 `killpg` 之前无新的归属发现、每个 SIGKILL 目标都
已登记）与 `_member_visibility` 四态分类。原 F1 6 项、F2 18 项、timeouts、nested/owned
horse、atomic bundle 共 **83 项**保持通过，`test/workers` 全层 **304 通过**。
`workers/timeouts.py` 同时属于 Bulk 与 Hi-TrAC 两个执行闭包，Bulk 走正式生成器
（`scripts/generate_bulk_rnaseq_execution_manifest.py`）、Hi-TrAC 走产品函数
`implementation_identity()` 重算：Bulk `965a9e78…`→`b4dd6488…`，Hi-TrAC
`fe775ae8…`→`88e81e8b…`，路径集不变（119/63），persistence contract
`ecf0e206…` 不变。冻结后同一工作区内，独立探针以同一棵真实进程树记录登记时序与
残余形态：`probes/cancel_freeze_abort_probe.py`（A 分支修复前后由
「Workflow tree was lost before freezing.」变为确认清理；C 场景的独立 session 孤儿
在清理后仍存活）与 `probes/registration_order_probe.py`（首个 killpg 之后无归属发现、
每个 SIGKILL 目标都在 `registered` 内、组内无存活进程、两条警告均已落日志）。

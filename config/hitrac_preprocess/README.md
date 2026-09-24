# Hi-TrAC H2 私有科学运行时

这不是公共 adapter 注册或 Bulk runtime。`conda-explicit.lock` 固定 Linux x86-64
的 124 个 conda artifact；`tools.lock.json` 同时记录 URL、SHA256、许可证和
已经真实验证的安装字节摘要。不得用相同版本号的另一构建或 PyPI cLoops2 替代。

准备时，先按 explicit 文件列出的准确 URL 获取 artifact，逐件核对
`tools.lock.json` 的 SHA256，再用任务局部 micromamba 的 `--offline` explicit
入口创建独立 prefix。安装原 cLoops2 Git commit
`de6cc732fa00b408551b9f4272933640c08447f1` 的正常 wheel，不修改源码；用固定
Python 执行原 `scripts/tracPre2.py`。wheel 内的 cLoops2 Python 文件必须逐件
匹配锁中的 `upstream.package_sources`。保留 BSD-3-Clause 及各 conda 包许可证。

运行时绑定文件仅声明 `schema_version: hitrac-runtime-binding-v1`、`prefix`、
`script` 和本锁的 `lock_sha256`。它不能自行选择其他 pin。资格入口重新核对
conda 包集合、version/build/artifact SHA、实际安装文件集合与字节摘要、原脚本、
cLoops2 包源码及工具 CLI 入口，包括 `bamToBed` 的原别名和 `cLoops2` 的原解释器。
`rm` 是明确固定的本机 `/usr/bin/rm` 字节，不继承任意 PATH。

安装摘要仅将实际 prefix 字节替换为 `${PREFIX}`，以容纳正常安装的解释器和
脚本 prefix；不修改运行文件。`.pyc` 与 `__pycache__` 不进入安装摘要。完整
重新安装/不同 prefix 长度的可重现性必须实际核验，不能仅凭这种归一化宣布
所有 relocation 都兼容；不匹配时拒绝，禁止现场刷新锁来制造通过。

参考使用单独经审核的构建记录：FASTA、完整六件 `.bt2`、contig 长度和逐件 SHA。
调用者须另行提供该记录的预期 SHA，不能从新提供的记录自行计算并据此宣称
受信。`.bt2l` 尚未验证。资格入口检验记录与各文件，staging 时再次核对每件字节。
快速测试中的合成索引字节只验证准入/漂移检测，不证明 Bowtie2 索引科学可用；
真实资格层使用真实 `bowtie2-build` 的六件输出与保存的成功构建记录。

输入暂采用完整四行记录的 gzip FASTQ（允许 concatenated gzip members），
逐对核对数量、ID、`/1`/`/2` 及存在的 Illumina mate 角色、SEQ/QUAL、字符及 gzip
校验。空数据可以完成格式准入，其后受 H1 政策 A 和原脚本执行结果约束。
输入文件可包含空格；原脚本只接触安全内部 token 的私有副本，不接触用户源路径。

已有运行时和真实索引的绑定准备入口（替换路径；不读取默认配置）：

```bash
TASK_PYTHON=/absolute/task/dev-env/bin/python
"$TASK_PYTHON" -I -S scripts/checkout_bootstrap.py --repository-root . verify-checkout
"$TASK_PYTHON" -I -S scripts/prepare_hitrac_bindings.py \
  --runtime-prefix /absolute/pinned/science-prefix \
  --tracpre2 /absolute/fixed-source/scripts/tracPre2.py \
  --reference-fasta /absolute/reference.fa \
  --index-prefix /absolute/index/genome \
  --index-build-record /absolute/bowtie2-build-command.json \
  --index-build-record-sha256 REVIEWED_BUILD_RECORD_SHA256 \
  --output-directory /absolute/new-private-bindings
```

构建记录须包含真实命令的 `argv`、`exit: 0`、`status: "completed"`；argv 首项为
该固定运行时的 `bin/bowtie2-build`，末两项为对应 FASTA 和索引 prefix。
入口拒绝复用已有 output 目录。它核对现有安装和参考字节，生成私有 binding 与摘要；
**它本身不执行索引构建，也不自动批准新参考**。输出的 reference binding SHA
须作为独立审阅对象，之后明确传入资格入口。

原源码归档SHA256为`e5f14c7f6723eb02293a2377b1f060ed9fd8f6f852a2d10b5d27f7aff79acdc5`；
本次正常构建、复用的cLoops2 wheel SHA256为
`d71c0175eb8bfa8e5847180395a5e60a437b8021486dafd715ed0194ed4e00ec`。
安装用`python -m pip install --no-index --no-deps /path/to/reviewed.whl`；
此为已有构建物摘要，不声称含时间戳的wheel重复构建后必然同SHA。
源码/包实际执行字节仍须满足本锁，不用PyPI同名版本替代。
固定上游与[BSD-3-Clause许可](https://github.com/YaqiangCao/cLoops2/blob/de6cc732fa00b408551b9f4272933640c08447f1/LICENSE)
保存在原源码归档；各conda包来源、artifactSHA及license字段在锁内，分发时保留各包许可。

私有资格入口请求示例（JSON文件权限0600，示例路径必须换成已核验的实际文件）：

```json
{
  "runtime_binding": "/private/runtime-binding.json",
  "reference_binding": "/private/reference-binding.json",
  "reference_sha256": "独立审阅后的64位参考binding摘要",
  "samples": [{"id": "sample 1", "r1": "/inputs/lane-merged_R1.fastq.gz", "r2": "/inputs/lane-merged_R2.fastq.gz"}],
  "threads": 2,
  "mapq": 10,
  "timeout": 300
}
```

```bash
"$TASK_PYTHON" -I -S -B scripts/qualify_hitrac_preprocess.py \
  --request /private/request.json --attempt /private/new-attempt
```

attempt须不存在，父目录已存在且canonical路径无空白或shell元字符；原始FASTQ路径允许
空白。成功条件是原子`complete.json`，不是stdout或summary出现。该文件仅供H2资格审核，
不会登记平台产物。公开范围仍仅原BEDPE和summary；0700目录中的argv、异常、日志、BAM、
trim FASTQ、SAM检查副本和SQLite核验文件都不作为公开artifact，成功/失败/取消均保留。
超时参数约束科学子进程（和单次BAM读取），不是包括身份散列的硬端到端上限。
Linux进程组负责原调用树；这不是任意恶意程序的沙箱，未证明主动setsid逃逸的清理。

测试分层与明确坐标（无本机硬编码路径）：

```bash
"$TASK_PYTHON" -I -S -B scripts/checkout_bootstrap.py --repository-root . pytest \
  test/adapters/test_hitrac_admission.py test/adapters/test_hitrac_outputs.py \
  test/adapters/test_hitrac_pairs.py test/adapters/test_hitrac_calls.py \
  test/adapters/test_hitrac_process.py test/adapters/test_hitrac_qualification.py \
  --basetemp /private/new-fast-basetemp

# 先运行 tiny_inputs.py，真实 bowtie2-build，再用上面的prepare入口生成并审核binding。
export HELIXWEAVE_HITRAC_RUNTIME_BINDING=/private/runtime-binding.json
export HELIXWEAVE_HITRAC_REFERENCE_BINDING=/private/reference-binding.json
export HELIXWEAVE_HITRAC_REFERENCE_SHA256=REVIEWED_REFERENCE_BINDING_SHA256
export HELIXWEAVE_HITRAC_TINY_INPUTS=/private/generated-tiny
"$TASK_PYTHON" -I -S -B scripts/checkout_bootstrap.py --repository-root . pytest \
  -m real_execution test/hitrac_qualification/test_qualification_real.py \
  --basetemp /private/new-real-basetemp --junitxml /private/real-junit.xml
```

快速层随现有CI分片自动收集。真实层必须显式选择marker和目录；不选marker时仓库默认
排除真实工具测试，不能把deselected算通过。环境缺失时失败，不添加skip。本轮本机
执行结果不等于远端资格job已接线或已通过。故障测试仅在经过真实runtime准入后替换
一项工具依赖，仍调用原科学工具及原脚本；它们与无替身科学基线分开计数。

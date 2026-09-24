# Hi-TrAC 预处理 adapter：部署、准入、升级与上游耦合

状态：与 H5 全链验收同时交付。本文记录实际部署坐标、运行时/参考准入、升级
路径和上游耦合点；科学契约与政策见
[hitrac-preprocess-contract.md](hitrac-preprocess-contract.md)，维护顺序见
[workflow-platform-agent-roadmap.md](workflow-platform-agent-roadmap.md)。
所有身份值来自固定锁文件与本轮实测记录，不是示例。

## 1. 固定身份

| 对象 | 身份 |
| --- | --- |
| workflow ID | `hitrac-preprocess`（独立 ID，不改 ENCODE/Bulk 注册） |
| 上游脚本 | cLoops2 Git commit `de6cc732fa00b408551b9f4272933640c08447f1`（tag v0.0.5）的 `scripts/tracPre2.py`，SHA256 `c3c4ef4e6287fa4ade97a6f5980d20c81ea345b8e8e13c88ae3b7c4bd3a67aec`；PyPI 同名 0.0.5 不是等价身份 |
| 运行时锁 | `config/hitrac_preprocess/tools.lock.json`，SHA256 `392e4a4ff93f2b14b380f16572d74a77dbca699112800bcafe50992d692090da`；124 个 conda 构建逐项锁定（Python 3.11.14、Bowtie2 2.5.4 `he96a11b_6`、samtools 1.23.1、bedtools 2.31.1、gzip 1.14、pandas 2.2.3、joblib 1.4.2、Biopython 1.85 等） |
| 工具入口 | `python`、`bowtie2`、`bowtie2-build`、`bowtie2-inspect`、`samtools`、`bamToBed`、`bedtools`、`gzip`、`cLoops2`、`rm`，逐入口 SHA256 锁定；Bowtie2 仅验证过 2.5.4 固定构建（两个 2.5.5 构建的日志形状不兼容原解析，见契约第 7 节） |
| 实现闭包 | `adapters/hitrac_preprocess/qualification.py` 的 `implementation_identity()` 聚合 adapter 全部模块、工具锁、三个入口脚本与共享平台/持久化/worker 文件；准入、执行前后与发布前重复核验 |

## 2. 部署与准入

平台侧不新增服务：复用现有 SQLite、Redis/RQ worker 与 API。Hi-TrAC 的执行
能力由两个服务器坐标门控，缺任一项或核验失败时 availability 保持
`not_configured`，创建/提交/执行全部 fail-closed：

- `HELIXWEAVE_HITRAC_RUNTIME_BINDING`：经审核的 runtime binding JSON 路径
  （`hitrac-runtime-binding-v1`），内容绑定运行时 prefix、固定脚本与锁摘要。
- `HELIXWEAVE_HITRAC_RUNTIME_SHA256`：该 binding 文件的 SHA256。
  装配入口 `adapters/hitrac_preprocess/deployment.py` 的
  `load_default_hitrac_adapter`；两处不一致即拒绝授予能力，没有用户侧开关。

参考由管理员经既有 Reference Profile 流程注册：

```bash
python -I -B -m encode_pipeline.cli.admin \
  --database-url sqlite:////path/platform.db \
  --reference-profile-config /path/reference-profiles.json \
  reference-profile register --safe-key <key> --display-name <name> \
  --organism <organism> --assembly <assembly> --config-key <key>
python -I -B -m encode_pipeline.cli.admin ... reference-profile verify <revision_id>
python -I -B -m encode_pipeline.cli.admin ... reference-profile enable <profile_id> --revision-id <revision_id>
```

私有配置 `reference-profiles.json`（0600）按
`helixweave-reference-profiles-v1` 编写，`hitrac-preprocess` 键的 payload 为
`hitrac-reference-profile-v1`，含参考 binding 文件路径与其 SHA256。参考
binding 固定六件 `.bt2` 小索引、FASTA 摘要与 contig 集合；仅 `.bt2l` 的索引
不被当前原入口接受。执行前、科学完成后与发布前都会重核这些字节。

worker 进程沿用既有入口（`python -m encode_pipeline.workers.cli` 或
`encode-worker`），启动时先执行迁移清单准入并打开既有 SQLite；API 与 worker
共用 `ENCODE_PIPELINE_DATABASE_URL`、`ENCODE_PIPELINE_WORKSPACE_ROOT`、
`ENCODE_PIPELINE_REDIS_URL`、`ENCODE_PIPELINE_QUEUE_NAME` 与
`ENCODE_PIPELINE_REFERENCE_PROFILE_CONFIG`。`ENCODE_PIPELINE_JOB_TIMEOUT_SECONDS`
是科学执行期限；RQ 外层兜底为其加 300 秒启动宽限与 30 秒清理宽限
（`workers/rq_queue.py` 的 `rq_job_timeout_seconds`），两者独立。
ProcessRunner 的 allowlist 只接纳经核验的 Hi-TrAC 运行时解释器（以及既有
ENCODE/Bulk 准入项），adapter 返回的任意 executable 不被信任。

## 3. 取消与超时链的固定依赖

取消确认依赖固定 RQ 2.10 停止消费点，升级 RQ 必须复核：

- `rq/command.py` 的 stop 命令先写停止标记再 kill horse；
- `rq/worker/worker_classes.py` 的 monitor 在 wait4 返回后消费该标记并调用
  停止回调；
- `rq/worker/base.py` 的原失败处理在同一 monitor 调用内存在第二消费点。

`workers/timeouts.py` 的 `DurableWorker` 在上述整个 monitor 上下文内要求每次
非空标记消费都匹配 job、horse PID/starttime 与完整嵌套清理证明；wait4 保持
锁外。未知、失根或失败一律拒绝确认，不承诺自动找回全部失根后代。平台取消
顺序为：API 先持久化 `cancellation_requested`（SQLite），再经 RQ pub/sub 发
stop；worker 收割完成后 `handle_execution_stopped` 才提交
`cancellation_acknowledged` 与终态 CANCELLED。HTTP 202、RQ STOPPED 或
kill_horse 返回单独都不构成取消完成。

超时走既有契约：ProcessRunner 到期杀死该 attempt 的进程树（含新建 session
的科学孙进程），run 终态 FAILED、reason_code `PROCESS_RUNNER_TIMEOUT`，不发布
任何结果；RQ 外层超时是兜底，不首先触发。

## 4. 结果、输出与通知契约

- 公开产物仅为每样本原 `*_all.bedpe.gz`、`*_unique.bedpe.gz` 和一份原
  `tracPre_summary.txt`；summary 为 15 个指标加样本索引共 16 个物理 TSV 列，
  索引列是内部 token（`s000001`…），平台仅在 metadata/QC 坐标映射回原用户
  样本名，不改原文件字节。BAM、裁剪 FASTQ、请求与诊断私有保留，不自动删除，
  不提供下载；原始 FASTQ 始终只读。
- 政策 A：任一样本 all/noBg 为空或需 cis 分母的 QC 集合无 cis，整批拒绝成功
  发布；不透传上游补 1 的计数，不做部分发布或 NA 映射。
- 发布是原子 bundle：artifact、QC、两者 generation、attempt、事件与
  append-only publication 同事务提交；回滚不留部分可见新结果。科学
  SUCCEEDED 与结果发布失败是不同状态；成功邮件只在本尝试精确完整 bundle 后
  发送（非 opt-in adapter 的既有通知政策不变）。
- gzip 字节不做跨执行一致性要求（时间戳字段）；科学比较使用契约规定的规范化
  PET 多重集合与指标。
- 通知为终端邮件：`HELIXWEAVE_TERMINAL_EMAIL_*` 与 `HELIXWEAVE_SMTP_*`；
  `local_plaintext` 模式仅接受 loopback 主机。成功邮件带当前 QC 摘要
  （最多渲染 12 条，其余聚合为计数行）。

## 5. 升级路径

- 运行时/工具：构建新 prefix 后生成新 runtime binding，更新两个
  `HELIXWEAVE_HITRAC_RUNTIME_*` 坐标并重启 API/worker。所有未完成任务在新
  准入下重核；漂移即拒绝执行（fail-closed），不产生混合版本结果。
  `tools.lock.json` 变化属于受控实现字节变化，必须经正式生成器同步执行身份。
- 参考：注册新 revision 并 enable；既有 run 绑定其快照时的 revision，
  不换绑、不重写历史结果。
- 平台实现：执行闭包内文件变化按 AGENTS.md 走正式 manifest/qualification
  生成；身份不一致时准入与发布前核验拒绝，run 失败而不是静默沿用旧字节。
- RQ/Redis 大版本升级前必须复核第 3 节的三个固定消费点坐标与停止证明边界，
  并重跑 worker 停止回归与 H5 取消链。

## 6. 上游耦合点（实测记录）

- `tracPre2.py` 的子命令退出码被 `cLoops2/utils.py` 忽略、缺工具/输入分支可
  退出 0、遇旧产物会跳过、summary 写两次：平台以全新 attempt、逐调用开始/结束
  记录与退出码、产物后验和 BAM 配对来源核验兜底；任何缺口整批拒绝。
- mapping ratio 依赖 Bowtie2 日志形状（固定脚本第 245 行读取、倒数第二行
  百分数）；换 Bowtie2 构建必须先核日志格式。
- 背景/Cis 判定与 QC 距离分组使用不同边界（浮点中点 `<1000` 与整数中点
  `<=1000/10000`）；PET 去重键为原 BEDPE 前六列字符串，QC 键另含 strand 与
  规范化端点；两者都不等同 Picard/samtools 去重。
- 无 cis/无 unique 时原 QC 将分母置 1：这些值不作可信科学数值发布（政策 A）。
- 科学并行：`threads`（2～8，默认 2）映射原脚本 `p`，samtools 另固定 `-@2`
  辅助线程；参数不是 CPU/RSS 硬上限。
- Redis/RQ：作业身份为 `run-execution-<sha256>`，停止经 pub/sub；通知 SMTP
  仅 loopback 明文或显式 TLS 配置。

## 7. H5 实测边界

本轮真实链（认证 HTTP API → validated snapshot → RQ/worker → 固定 tracPre2
与真实工具 → SQLite 原子发布 → 原 API/桌面/移动页面 → 实际下载）及取消、
超时、政策 A/原脚本失败拒绝均已通过；证据见本轮交付目录。以下不是本轮结论：
完整 Protected Bulk Gate（按用户决定留 H5 后统一执行）、真实 SMTP 投递
（仅 loopback 捕获）、远端 CI、真实实验数据的生物学有效性、科学 tiny 之外的
输入规模。F2 发布失败的真实 worker 内故障注入无安全确定性接缝：该精确顺序由
原服务/仓库加真实 SMTP 客户端到 loopback 捕获的定向集成覆盖，并在本轮报告
明确标注层级。

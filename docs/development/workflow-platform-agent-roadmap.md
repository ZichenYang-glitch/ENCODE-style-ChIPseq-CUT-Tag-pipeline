# HelixWeave Product Roadmap

This roadmap records the maintained product boundary, delivered baseline, and
explicit follow-ons. It does not promise dates. Detailed implementation history
belongs in Git history, the changelog, and release evidence.

## Product boundary

HelixWeave is a workflow-neutral omics platform for one laboratory operating
one Linux x86_64 or systemd-enabled WSL2 host on a trusted LAN. SQLite and the
local filesystem are canonical stores; Redis/RQ is an execution handoff, not a
second lifecycle authority.

The bundled registry contains:

- ENCODE-style ChIP-seq, CUT&Tag, ATAC-seq, and MNase-seq through Snakemake;
  and
- `bulk-rnaseq` through the pinned nf-core/rnaseq 3.26.0 Nextflow runtime.

Workflow schemas, commands, scientific behavior, output paths, and artifact/QC
extraction remain adapter-owned. Authoring may be available without runtime
assets; execution remains fail-closed until the complete operator binding is
admitted.

HelixWeave is not a hosted multi-tenant or high-availability service. Trusted-
LAN authentication is an administrator/member boundary, not general RBAC or
multi-tenant isolation.

## Delivered baseline through PR #179

The current baseline includes:

- workflow-neutral schema authoring, structured validation, immutable
  snapshots, durable lifecycle/events/logs, cancellation, restart recovery,
  artifact/QC extraction, safe downloads, and responsive browser journeys;
- Project, Sample, SampleRevision, StoragePool, InputFile, InputFileRevision,
  and immutable input-use provenance;
- administrator-managed revisioned Reference Profiles bound to snapshots and
  runs;
- cross-run artifact publication and provenance;
- administrator diagnosis, exact-identity fail, and one-use requeue recovery;
- an offline Linux/systemd deployment operator with independently managed
  platform and scientific runtime slots;
- trusted-LAN accounts, sessions, CSRF protection, administrator/member roles,
  audited privileged actions, and terminal email preferences;
- staged FastQC/MultiQC evidence for the ENCODE-style workflow and dynamic
  current-generation QC notification summaries;
- the compact Rosemary desktop/mobile browser experience;
- pinned Node 22.23.1/npm 10.9.8 authority and independent CI, Lint, Lock
  Check, coverage, frontend, browser, platform-real, scientific-real,
  container, and Protected Bulk evidence tiers; and
- a read-only Omics Intake Bundle 0.2 inspection boundary and a read-only Agent
  explanation boundary.

Controlled tiny and synthetic evidence demonstrates integration, execution,
lifecycle, offline reuse, and result contracts. It is not biological
validation, public-data reproducibility, or production-scale performance.

## PR #180 release closure

PR #180 closes the v0.4.0 release candidate without adding scientific or UI
features. It synchronizes the release identity, canonical generated assets and
contracts, released-schema migration evidence, public release documentation,
and the three-asset publication boundary. It also provides the narrow public
offline producer that composes the existing platform, ENCODE runtime, and Bulk
RNA-seq deployment bundle formats from release artifacts and explicit
operator-owned caches.

The candidate is complete only after ordinary exact-HEAD checks, wheel/sdist
and two no-checkout installs, release-shaped bundle production, the released
schema-08 to schema-15 migration proof, one separately authorized supported-
host deployment window, and—if the final implementation manifest requires
it—one exact-HEAD combined Protected dispatch. Date freeze, merge, tag, and
publication remain separate human gates.

## Explicit follow-ons and non-goals

The following do not block v0.4.0 and require separately scoped decisions:

- further Agent capabilities or any Agent write action;
- PostgreSQL, object storage, Kubernetes, multi-host or high-availability
  operation, Slurm/HPC/cloud executors, or remote workspace semantics;
- Singularity/Apptainer support or published Docker/OCI images;
- IGV.js or other new scientific/browser capabilities;
- dependency-major maintenance, CSP/proxy hardening beyond the trusted-LAN
  contract, and additional UI polish; and
- redesigning Bulk qualification, CI, runner, release automation, or the
  existing deployment bundle/state contracts.

Do not change the `encode_pipeline` import namespace, compatibility CLI names,
repository slug, workflow identities, or artifact URI scheme without an
explicit compatibility decision.

## Post-release maintenance priorities

Sequenced follow-ups identified during the 2026-09 local bring-up and the
accompanying full-stack code review. Defect evidence and root-cause detail
live in `docs/development/bug-log-2026-09-local-bring-up.md`; this section
holds only outcome statements and their exit evidence.

Tier 1 — close the real-path gaps first:

- Completed implementation: real-submission smoke coverage in `4c64f0b`.
  The browser CI tier materializes the served authoring schema through the
  workbench's rjsf logic, validates bulk defaults with the real adapter, and
  executes the controlled ENCODE fixture with `cores: 2`, asserting both argv
  and Snakemake's allocated threads. The tests guard bug #1 and bug #8; they
  do not replace a real bulk scientific run or the Protected Bulk Gate.
- Local deployment documentation now records the rootless/containerd bring-up,
  the two-daemon workaround, explicit container identity, Reference Profile
  registration and local verification toolchain in
  `docs/development/local-platform-runtime.md`. Remaining exit evidence: a
  fresh host reaches `doctor` green from the document alone, with bulk
  availability checked separately; shim retirement and a real rootless run
  remain operator verification tasks. The unresolved staging/admission storage
  mismatch is tracked in bug #5.
- Dual run-repository conformance: one parametrized behavior suite executed
  against both the in-memory and SQLAlchemy run repositories. Exit evidence:
  shared suite runs in CI; drift between the two implementations fails the
  build.
- Observability for deliberate silence: structured, payload-free failure
  breadcrumbs in the `_safely` notification and doctor paths. Exit evidence:
  a killed notification channel is diagnosable from logs without exposing
  private payloads.

Tier 2 — batched contract work, one Protected gate:

- Completed implementation: stage-naming retirement per
  `docs/architecture/stage-naming-retirement-plan.md`, container uid/gid as
  explicit deployment coordinates (bug #7), ENCODE command ownership moved
  back into the adapter, and the ENCODE `--cores` fix (bug #8), delivered in
  `2c26a9e`, `55ef928`, `02303b2`, and `0f01241` respectively. A combined
  Protected Bulk Gate for the batch remains pending; completion here records
  the implementation and targeted checks, not an unrecorded gate result.
- Completed implementation: frontend schema-driven master-switch cascade
  (bug #4), including QC sub-flags and the declared neutral states for UMI/rRNA
  removal. Form and desktop/mobile browser tests cover the cleared request data;
  backend fail-closed validation remains unchanged.
- Upstream coupling ledger: documentation completed, pending acceptance;
  [13 entries](upstream-coupling-ledger.md) record MultiQC sample-name handling,
  the nf-core parameter allowlist, and consumed QC formats, with pinned sources,
  local consumers, test boundaries, and concrete upgrade checks.

Tier 3 — investigations:

- Unify Docker storage semantics between staging and runtime admission so a
  single daemon can serve both (bug #5).
- Gradual slimming of the known hotspots (`persistence/repositories.py`,
  `services/run_repositories.py`, `adapters/bulk_rnaseq/runtime_assets.py`,
  `adapters/encode/manifest/make.py`) by extracting the touched cluster
  whenever work lands there — no standalone rewrites.

## 已确认的维护顺序与审核边界（2026-09-22）

本节记录用户接受第二轮独立复核后的计划，作为当前交付顺序；上方
Tier 分组保留为维护背景，不再作为重排依据。保持 PR-2～PR-17 原顺序，
不新增 PR。PR-2 的调查交付已获用户接受；PR-3 的实现与定向测试已获用户
接受，新身份的 Protected Bulk Gate 仍待补验；PR-4 实现和定向验证已通过
独立审核，用户允许继续，新身份 Gate 仍待补验；PR-5 文档返修已通过独立复验，用户允许继续；PR-6
返修的实现与定向验证已通过独立复验，用户允许继续 PR-7；完整在线 staging、
部署 admission、科学小样本及各身份 Gate 尚未通过。PR-7 已按用户批准的固定 hint 方案完成实现与定向验证；快照失效后成功提示的限定小返修已通过独立复验，用户已接受；
PR-7 新身份 Gate 仍待补验。PR-8 的规则依赖修复与定向验证已通过独立复验，用户已接受。
PR-9 的实现与定向验证已通过独立复验，用户已接受；Usage 文案收尾完成。
PR-10、PR-11、PR-12、PR-13、PR-14、PR-15 的实现与定向验证已通过独立复验并获用户接受；PR-15 新身份 Gate 仍待补验。PR-16 的实现与定向独立复验已接受，新身份 Gate 欠验单列。PR-17 主体及三处文案返修均已通过独立复验；fastp 字段修复已接受。成功流式下载超时调查已通过独立复验：当前受限执行环境阻断 asyncio self-pipe 写入；允许本机 socket 后原测试及真实 ASGI/loopback 下载通过，无需生产或测试代码修改。最终集成与候选准备状态见后文第四收尾单元，不能据此关闭历史 Gate 欠验。
编号为本次维护队列编号，不是远端仓库的 pull request 编号。
2026-09-23 用户确认先完成剩余事项再收尾；后续 Hi-TrAC adapter 单独推进，
具体安排见后文“当前批次收尾与 Hi-TrAC adapter 接入”。

实施方负责调查、修复及提交验收证据；本审阅会话只负责独立审核，
不承担 PR 实施。一次只推进一个 PR，完成后停下等待用户验收，
不自动提交或推送。PR-2、PR-6 先调查；未确认需要代码修复时，
以调查结论、原始证据和未验证边界收尾，不为产生 diff 而改代码。

**PR-2 状态（2026-09-22）：调查完成，用户已接受，当前未复现。** 原 MACS3
规则的小输入执行、两条命令路径的传递及托管激活器组件验证未确认生产缺陷，
因此没有实现修复。历史失败 run/traceback 和已确认部署的 ENCODE 托管运行时
均缺失；未验证完整托管 Snakemake/RQ 执行，不据组件通过宣布历史问题不存在。
历史 bug #3 保持 **OPEN (unconfirmed)**，没有改为已修复。
证据及边界见 [bug log #3](bug-log-2026-09-local-bring-up.md#3-open-unconfirmed--encode-run-macs3-conda-environment-activation)。

**PR-3 状态（2026-09-22）：实现与定向测试已接受，新身份 Gate 待补验。** 用户已授权
有限范围的异常契约对齐：SQL `_insert_event` 在转换 context 前检查 Mapping，
非 Mapping 抛 `ValueError("context must be a mapping")`，合法非 dict Mapping
仍可写入。契约限于其他调用前提满足时，本次覆盖的创建、更新、preflight
事件写入拒绝非法 context 且状态不变；不前移到 Draft 构造，不统一所有内部
异常，不新增 JSON 限制。没有证据证明此前合法业务请求失败或事务损坏。

原有断言保留，补齐两后端相同的创建/读取、重复标识、状态 CAS、缺失对象、
事件/日志顺序与分页，以及经 RunService 校验的生命周期和终态重试覆盖；
复用既有 history/recovery 与数据库专项并发、重开测试。共享测试 62 passed，
相关持久化及执行身份测试 98 passed；现有 CI shard 2 收集共享测试，无需
修改 CI YAML，未宣称远端 CI 已通过。

`repositories.py` 属于执行闭包，已重建 111 文件 manifest 及配套 qualification；
仅该受控源文件条目变化，持久化契约/schema projection 不变。新 aggregate 为
`4baa74b571390d7f5faf0803bbc768fe28c21489ed2c9c665f7ec46eee943f66`。
Protected Bulk Gate 原平台验收入口已尝试，缺少显式 runtime root 而在前置
检查失败；fixture manifest、测试 Redis、Docker executable/socket 坐标同样
未配置，且当前未提交工作区不满足 protected CI 的精确已审阅 clean SHA 要求。
未进入真实执行，不能视为 Gate 通过或所有验证完成；用户接受实现和定向
测试并允许带此欠验继续 PR-4。后续须绑定 PR-3 精确身份补验，不能用
PR-4 或后续身份的通过结果追认。
证据：`/tmp/helix-pr3-align-ohpg1w_q/`；上一轮失败证据保留在
`/tmp/helix-pr3-repository-hi7wfzmq/`。G0 的旧身份 Gate 欠账仍独立跟踪。

**PR-4 状态（2026-09-22）：实现和定向验证已接受，新身份 Gate 待补验。**
修改前已将 PR-3 的 111 个受控文件及原 manifest/qualification 逐字节保存到
`/tmp/helix-pr4-logs-ply3x7gf/pr3-recovery/closure/`，并独立核对大小、SHA-256
和 aggregate；`identity.json` 记录 Git 基线及每个文件的哈希和权限，
`README.md` 说明恢复方法。该副本恢复的是受控源码闭包与资格记录，不包含
私有部署配置或外部运行时，不是 Gate 通过证据。

本 PR 使用标准库 logging 的固定 message 与 `component`、`phase`、
`reason_code` 字段，覆盖通知 `_record_outcome` 内部事件记录失败、
RunService 外层 notifier 异常及 worker 成功通知/超时保护分支。
普通邮件发送失败仍仅产生已有持久事件，不重复补写、不重试、不改终态。
当前 doctor 没有名为 `_safely` 的函数；经用户确认，覆盖 deployment
doctor 的 `_run_probe` 及 local-platform JSON doctor 的 environment、
workflows、recovery 三个异常降级点，保留原 ProbeResult、JSON 和退出码。
日志不传入异常对象、异常文本、traceback、载荷、配置、环境值或私有路径；
阶段只来自固定调用点或固定 doctor 检查清单，日志写入异常仍为 best effort。

故障注入与正常对照 24 passed；相关通知、doctor、RunService、恢复、worker
通知顺序及执行身份回归合计 248 passed，0 failed/0 skipped。修复前的同组
日志断言有 14 failed/8 passed（后补两项已知 DeploymentError fallback 对照）；
失败原因是缺少诊断日志。实际格式化日志、结构化字段、异常信息缺失和正常
路径无失败日志均已验证。两个使用 PYTHONPATH 的旧导入探针未运行；doctor
打包/stage 测试、真实 SMTP、Redis、生产数据库及大型科学流程未运行。
现有 CI 分片收集这些测试，无需改 CI YAML；未宣称远端 CI 已通过。

实际闭包中变动的是 `services/runs.py`、`services/terminal_notifications.py`、
`workers/jobs.py`、`workers/terminal_notifications.py`，已按规定入口更新
111 文件 manifest 和 qualification；两个 doctor 模块不在该受控清单中。
PR-4 aggregate 为
`4d23e74017f46339d5fea7f3b491641c7cfe279867776be1edb2722348119962`，
持久化契约及 schema projection 未变。由于执行闭包改变，需要为本身份
保留 Protected Bulk Gate 验证项；runtime root、fixture manifest、测试 Redis、
Docker executable/socket 坐标仍缺失，本轮未重复运行已知会在前置检查失败
的 Gate 入口。**PR-4 新身份欠验、PR-3 的 4baa74…943f66 欠验与 G0 是三项
独立记录**，后续结果不可互相追认。ruff 两项及限定文件 diff 检查通过。
本轮证据、独立 diff 和原始命令输出：`/tmp/helix-pr4-logs-ply3x7gf/`。
独立审核证据：`/tmp/helix-pr4-review-ncvkzk7m/logs/diagnostics-command.json`、
`identity-and-recovery-check.json`。用户允许推进 PR-5，不代表上述 Gate 已通过；
已有 PR-3 恢复材料与 PR-4 证据保持原样。

**PR-5 状态（2026-09-22）：上游耦合台账文档返修已通过独立复验，用户已接受。**
新增 [上游耦合台账](upstream-coupling-ledger.md)，沿用既有 adapter 决策与
固定契约文件，记录 UC-01～UC-13：bulk MultiQC 的清洗/替换与公开机器表、
nf-core 参数分类及转换、FastQC/Trim Galore/fastp/STAR/Salmon/featureCounts/
Picard/RSeQC 格式、阈值表接线，以及独立的 ENCODE MultiQC 1.35 配置。
每项列出版本与官方来源、本地函数/规则行号、已有测试的实际覆盖、升级命令
和证据边界；本地有意限制或与上游的差异不自动登记为缺陷。

初次交付核对固定 NF commit 的下载源码与 source manifest，以及 vendored schemas、
MultiQC 1.33 清洗/分组源码和 RSeQC 5.0.4 源码包哈希；只下载少量官方源码，
未安装或执行上游工具。经 fixture 检查及受保护路径/数据库/网络访问审计，
canonical bootstrap 执行 23 个已有合成契约测试函数，展开为 **51 passed、
0 failed、0 skipped**；不是科学流程或真实 MultiQC 执行证据。
仍未核实部署 OCI digest、worker/本机工具版本、Subread 2.0.6 底层源码
（本轮官方 URL 返回 404）、NF CSI utility 接线及完整报告/多 library 实测。

本轮仅修改台账与 roadmap，不改执行行为、依赖、manifest 或 qualification，
不重建身份、不打包、不运行 Bulk Gate；PR-3、PR-4 与 G0 欠验继续分别保留。
本轮原始来源、命令输出、引用核对及独立 diff：`/tmp/helix-pr5-ledger-si4g0_12/`。

独立审核要求的两项文档返修已记录：结果提取属于 `BulkRnaSeqResultsWorkflowAdapter`
子类；UC-06 明确 fastp 1.0.1 将版本写在 `summary.fastp_version`，本地解析器
却要求顶层字段，既有合成夹具存在同样盲区。审核的原函数对照只移动版本键，
合法 fastp 配置校验均通过；上游布局被解析器和公开 QC 提取拒绝，探针结果
**1 passed / 1 failed，退出码 1**。既有 51 项测试通过不能否定该兼容差异。
artifact 消费链仅源码核对；默认 Trim Galore 路径不因该差异触发。
本次返修未复跑行为测试、未运行真实 fastp/容器/完整流程，部署版本仍未核实。
撤回初次交付报告中没有确认生产问题的笼统结论；用户决定延期修复 fastp，修复归属尚未决定，
本 PR 只改文档，不新增 PR、不调整顺序。审核证据：`/tmp/helix-pr5-review-6b0n09oq/`；
返修报告、检查与相对本轮基线的独立 diff：`/tmp/helix-pr5-doc-revision-3e8r8g68/`。
用户已允许继续 PR-6；fastp 不并入本轮，也不新增 PR 或调整队列。

**PR-6 状态（2026-09-22）：返修实现与定向验证已通过独立复验，用户已接受；完整验证仍待补验。**
真实 Docker 29.1.3（Engine `fbf3ed2`）对照确认：classic 存储的 inspect ID 为
config digest，containerd 为目标 manifest digest；保存同一合成镜像时，classic
转为未压缩 OCI 归档，containerd 保留 gzip 层及原 manifest。旧 staging 的
config-ID 断言和仅接受未压缩 OCI 的归档读取器，不能配合现有 containerd
admission。这里的归档身份是 `index.json.manifests[0].digest` 指向的 manifest，
不是 index 文件自身摘要；archive 文件 SHA 与 config SHA 又分别独立。

单 daemon 方案沿用操作者显式选择的 containerd endpoint，不增加环境 fallback：
staging 严格核对选定的 linux/amd64 manifest ID、按裸 digest 导出；归档消费增加
Docker v2/OCI gzip 表示支持，保留 config 字节、压缩层 SHA、解压 diff ID、大小上限、
禁止 RepoTags 等校验。旧未压缩归档仍受支持；admission 的精确 image ID/rootfs
检查、worker endpoint 匹配及生命周期不变。只改两处生产文件及对应测试，
没有修改 daemon 配置、迁移存量数据或调整公开 schema。

首轮交付回归先有 **3 failed**，修复后转绿。当轮 canonical bootstrap 定向集合
**341 passed / 0 failed / 0 skipped**；真实归档校验、原 Docker admission 组件及
ProcessRunner 的隔离小镜像执行 **1 passed**，输出固定标记，错误 rootfs/daemon
仍拒绝。该实验不是完整 Nextflow/RQ 或 nf-core 小样本验收。在线 pull 在 daemon
现有代理连接阶段失败，未完成真实在线 staging → 完整 runtime admission；
未改代理或 daemon 配置。两个自建测试镜像已移除，既有镜像 ID 集合与开工一致。

`runtime_assets.py` 改变受控闭包，已按规定入口重建 111 文件 manifest/qualification，
首轮 aggregate 为 `d6b44c855d99d5bf4016c8cc63811abb9a059f737a61cc58a596a27f5a6bb6ee`。
受控条目仅此文件变化，32 项身份测试在上述集合中通过；生成资格记录不等于 Gate 通过。
五项 Gate 运行坐标仍未配置，也没有本轮对应的已审核 clean commit；未进入完整
runtime admission、Redis/RQ、科学/产品 Gate。**PR-6 本身份欠验独立于 PR-3、
PR-4 与 G0**，不能互相追认。旧 111 文件与原 manifest/qualification 的逐字节副本
已保存在本轮 `baseline/`，`closure-recovery.json` 记录哈希；旧证据未覆盖。

交付、原始命令/退出码、探针纠错、身份对应表、限定文件检查及本轮独立 diff：
`/tmp/helix-pr6-docker-w93xj3bx/`。该首轮交付当时停下等待审核；后续返修验收状态
如下，不声称历史完整部署已复验，不自动暂存、提交或推送。

**本次返修：** 独立审核确认新增的 archive descriptor/layer `mediaType` 集合判断
对 JSON `[]`/`{}` 抛出未捕获 `TypeError`，使非法归档的结构化拒绝退化；这不是
正常 Docker 产物必然失败，也不是镜像完整性绕过。仅在
`runtime_assets.py:1649,1715` 增加字符串类型守卫，继续返回
`BULK_RNASEQ_RUNTIME_ASSET_CONTRACT`；不改既有 distribution 边界、daemon
传递、精确身份、gzip SHA/size、解压 diff ID/大小限制或其他拒绝规则。

正式测试 `test/adapters/test_bulk_rnaseq_runtime_assets.py:2456` 通过公开
`verify_runtime_asset_closure` 验证四种组合，保持归档登记摘要及目标 blob 身份一致。
修复前 **4 failed**（均为 TypeError），修复后 **4 passed**；最终三个指定文件
合计 **163 passed / 0 failed / 0 skipped**（staging 20、runtime assets 111、
execution identity 32），保留 OCI/Docker v2 gzip 正向、损坏及大小限制覆盖。
Ruff 和限定六文件 diff 检查通过。canary 的显式 `/tmp` 临时目录仅由验证 harness
按固定前缀重定向到本轮目录，未修改生产 canary 或测试断言。本轮未运行真实 Docker、
在线 staging、部署 admission 或科学小样本；首轮记录的完整验证边界继续保留。

返修后按原入口生成 111 文件 manifest/qualification，aggregate 为
`15646a6a76bad8d8ad7df6855d2abe0b179fcd186c233000311633ab8892e27c`；受控源文件
仅 `runtime_assets.py` 变化，持久化契约字节不变。五项本地 Gate 坐标仍缺失，
未重复启动已知缺前提的 Gate；新身份欠验与 PR-6 首轮身份、PR-3、PR-4、G0
分别保留，不跨身份追认。返修前精确闭包、manifest/qualification 及恢复哈希位于
`/tmp/helix-pr6-media-fix-s63ji4z_/baseline/`（`recovery.json`），未覆盖旧证据。
本次报告、红绿日志、身份检查和相对本轮开工基线的独立 diff 位于
`/tmp/helix-pr6-media-fix-s63ji4z_/`。返修现已通过独立复验，用户接受实现与定向验证，
允许继续 PR-7；完整在线 staging、部署 admission、科学小样本及 Gate 的欠验仍保留。

**PR-7 状态（2026-09-23）：实现与定向验证（含限定小返修）已通过独立复验并获用户接受；新身份 Gate 待补验。**
用户已批准 `/tmp/helix-pr7-conflicts-4u30tn6r/proposal.md` 的限定方案。
仅在 adapter 既有六处拒绝分支填入固定 `Issue.hint`：UMI/trimmed-read outputs、
`min_trimmed_reads`、`rseqc_modules`、`deseq2_vst`、featureCounts biotype 参数。
错误码、path、拒绝条件、用户配置和公共字段均保留；schema、OpenAPI、生成客户端未改。
QC master 重新开启后，子开关仍可能关闭；`deseq2_vst=false` 仍是存在的参数。
提示准确说明这些前提，不承诺仅打开 master 即可消除全部冲突。

通用前端按已服务 schema/path 定位字段，无法直接聚焦时回退到已知 section 或
YAML 编辑器；提示在页签外可发现，并显示在对应表单位置。hint 仅作文字展示，
不解析其内容，不硬编码 bulk 依赖。编辑使旧校验诊断失效，延迟响应不覆盖新配置；
创建结果提示独立保留，包括 `RUN_CREATE_INPUTS_CHANGED` 与 `RUN_CREATE_UNCONFIRMED`，
既有快照失效、创建请求和不确定结果的重试语义未改。无自动跨 section 改值。

验证证据位于 `/tmp/helix-pr7-implementation-b68brofy/`：

- 新增 adapter 定向回归 **7 failed → 7 passed**；相关全文件 **329 passed**。
- 原 HTTP 路由、真实 adapter/reference resolver、隔离 SQLite：六种冲突与显式纠正
  **6 passed**；临时应用没有已准入科学运行时，合法输入成功但不签发可执行快照。
- 通用 UI、schema/draft 与既有 workbench 路由 **120 passed**；含延迟成功/失败响应、
  创建结果保留、真实 RJSF 禁用控件与 YAML 回退。typecheck 通过。
- 从 `frontend/` 工作目录生成的完整样式构建，桌面 1440×900、移动 390×844
  **各 3 passed**：真实 HTTP 错误、QC/UMI 关闭、键盘定位/纠正、配置保留和延迟响应。
  无横向溢出；已复核移动端错误及 outputs 截图。专用 `@validation` 场景使用
  `test/browser/workbench_validation_runtime.py`，本轮通过证据不等于远端 CI 通过。
- 依照现有 CI 资产闭包要求同步预编译资产；先生成到本轮 `/tmp`，确认原资产与
  开工基线逐字节一致后更新。最终资产身份 `sha256-9721877881da9eb52f08490b4ef16b08160e8e0fd8aca15d774aba08e5755bef`，
  前端部署资产测试 **21 passed**。依赖和 API 契约摘要不变。
- 早期探针错误（options/tool 名称、schema/DOM 测试假设、只读 SQLite URI 守卫、
  回环服务网络隔离、登录限流、CodeMirror 替换、错误构建 cwd 漏载 Tailwind）均保留
  原失败记录并纠正。无样式构建的浏览器通过不作为最终布局证据。

`validation.py` 是本轮唯一变化的既有 111 项执行闭包源码；已按规定入口更新
manifest/qualification，**32 项身份测试通过**。新 aggregate 为
`b6751715f1380ad5ae7ae993cf4a2c417602b24de81e37fdf6fed982492b3cf8`，
持久化契约字节不变。旧 `15646a…e27c` 的 111 项源码及 manifest/qualification
精确副本在本轮 `baseline/`，恢复哈希已核对。
五项 Gate 坐标（runtime root、fixture manifest、测试 Redis、Docker CLI/socket）
仍缺失，没有本轮已审核 clean commit；未启动已知缺前提的 Gate，未进入真实科学执行。
**PR-7 新身份欠验独立于 PR-3、PR-4、PR-6 各身份及 G0，不跨身份追认。**
fastp 继续延期，修复归属待用户裁决；PR-7 交付时未启动 PR-8，维护顺序不变。
完整命令、退出码、原始输出、构建探针边界及本轮独立 diff 见本轮 `report.md`、
`logs/`、`this-round.diff`；此前调查证据保留在 `/tmp/helix-pr7-conflicts-4u30tn6r/`。

**PR-7 返修（2026-09-22）：** 独立复验确认创建请求收到
`VALIDATED_SNAPSHOT_EXPIRED`/`VALIDATED_SNAPSHOT_STALE` 后，父组件成功校验提示未撤销，
与创建错误并存（提示回归，非创建绕过或重复创建）。返修在 `ValidatedSubmission.tsx`
失效分支按 snapshot 归属替换对应成功 notice，明确要求重新校验；创建错误诊断、
`RUN_CREATE_UNCONFIRMED`/`RUN_CREATE_INPUTS_CHANGED` 保留语义、用户草稿与参考选择均不变；
未改后端拒绝条件、schema/OpenAPI/生成客户端与执行闭包，未重建 manifest/qualification，
不新增 Gate 要求。正式回归：`input-workbench-route.test.tsx` 扩展 EXPIRED 用例并新增
STALE 用例（创建前成功提示出现；失效后错误诊断保留、旧成功提示消失、按钮禁用、
草稿与参考保留、重验后恢复有效状态）；前端单测 **392 passed**、typecheck 通过；
预编译资产按既有流程从 `frontend/` 重建并同步，新资产身份
`sha256-ec23dc88c3e60a081c6148db159973bc99049c5ed853036cb3918bd1d51f60bf`，
资产测试 **21 passed**、执行身份测试 **32 passed**（aggregate 不变）；`@validation`
桌面/移动浏览器 **6 passed**（真实后端加载本轮资产；该栈无 admitted runtime，
EXPIRED/STALE 创建路径由 jsdom 模拟响应用例覆盖，不称为真实后端过期验证）。
证据与本轮独立 diff 位于 `/tmp/helix-pr7-notice-fix-round1/`；复验证据保留在
`/tmp/helix-pr7-independent-review-5yscsjvg/`。PR-7 及各历史身份、G0 的既有 Gate
欠验原样保留。

**PR-7 最新验收（2026-09-23）：** 用户已接受限定小返修。
独立复验 `/tmp/helix-pr7-notice-retest-f65cu99r/report.md` 实跑相关前端
**121 passed**、原独立探针 **8 passed**、资产测试 **21 passed**，typecheck 通过。
该复验未重跑浏览器、科学流程或 Bulk Gate；PR-7 及各历史身份、G0 欠验保留。

**PR-8 状态（2026-09-23）：依赖修复与定向验证已通过独立复验，用户已接受。**
修复 `workflow/rules/qc.smk:927`：`project_qc_summary` 的输入从
`TREATMENT_SAMPLE_IDS` 改为 `PEAK_SAMPLE_IDS`，与单样本目标及产物目录的
`peak_centric` 契约一致。mixed treatment 输入仍合法，summary 默认值和 MNase
MACS3 guard 保留；未改汇总脚本、输出格式或 MNase 专用 QC/summary，未涉及 PR-14。

正式回归 `test/workflow/test_project_qc_summary_dag.py` 运行原 Snakefile：
修复前 **2 failed / 4 passed**（mixed 默认与显式 summary=true 均触发 MNase
MACS3 guard），修复后 **6 passed**。六组包括 mixed 默认、显式开启、关闭，
MNase-only、两个 ChIP treatment 的 peak-only、缺 FASTQ。原样重放修复前保存的输入：
前五组分别 **56 / 56 / 52 / 27 / 59 jobs**；缺 FASTQ 仍单独返回
`MissingInputException in rule trim_galore`。非 strict-inputs 配置与样本校验六组均通过。
测试断言 peak caller、peak counts、FRiP、通用汇总只覆盖 peak treatment，项目汇总
准确消费其单样本汇总；MNase 专用汇总、片段/信号目标、完成标记及结果清单依赖保留。

本轮使用非空合成 FASTQ、独立 experiment、临时 blacklist/chrom sizes，整个 qc 块
省略或仅设 summary。为隔离调度关闭 MultiQC、replicate analysis、control，非“全默认流程”；
Bowtie2 prefix 仅作为 DAG 参数，未构建索引、未运行比对或科学 shell。
相关既有 validator、MNase/pooled DAG、汇总 TSV、产物目录与构建身份测试合计
**123 passed / 0 failed / 0 skipped**（含新增六组）。首次扩展收集因既有
`scripts` namespace 导入顺序失败；按产物目录测试先行的顺序复跑同组测试通过，
未改断言或用 PYTHONPATH 绕过。Ruff 与限定文件 diff 检查通过；单文件 snakefmt
因本地缺 shfmt 返回 123，完整格式检查未通过，不将其记为成功。

`qc.smk` 属于 ENCODE 动态源码指纹范围；开工旧字节与新字节摘要已保存。
未调用会读取受保护 default profile 的完整 ENCODE 指纹入口，未宣称已获得部署构建摘要。
Bulk 清单 111 项均逐字节匹配，aggregate 仍为
`b6751715f1380ad5ae7ae993cf4a2c417602b24de81e37fdf6fed982492b3cf8`，
本轮三个改动路径均不在该闭包；未改 Bulk 配置、引用或 artifact/QC 语义，
依据 AGENTS.md 的影响判定不新增 Bulk Gate，也不重建 manifest/qualification。
PR-3、PR-4、PR-6、PR-7 各身份及 G0 的 Gate 欠验不变，不跨身份追认。

证据、命令与原始输出、限定开工基线及独立 diff 位于
`/tmp/helix-pr8-mixed-summary-4s5paalk/`。所有执行为定向测试或 DAG 证据，
不声称完整科学流程或远端 CI 已通过；新测试由既有确定性分片收集，无 full_main 排除。
独立复验 `/tmp/helix-pr8-independent-review-2c7fol1o/report.md` 实跑
**123 passed**，六组 DAG 回归通过；未执行真实测序流程。此次接受不代表历史 Gate
已补验。fastp 继续延期，PR-2～PR-17 顺序不变。

**PR-9 状态（2026-09-23）：S6/S7 实现与定向验证已通过独立复验，用户已接受。**
`nrf_pbc` 与 `preseq_complexity` 均改用已有的坐标排序、MAPQ/flag 过滤后且
去重前 `{sample}.{MAPQ_TAG}.bam`；preseq 按已验证 layout 仅为 PE 添加 `-P`。
NRF/PBC 仍默认开启，preseq 仍默认关闭。未改变过滤或去重政策、duplicate flag
保留行为、片段键定义、零分母 `NA`、FRiP/峰/信号及 dup_metrics 派生指标的输入。
preseq 不增加 `-Q` 或外推参数，工具失败仍导致原规则失败，不生成替代曲线。

修复前共享 DAG 断言 **8 failed / 1 passed**；原规则实际执行的 SE/PE 用例均因
去重丢失频数而在 preseq 失败（同轮其余 **7 passed**）。修复后相关 validator、
DAG、汇总、产物和身份回归 **138 passed / 0 failed / 0 skipped**，包含 PR-8 六组。
真实工具层 **15 passed / 0 failed / 0 skipped**：samtools 1.23.1 的真实去重链、
SE/PE 手算数值、空/无重复/clipping 边界、同 PE BAM 加减 `-P` 对照及无法外推边界。
八组原规则执行覆盖 SE/PE × yes/no/auto narrow/auto broad，分别单独请求两项指标，
确认不需要 final BAM 或 duplicate_handling 创建目录。真实过滤保留 0x400，排除
低 MAPQ、secondary、supplementary 对照；DAG 同时核验非默认 MAPQ 17/42。
preseq 3.2.0 在 1070 fragments/658 keys 的夹具上用生产默认外推参数成功；
诊断对照只另加 `-v`，没有用旧探针的 `-Q/-e/-s` 冒充生产规则。

新 DAG 测试进入既有 fast 分片；真实工具测试进入既有 `real_execution` 层，
为该层接入已有 `preseq.lock` 和显式 PRESEQ 路径，未改锁文件、未本地安装依赖。
该 CI 层仍按既有 dispatch/schedule/release 触发；未声称远端 CI 或 conda 激活已验证。
Ruff 和限定文件 diff 检查通过；Snakemake lint 保留全部 65 条警告，仅同步三处行号，
无新增警告。单文件 snakefmt 因本地缺 shfmt 返回 123，完整格式检查未通过。

Bulk 清单 111 项匹配，aggregate 仍为
`b6751715f1380ad5ae7ae993cf4a2c417602b24de81e37fdf6fed982492b3cf8`；
本轮改动不在该闭包，未改 Bulk 配置、参考或 artifact/QC 契约，不新增 Bulk Gate，
不重建 manifest/qualification。ENCODE 的 `qc.smk` 动态指纹输入已改变，旧字节与
前后摘要已保存；完整入口会读取保护 profile，本轮未调用，未获得完整 ENCODE 构建摘要。
PR-3、PR-4、PR-6、PR-7 各身份及 G0 欠验分别保留，不跨身份追认。

证据、原始命令输出、夹具、限定开工基线与独立 diff：
`/tmp/helix-pr9-complexity-onf9kg2b/`。实测限于小型合成 BAM、原过滤/指标规则与
真实 samtools/preseq；未运行比对、Picard、完整测序流程或验证真实数据预测质量。
预置排序 BAM 不冒充上游比对证据；其他 QC 开关关闭的隔离配置不称为全默认流程。
独立复验 `/tmp/helix-pr9-independent-review-uhcqpzj_/report.md` 重跑定向测试
**138 passed**、真实工具 **15 passed**，另有未标记重复的 SE/PE 独立对照；
首轮精简 PATH 缺 conda 的探针错误修正后通过，不算产品缺陷。

**PR-9 Usage 收尾（2026-09-23）：** 按用户要求仅修改
`scripts/calc_nrf_pbc.py` 模块 docstring：示例改为 `sample.mapq30.bam`，注明
MAPQ 只是示例且由配置决定，输入应坐标排序、过滤后并保留重复。AST 去除模块
docstring 后与开工字节解析结果一致；Ruff 与限定 diff 检查通过，未重跑科学套件。
阶段基线、SHA、独立 diff 和检查位于
`/tmp/helix-pr9-usage-pr10-lbwpfx6w/pr9-usage/`。脚本文案仍改变 ENCODE 动态
身份输入字节；完整身份捕获会读取保护 profile，未执行。Bulk 受控闭包及契约未变，
不重建清单；各历史身份与 G0 欠验保留。此阶段完成后才建立 PR-10 基线并继续，
fastp 延期及既定顺序不变。

**PR-10 状态（2026-09-23）：S2 实现与定向验证已通过独立复验，用户已接受。**
`workflow/rules/consensus.smk:300` 仅将 narrow 调用改为
`--final-output={params.final_output:q}`。空值成为合法的 `--final-output=`，
非空值仍受 shell 引号保护；没有改变 helper、final_method、summary 字段或算法。
broad 和 SEACR 的合法路径返回非空，此轮未改其调用，也不将其认定为同一缺陷。

正式回归 `test/workflow/test_consensus_execution.py` 使用两个独立 treatment bioreps、
合法十列预置峰，运行原 Snakefile、原 helper 和原消费者。修复前 **2 failed / 8 passed**：
ChIP-seq/ATAC narrow 的原脚本均在 argparse 报缺少 final-output 值，脚本退出2，
Snakemake退出1；同一组产物/summary断言修复后 **10 passed**。
验证空 final_output 保持空、final_method 保持 none；CUT&Tag narrow 的非空路径
及实际值含空格的路径均完整进入 summary，峰内容逐字段一致，未要求生成该 final 元数据路径。
非法三列 narrowPeak 仍被原脚本及原规则拒绝，不放宽输入；合法空峰由原脚本测试保留。

新增测试还区分父级 reproducibility 省略/关闭时默认目标不调度 consensus，
与显式请求产物仍可进入该规则；覆盖 replicate_analysis=false。原有配置、
consensus/IDR DAG、脚本、catalog、身份及 PR-8/PR-9 DAG 定向回归另为
**164 passed / 0 failed / 0 skipped**。本轮独立有效用例合计174，未重复累计复跑。
九项无科学工具依赖的用例进入既有 fast 分片；失败任务诊断需要 conda 的一项负对照
进入既有 real_execution 层，CI 显式收集该文件，工具缺失即失败，不增加 skip。
精简 PATH 探针曾在九项完成后停于 Snakemake 的失败诊断；已停止，归为环境边界，
不记为 S2 缺陷或通过。未修改 Snakemake/conda 或安装依赖。

Ruff、限定文件 diff 检查通过；隔离 Snakemake lint 与原65条警告基线完全一致，
工具退出1（已有警告），未修改基线。snakefmt 因本地缺 shfmt 退出123，未宣称格式检查通过。
测试只执行 Python stdlib consensus 消费者，未激活/构建 conda 环境；不验证 MACS3/SEACR
生产者、比对、完整测序流程或部署环境。预置峰不冒充上游科学产物验证。

PR-10 基线在 PR-9 Usage 检查完成之后建立；四路径独立 diff、原始输入和命令证据在
`/tmp/helix-pr9-usage-pr10-lbwpfx6w/pr10/`，PR-9 文案 diff 不混入其中。
consensus.smk 字节改变 ENCODE 动态身份输入，前后字节/SHA已保留；完整入口涉及保护
profile，未捕获完整摘要。Bulk111项仍匹配原清单，本轮改动不在该闭包，亦未影响
Bulk 配置、参考或artifact/QC契约；不新增 Bulk Gate、不重建 manifest/qualification。
PR-3/4/6/7各身份与G0欠验分别保留，不跨身份追认。
独立复验见 `/tmp/helix-pr10-independent-review-yxbnu5zj/report.md` 与同目录
`commands.md`：164项定向回归、9项 fast 原规则用例、1项真实工具层负对照，
共174 passed；环境、工具和科学执行边界仍以该报告为准，不把本地通过写成远端 CI 已通过。

**PR-11 状态（2026-09-23）：S4/S13 实现与定向验证已通过独立复验，用户已接受。**
`workflow/rules/qc.smk:782` 只移除固定 `-speak=0`，保持 final BAM、
`-x=-500:15`、`-rf`、线程及失败处理。实际 `run_spp.R` 先计算曲线并按相关系数排序，
再由 speak 覆盖选峰；NSC/RSC 使用第一候选，不能称为“停止整条曲线计算”。
`scripts/parse_cross_correlation.py:_primary_fragment_length` 为表头/无表头
estFragLen 增加专用解析：保留工具序列的第一项；任一空、NA、非数值或非有限 token
使该字段为 NA，不跳过坏首项。七列 summary、其他标量解析、comment fallback、
样本名及 basename 脱敏不变，quality_flag 仍仅依赖 NSC/RSC。完整候选保留在原 `.cc.qc`。

证据目录 `/tmp/helix-pr11-cc-h_oa_bgg/`，入口为 `report.md`、`commands.md`、
`this-round.diff`。固定随机种子的合法 SE 合成 BAM 分别含20000/60000条50M记录。
同 BAM 同窗口的串行单变量对照仅增减 speak：主候选0→150；另一真实 R 输出
`200,300,1115`，旧 parser 为NA，修复后为200.0。具体 NSC/RSC 是本夹具的观测，
不固定为通用阈值或改善保证。真实单峰、多峰原规则和原 summary 消费者均保留 `-p=1`：
修复前 **3 failed / 1 passed**，相同断言修复后 **4 passed**。
初次沙箱 socket 拒绝已单独记录并停止；允许本地 socket 后完成原规则验证，未修改 R 或线程。

parser新增列表/CLI回归修复前 **10 failed / 49 passed**，修复后 **59 passed**。
相关 parser、QC配置、optional QC 与 summary/MultiQC 接线，以及PR-8～PR-10必要DAG回归
合计 **199 passed**；加上真实工具层4项，本轮独立有效用例共 **203 passed / 0 failed / 0 skipped**。
默认 cross_correlation 省略/false 均不调度，true只为treatment调度；显式目标可达性单独覆盖。
CI沿用现有测试分层，新真实工具测试由 `test/real_execution/` 收集；科学job新增已有
`phantompeakqualtools.lock` 环境，显式 RUN_SPP，fixture将其bin置于PATH以选择同环境Rscript。
未运行远端CI或重新安装锁定环境。本机phantompeakqualtools1.2.2、spp1.16.0、
caTools1.18.3、snow0.4.4；R4.4.3的本机构建与当前lock不同，具体build、摘要和模块路径已记录，
不能把本地通过当作完整锁定环境复验。

Ruff检查、格式检查及限定文件diff检查通过；隔离lint保留原65条警告，正文与定位完全一致，
不更新基线。初次关闭MultiQC的探针只产生64条，改用启用MultiQC的隔离配置后比较通过。
snakefmt因缺shfmt退出123，该项格式工具检查未完成。
本轮只验证小型输入的选峰、解析、原规则执行及DAG/接线；未执行比对、完整测序、
实际MultiQC HTML生成、conda部署或生物学代表性验证。S14窗口政策、fastp延期与后续顺序不变。
qc.smk和被引用parser的字节改变ENCODE动态身份输入；受保护profile使完整身份入口未运行。
Bulk111项闭包及其配置、参考、artifact/QC契约未受影响，不重建manifest/qualification，
不新增Bulk Gate要求；各历史身份及G0欠验仍分别保留。
独立复验 `/tmp/helix-pr11-independent-review-y0eknqei/report.md`、同目录 `commands.md`
重新执行199项回归＋4项真实工具测试，全部通过；用户已接受此范围。
完整锁定R环境、远端CI、snakefmt及完整测序的原有验证边界不因验收而消除。

**PR-12 状态（2026-09-23）：S16 限定依赖修复与定向验证已通过独立复验，用户已接受。**
`workflow/rules/common.smk:_bamcoverage_inputs` 沿用现有ChIP-seq的MACS3峰目录依赖，
新增条件仅为treatment、SE、CUT&Tag broad、extend_reads=auto/yes（省略默认auto）。
MACS3在 `peaks.smk` 声明的是峰目录output，日志仅为log；没有将日志单独塞入input。
保留 `metadata.smk` → `cuttag.smk` → `chipseq.smk` 的真实延伸委派、日志提取及原bamcoverage。
PE、no、固定整数、control和narrow CUT&Tag未新增峰依赖；narrow原有
`--nomodel --shift -100 --extsize 200` 保留；ChIP-seq已有行为和ATAC/MNase对照不变。

正式测试 `test/workflow/test_bamcoverage_dependencies.py` 运行原Snakefile、原MACS3规则、
原消费者及真实dispatch/helper；仅macs3、bamCoverage可执行文件为替身。
24项相同断言：修复前 **9 failed / 15 passed**，修复后 **24 passed**。
全新工作目录下、同一合法输入，旧CUT&Tag消费者优先argv为200、生产者优先为150；
修复后两种优先级都先完成生产者并使用150。单独请求bigWig也会拉入依赖。
无预测或预测小于60时，生产者完成后仍回退200并保留警告。
已记录150后仅改模型日志为250，下一轮检测params变化并重跑；再运行无变化时无任务。
这区分首次依赖与增量检测，不再声称“永不重建”。DAG初次求值仍可能打印既有回退警告，
测试以生产者完成后真实消费者argv判定，未修改警告策略。

其他配置/信号轨道及PR-8～PR-11必要DAG回归 **147 passed**，本轮独立有效用例
合计 **171 passed / 0 failed / 0 skipped**。CI既有deterministic分片自动收集新测试，无需改YAML；
本机显式使用已有conda支持Snakemake正常metadata/params检测，未激活或创建规则环境，
未运行远端CI。输入为可复现、真实samtools生成并quickcheck的微型SE/PE BAM/BAI；
调度测试的峰文件和`.bw`内容来自工具替身，只证明调度、argv和回退，不证明科学bigWig或MACS3模型。

证据 `/tmp/helix-pr12-dependency-sqepa9qp/` 的 `report.md`、`commands.md`、
`this-round.diff`、`logs/scheduling-results.json` 保留基线、全部红绿命令和原始输出。
Ruff及限定diff检查通过；隔离lint与原65条警告（正文及定位）完全一致，无需修改基线。
snakefmt仍缺shfmt，退出123，未完成该项格式工具检查。
规则字节改变ENCODE动态身份输入，完整捕获会触及保护profile，未执行、未声称获得完整摘要。
只读核对Bulk111项大小/SHA仍匹配，本次修改不影响该闭包及Bulk科学契约，不重建
manifest/qualification，不新增Bulk Gate。PR-3/4/6/7及G0各历史身份欠验分别原样保留。
独立复验 `/tmp/helix-pr12-independent-review-c1afxfc1/report.md`、同目录 `commands.md`
重跑 **171 passed**；用户已接受该范围。替身调度实验不等于真实科学bigWig或MACS3建模验证，
远端CI、snakefmt及各历史身份Gate边界继续保留。该轮未修fastp、S14、N2或MNase/Picard依赖。

**PR-13 状态（2026-09-23）：N1 限定判级修复与定向验证已通过独立复验，用户已接受。**
`scripts/chipseq_idr_summary.py:73–81` 与 `scripts/idr_reproducibility_summary.py:114–122`
改为在两个分母均为正时，用原始整数计数严格比较 `分子 < 2 * 分母`。
`compute_ratio` 仍返回原三位小数字符串、NA或inf；1.9996展示2.000但判pass，
恰好2、2.0002、NA、inf仍fail。只改判级代码块；计数、CLI、八列/十五列summary、
metadata及复制行为保持。ChIP-seq true→conservative、pooled→optimal，统一脚本复制true；
Nt>Np对照不改变这一N2政策，三个wrapper继续委派统一脚本。

正式回归扩展 `test/scripts/test_idr_reproducibility_summary.py`，复用原峰生成器：
唯一坐标合法十列narrowPeak覆盖两入口18组计数，各自真实CLI执行；三个wrapper另做代表性CLI。
相同测试字节红灯 **13 failed / 37 passed**，最小生产修复后 **50 passed**；
其中原11项测试及broad十七列产物复制用例保留。39组红绿CLI输入和复制产物逐字节相等，
summary仅13组status由fail变pass，其余字段及26组status不变。
既有配置、IDR路径及原规则DAG回归 **151 passed**，合计 **201 passed / 0 failed / 0 skipped**。
DAG的空FASTQ占位仅证明调度接线；没有运行IDR拟合、上游峰生产或完整测序。

单文件首次收集发现既有 `from scripts import ...` 依赖其他测试添加根路径，
已改为按当前checkout的明确脚本路径加载；canonical bootstrap/provenance仍启用，未改既有断言。
CI现有deterministic shard2收集全部50项，无需修改YAML，远端CI未执行。
Ruff check/format及限定四文件diff检查通过；没有改规则或lint基线，未重跑snakefmt。
证据及独立diff：`/tmp/helix-pr13-idr-bgzdq0j8/` 的 `report.md`、`commands.md`、
`this-round.diff`、`logs/cli-results.json`，保留首次收集错误与证据检查器重复统计pytest别名的纠正记录。

两个脚本均由原规则字面引用，字节改变ENCODE动态身份输入；完整入口会读取保护profile，
未运行、未声称获得完整摘要。只读核对Bulk111项大小/SHA及framed aggregate一致；
本次四路径不在闭包，未影响Bulk配置、参考或artifact/QC契约，不重建manifest/qualification，
不新增Bulk Gate。PR-3/4/6/7各身份及G0欠验分别保留，不跨身份追认。
用户确认独立重跑50项专项及151项相关回归通过；严格<2、N2复制政策及上述验证边界保留。
该轮未改变borderline、IDR科学参数或后续顺序。

**PR-14 状态（2026-09-23）：S17 限定依赖修复与定向验证已通过独立复验，用户已接受。**
`workflow/rules/mnase.smk:454–458` 在原 `mnase_qc_summary.input` 增加受规范化
`QC_CONFIG.picard_metrics` 控制的insert_size_metrics文件依赖，指向原Picard规则声明的输出。
Picard默认false；省略、false和字符串false不新增依赖，不增加关闭路径参考要求。
保持params、原summary脚本、23列表头及字段含义；该字段是路径或NA。关闭时若旧文件存在，
仍报告路径。依赖不受通用qc.summary开关控制，缺参考继续被既有校验拒绝。

正式测试 `test/real_execution/test_mnase_summary.py` 保留完整原Snakefile、Picard生产者、
MNase路径/范围/caller helpers及原summary消费者；仅外部Picard为替身，并在正常臂产生全部四个声明输出。
真实samtools 1.23.1生成、排序、索引和quickcheck合法小型PE BAM，原summary得到2/4/6条
read records（不是PE片段数）。预置bigWig是路径占位，不是科学产物。
同一组正式断言：修复前 **8 failed / 6 passed**，修复后 **14 passed**；
旧summary优先写NA、生产者优先写路径，修复后两种优先级均先完成Picard并写路径。
单独请求summary会拉入Picard；生产者失败或故意缺输出时不再抢跑summary；无变化再次请求无任务。
覆盖qc.summary=false、布尔/字符串开关、关闭但旧metrics存在、缺参考和mixed treatment/control默认目标资格。
红绿测试函数及断言AST相同；最终仅删除Ruff发现的未用顶层导入，并以最终文件重跑14项。

相关QC/MNase/参考配置、optional QC、MNase及PR-8项目DAG、路径catalog回归 **213 passed**，
独立有效通过合计 **227 passed / 0 failed / 0 skipped**，不累计重复绿灯。
新增测试进入现有scientific real-execution层，该层已准备锁定samtools并收集test/real_execution；
本机显式使用已有工具，未安装/物化规则环境，未运行远端CI或真实Picard分布计算。
Ruff及限定diff检查通过；隔离lint与既有65条警告正文/定位完全相同，基线未改。
snakefmt针对mnase.smk因缺shfmt退出123，该项未完成，不修改格式基线或安装共享依赖。
证据、命令、红绿产物、限定开工基线和独立diff：`/tmp/helix-pr14-mnase-u43pwmqv/`。

规则字节改变ENCODE动态身份输入，完整入口会触及保护profile，未捕获完整摘要。
只读核对Bulk111项大小/SHA和aggregate匹配，本次路径不在闭包，未改Bulk配置、参考或artifact/QC契约；
不重建manifest/qualification、不新增Bulk Gate。G0、PR-3/4/6/7等各历史身份欠验分别保留。
用户确认独立重跑14项专项及213项相关回归通过；真实Picard、科学bigWig、远端CI、snakefmt
及历史Gate边界继续保留。该轮未改变科学口径、S14、N2、fastp或后续顺序。

**PR-15 状态（2026-09-23）：P2 限定HTTP分页实现与定向独立复验已通过，用户已接受；新身份Gate仍待补验。**
`src/encode_pipeline/api/routes/runs.py:771、:824` 仅为events/logs的Query添加`le=100`，
保留default=50、ge=1及limit+1前瞻；服务/repository继续允许内部查询101条。
超限返回既有400/API_REQUEST_INVALID，不改envelope、认证、stream、游标或PR-16响应声明。
上限依据原run-history的1～100契约及前端PAGE_LIMIT=100，不扩展其他分页接口。

正式测试 `test/api/test_run_progress_pagination.py` 使用原应用、路由、RunService、
真实临时文件SQLite/内存repository和真实测试认证session；空registry隔离科学运行创建。
同一测试字节红灯 **24 failed / 46 passed**，修复后 **70 passed**。
SQLite旧版2^63−2返回200，2^63−1、2^63、2^100返回500；6次诊断另保存真实OverflowError异常链。
内存旧版超大正整数可返回200，不能混写为两个后端均溢出；新HTTP边界两后端一致拒绝。
覆盖省略、1/50/100、101及整数边界、非正数/非整数，99/100/101条分页及末页，
101条内部前瞻、顺序、cursor、stdout/stderr隔离、未知run/游标和未认证请求。
观察包装器委派真实查询，非法limit不进入分页查询；各fixture检查会话归还并释放engine。

API/history/OpenAPI回归41项、服务/双repository分页13项、SQLite重开及并发2项、
执行身份32项、前端资产3项通过，连同专项为 **161 passed / 0 failed / 0 skipped**；
Gate前置失败另列，不混入通过数。OpenAPI中会改HOME的既有测试按本轮纪律未选取，
会捕获受保护ENCODE profile的创建流程未运行；没有新增skip或修改原断言。
新增文件由现有deterministic shard1收集；本轮未运行远端CI。

执行现有`npm --prefix frontend run openapi:regenerate`，公开schema仅新增两个maximum=100，
生成客户端仅两个参数模型增加@maximum注释，number本身不执行范围校验。
前端runClient/RunProgressPanel **53 passed**，typecheck/build通过。通过正式资产工具同步
API摘要和asset manifest身份；正确cwd重建的8个静态文件与旧包逐字节相同。
Ruff check/format及限定diff检查通过。保留Orval import.meta/CJS、Vite大chunk等工具警告。
初次测试误设technical_message应省略（原契约为null）、Gate未指定marker被排除、
首轮Vite错误cwd导致Tailwind未命中，均已记录并纠正；未将这些搭建问题计作产品缺陷。

`api/routes/runs.py`在Bulk111项受控闭包内。本轮先保存111项及manifest/qualification精确副本，
再由正式工具生成；路径和文件数不变，仅该源码条目变化，persistence/schema projection不变。
旧aggregate `b6751715f1380ad5ae7ae993cf4a2c417602b24de81e37fdf6fed982492b3cf8`；
新aggregate `3535fdb864cf0fc5b32eb31c251d76e2282d136446b8812c20c88ecbe8160e6b`，
manifest SHA `bca78075eb62c32007c6ae098c14869cb1e2434d686888d4b70938ae738a3de7`。
本轮未改变ENCODE动态指纹选择的科学文件，也未运行会触及保护profile的完整捕获。
正式Protected Bulk Gate入口已显式选中并启用，在缺少
`HELIXWEAVE_BULK_RNASEQ_RUNTIME_ROOT`的前置检查退出1；fixture、测试Redis及受控Docker坐标
同样未提供，当前未提交工作区也没有对应的已审核clean SHA。未进入runtime admission、
Redis/RQ、科学/生命周期执行，生成qualification不是Gate通过。PR-15新身份欠验独立保留，
不追认G0、PR-3/4/6/7等历史身份。
证据、命令、旧新身份恢复材料及本轮独立diff：`/tmp/helix-pr15-pagination-t130y1uj/`。
用户确认该轮无需返修，授权继续PR-16；上述Gate欠验及边界不因推进而消除。
fastp延期和PR-2～PR-17顺序不变。

**PR-16 状态（2026-09-23）：P3 的500响应声明、实现与定向独立复验已接受；新身份Gate待补验。**
逐项核对30个operation，在auth、agent、preflight、runs、workflows五个路由文件中
补22处缺失的500声明；保留已有8处专用模型。createRun声明RunResponse；getRun声明
RunResponse与ValidationResponse联合，分别对应引用证据读取失败和全局异常。
未改main.py异常处理器、公共模型、成功响应、认证、生命周期或PR-15分页逻辑。
新增正式测试使用原应用、真实临时SQLite及认证session，窄依赖故障在写入/外部调用前阻断；
30个全局异常路径及2个显式引用错误路径均校验实际HTTP内容、脱敏和JSON Schema。
同一测试字节红灯24 failed/10 passed，修复后34 passed；32个故障响应前后逐项相同。
另覆盖14个正常读取响应。故障注入不代表正常请求必然失败，也不是科学执行验证。

PR-15分页70项、API/OpenAPI/auth/history/agent定向回归50项、执行身份32项、资产3项通过，
连同专项共189项Python通过；前端runClient/fetcher/RunProgressPanel 60项通过。
命令和原始结果见本轮证据，不将红灯、重复批次、Gate前置失败计入通过数。
构建内typecheck、Ruff check/format及限定diff检查通过；远端CI未运行。
既有成功流式下载测试在本机验证中未完成，单独记录，不修改原断言；其500准备失败路径已覆盖。
会改HOME的既有OpenAPI测试未选取。新增测试由既有deterministic shard1收集，无需改CI YAML。
正式openapi:regenerate只增加22个500响应，生成GetRun500联合类型及其导出，调用签名不变。
正确frontend cwd构建并经正式资产工具生成，8个静态文件逐字节不变；仅同步API摘要及资产身份。
API SHA为`a26661fd86823c8ac7fee5d62294bcab9d0df6100dc77ab6e8248dc16df82265`，
资产身份为`sha256-d69d10c27df371ff46ee04d6471e5d5a7eea6aea17e0c15be887045a086b3c24`。

Bulk仍为111项，仅preflight.py、runs.py、workflows.py三个受控条目改变。
开工前保存全部111项及manifest/qualification精确恢复副本，正式工具同步后逐条核对大小/SHA，
并独立重算aggregate；persistence/schema projection不变。
旧aggregate `3535fdb864cf0fc5b32eb31c251d76e2282d136446b8812c20c88ecbe8160e6b`；
新aggregate `1016090c76abc2012b1b7f875d9262851dbf973c0bdad8b50a47d7098a282dff`。
未改变ENCODE动态指纹选择的科学文件，未调用会读取保护profile的完整捕获。
正式Gate入口已启用并选中，在缺少HELIXWEAVE_BULK_RNASEQ_RUNTIME_ROOT前置检查退出1；
fixture、测试Redis及受控Docker坐标同样缺失，未进入admission、队列、生命周期或科学执行。
qualification生成不是Gate通过，PR-16欠验独立于PR-15、PR-3/4/6/7与G0，不跨身份追认。
本轮证据、operation对照、旧新身份、红绿日志和独立diff：`/tmp/helix-pr16-contract-31bkkhz8/`。
独立复验记录：`/tmp/helix-pr16-independent-review-0q29unx0/{report,commands}.md`。
主审重跑 **155 项 Python、60 项前端测试通过**，22 个新增 500 声明、getRun
联合模型、生成客户端及身份同步核对通过，用户已接受并允许继续 PR-17。
PR-16 当次复验中成功流式下载仍超时，原因当时未确定；没有证据认定为 PR-16 引入的回归。
该次接受不代表新身份 Gate、成功流式路径或远端 CI 已通过；后续调查及成功证据见下文收尾第三单元。

**PR-17 状态（2026-09-23）：主体及三处文案返修均已通过独立复验，用户允许继续当前批次收尾。**
仅修改 configuration、qc-interpretation、reproducibility-policy、defect-audit 与本 roadmap
五份 Markdown。按当前原函数/规则澄清 raw/blacklist-filtered 并存与消费者、
FRiP 的 read-record/未平移口径及配置依赖选择、已有条件 BigWig 转换、consensus
支持比例 score/max signal/p/q=-1 和一行 13 列 summary；说明连通分量而非逐碱基交集。
保留 PR-9/11 已修的复杂度输入与主候选解析，不修改科学政策。
审计保留 S1–S20、P1–P6、D1/D2 共 28 个编号；核对原正文 V/R 各 13，D1/D2
无标签、原表遗漏 D2。逐条区分原主张、复核纠正、已接受修复与待决/证据不足，
N1/N2/N3 单列，不计入原统计；旧建议顺序明确由本路线图取代。

本轮主要为源码与历史记录核对；通过 canonical bootstrap 后用原 consensus
脚本实际执行合法三 biorep narrowPeak/broadPeak 两组小输入，均退出 0，
输出 score=1000/667、最大 signal=30/40、p/q=-1、13 列汇总与 support_distribution
断言通过。不是测序、峰调用或 IDR 拟合验证；没有重跑历史回归/浏览器/完整 HTML。
本地链接、QC 11 项 7 true/4 false、原 28 条统计与限定文件 diff 检查通过。
在线核对 UCSC 格式与固定 ChIP v2.2.2 overlap/WDL；xcor wrapper 在线获取失败，
使用已归档版本源码并核 SHA，未声称补齐原先缺少的 overlap 归档。

开工字节/SHA、原始命令、独立 diff、旧新说法对照与检查结果在
`/tmp/helix-pr17-docs-mbv04ttu/`。仅文档字节变化；只读核对 Bulk 111 项、
manifest/qualification、相关 ENCODE 科学文件和现有 API/资产文件与开工基线一致。
未调用完整 ENCODE 身份捕获，未重建清单/资格/客户端/资产，不新增 Bulk Gate。
G0 与 PR-3/4/6 各身份、PR-7、PR-15、PR-16 欠验分别保留；PR-2 历史/部署边界、
PR-6 完整执行边界、fastp 延期及修复归属未定、科学政策和证据不足项均不关闭。
本轮停止待独立验收，不提交、推送或自动开始新任务。

后续独立复验记录：`/tmp/helix-pr17-independent-review-67O3U3Bg/`。
该次复验的两组原 consensus CLI、链接、统计和限定身份核对通过。按其意见，
本次仅修改审计 P4、S12、S15 正文及对应总表的历史转述：分别恢复为本次取消与
既有终态的区分、p/q=-1 的下游用途及全支持例、线粒体指标影响与先定政策。
保留原处置、28 条编号与 V/R 归属及全部历史欠验；本次只做文档核对、链接检查
及两份授权文档的限定 diff 检查，未重跑数值测试或执行 Gate。
返修基线、独立 diff、命令及报告：`/tmp/helix-pr17-wording-_16l_4sz/`。
返修已通过独立复验；用户已授权下文当前批次收尾的 fastp 单元。
下载超时调查、Hi-TrAC 实施及历史 Gate 欠验仍按各自边界推进。

| 顺序 | 范围与交付结果 | 验收要点 |
| :--- | :--- | :--- |
| PR-2 | MACS3 环境调查（bug #3）；区分继承 worker 环境的 legacy 路径与显式启用 conda 的托管路径。确认缺陷后只修实际出错的调用点。 | 记录实际命令、conda 激活、PATH、entrypoint/shebang 和解释器；相关 DAG 与受控小样本执行通过。不能复现时记录原失败身份或日志缺失等边界，不能把 dry-run 当作执行成功。 |
| PR-3 | 内存与 SQLAlchemy run-repository 使用同一套参数化行为测试。 | 生命周期、事件顺序、重复操作、冲突处理通过相同断言并进入 CI；发现生产差异先报告，不顺带重构。 |
| PR-4 | 通知与 doctor 静默失败的结构化日志；分别覆盖 `_record_outcome` 内层事件记录失败、外层 notifier 异常、doctor 静默路径。 | 故障注入能定位失败环节；不改变生命周期和失败处理语义。日志不得包含载荷、环境值、私有路径或原始异常文本。普通发送失败已有持久事件，不能再描述为所有通知失败无痕。 |
| PR-5 | 上游耦合台账，仅文档：MultiQC 样本名处理、nf-core 参数白名单、QC 输出格式等。 | 每项注明上游版本、实际源码/调用关系、本地实现与测试、升级复核事项；存在 wrapper 不等于实际 WDL 接线。 |
| PR-6 | Docker staging/admission 存储语义调查（bug #5）；先核 config digest、archive/index digest、加载后 image ID 的关系，提出单 daemon 方案后再实施。 | 同一 daemon 完成 staging → admission → 受控小样本；保留镜像完整性检查。若改变执行闭包，由实施方更新相应身份清单并完成适用的 Protected Bulk Gate，不靠放宽校验通过。 |
| PR-7 | 跨 section 冲突的表单内联提示，覆盖 QC/UMI 关闭后 advanced/outputs 依赖冲突。 | 优先消费既有校验结果，保留用户配置及后端拒绝规则，不在通用前端硬编码 bulk 语义；单元和桌面/移动浏览器验证通过。需扩展 schema 契约时先报告方案。 |
| PR-8 | S1：修复混合 treatment 项目默认 summary 将 MNase 拉入 MACS3 的依赖路径。 | 合法输入下覆盖 mixed 默认、summary-off、MNase-only；缺 FASTQ 时的 MissingInput 单独记录，不与 MACS3 guard 混淆。 |
| PR-9 | S6/S7：复杂度 QC 使用保留重复的适当输入阶段，并正确处理 preseq 的 PE 模式。 | 真实 samtools 去重前后与合法 PE 加减 `-P` 对照；允许报告 preseq 无法外推，不能预设总会输出低曲线。`distinct == total` 本身不是错误，去重后也非恒为 1/1/NA。S7 为 opt-in。 |
| PR-10 | S2：修复 opt-in narrow consensus 的空 `--final-output` 参数渲染。 | 主开关默认 false 的对照保留；开启后的受影响 ChIP/ATAC narrow 原规则实际执行通过，不仅是 dry-run。 |
| PR-11 | S4/S13：修复固定 `-speak=0` 覆盖选峰及真实多候选 estFragLen 解析。 | 同 BAM 单变量对照，真实 R 单峰/多峰输出驱动原 parser；不声称曲线停止计算或移除 speak 后必然多峰。S14 窗口政策不纳入本 PR。 |
| PR-12 | S16：补 SE CUT&Tag 自动片段延伸所需的首次执行依赖。 | 原消费者/helper 下先后调度得到一致 argv；单独覆盖下一轮 params 变化检测，撤回“永不重建”的主张。 |
| PR-13 | N1：两个 IDR summary 使用未格式化比率判级，展示精度与判级分离。 | 唯一合法十列峰构造 1.9996，应按现有 `<2` 契约通过；恰好 2 等边界保持原契约。不夹带 borderline、等号或 optimal 集合选择政策，不重复立项。 |
| PR-14 | S17：启用 Picard 时补齐 MNase summary 的生产者依赖。 | 原规则/helper 验证首轮顺序不再造成 NA/路径差异；`insert_size_metrics` 按文件路径字段验收，不当成数值。Picard 默认关闭。 |
| PR-15 | P2：events/logs 分页上限及 `limit+1` 边界处理。 | 使用真实临时 SQLite/HTTP 覆盖正常值、`2^63−1` 及更大输入，避免合法通过参数校验后因整数绑定溢出返回 500。 |
| PR-16 | P3：补齐实际受影响 operation 的 500 错误响应契约。 | 故障注入对照响应 envelope 与 OpenAPI；按既有流程再生成并检查客户端，保留脱敏。不能从故障注入推出正常请求必然失败。 |
| PR-17 | 文档同步：S5/S8/D2 的 raw/filtered 消费者、D1 过时文字、N3 consensus 字段、启用条件及审计措辞/统计。 | 说明过滤 BAM 与过滤峰的条件；peak_counts 只消费峰。score 为支持比例，p/q=-1 合法，summary 不提供声称的逐峰字段。修正文档不能暗中改变科学输入或输出契约。 |

### 证据与决策边界

- 后续验收采用补强后的合法 PE、真实 R 多峰、唯一合法十列峰及原
  helper 调度证据。分别标注真实工具/原脚本执行、隔离调度、仅源码推断；
  假生产者和假 bamCoverage 只能证明调度或 argv，不能证明科学产物正确。
- QC 默认真源为 `src/encode_pipeline/config/qc.py`：11 项，7 true、4 false。
  reproducibility、IDR、MultiQC 等另追各自配置解析，不能由 QC 子开关
  类推主开关。S2/S3/S9/S10/S11/S12 不应统称默认执行。
- 原上游证据归档缺 overlap 源码；引用与归档完整性应分别核对。
  工作区 TSV 存在不能证明平台列表、下载或浏览器可达。
- N2 为兼容政策：`docs/idr-contract.md:50–51` 明确约定
  true → conservative、pooled → optimal。是否改为上游较大集合选择
  需要单独科学决策，不能随 PR-13 舍入修复改变。
- S3/S9/S10/S11/S14 的算法或 ENCODE 兼容选择，以及 S19/S20 的额外
  汇总展示，保留决策边界；与上游不同不自动构成需要修复的缺陷。
  P1 不可信输入越界、S15 科学损害、S19 完整 HTML 可见性继续标为证据不足。
- S12、P4、P6 的原实现缺陷指控撤回；S18“完全没有分布生产路径”的
  中心指控撤回。原审计正文为 28 条，V/R 各 13，D1/D2 正文无标签；
  汇总表漏 D2。V 包含仅源码阅读，不等于实际复现；部分确证不能全部计作
  正确或误报。

### G0 与 Protected Bulk Gate

G0 独立跟踪此前维护批次欠缺的 Protected Bulk Gate，不计新功能 PR。
验收须绑定待验证的旧 revision、执行身份和当时的闭包；未来构建通过
不能追认旧身份。缺少旧身份或环境证据时明确记录缺口，不宣布已补验。

每个实施 PR 按最终 diff 对 pinned runtime、科学配置、参考身份、执行闭包
及 artifact/QC 契约的影响，依据 `AGENTS.md` 评估相应 Gate。ENCODE 自身的
科学改动不自动等同 Bulk Gate；共享平台/通知文件也不能按目录豁免。
仅同步文档且不改变闭包时不要求 Bulk Gate。只读身份核对不重建清单，
也不构成 Gate 通过证据。

## 当前批次收尾与 Hi-TrAC adapter 接入（2026-09-23）

用户已确认：先完成当前批次收尾，再实施独立 Hi-TrAC 预处理 adapter。
Hi-TrAC 的只读上游与契约调查可并行；不将新 adapter 的代码、运行时或身份变更
混入当前收尾版本。本节更新此前 fastp 延期安排，不重排 PR-2～PR-17，也不另编
维护 PR 编号。实施方逐个完成下列单元，每个单元停下等待独立验收；
本审核会话负责规划记录与独立复验，不承担生产代码实施。

### 当前批次收尾顺序

| 顺序 | 工作单元 | 退出条件 |
| :--- | :--- | :--- |
| 1 | PR-17 文案返修：仅修审计 P4/S12/S15 的历史转述及对应总表，同步验收状态。 | 对照原审计开工副本与复核证据；保持裁决、科学政策及编号不变，限定文档检查通过。PR-17 返修已通过独立复验，用户允许继续第二单元。 |
| 2 | fastp 固定版本字段兼容修复：核对上游 `summary.fastp_version`，修正原解析器及相应契约测试。 | 合法上游布局可提取 QC；核对共享 artifact 消费路径，保留版本、计数、重复键和非有限值拒绝；用相同断言保存修复前后证据。按实际闭包影响同步必要身份。332 项独立复验通过，用户已接受；完整 Gate 留待集成收尾。 |
| 3 | 成功流式下载超时定位：区分测试客户端、资源清理、运行环境与产品行为。 | 验证正常下载及相关审计路径；有证据的缺陷做最小修复，无缺陷则以原因及可复核成功证据收尾。不得靠放宽断言、跳过测试或单纯增加超时宣布通过。已独立复验并获用户接受：当前受限执行环境的 socket 唤醒问题，无需生产/测试修复；原测试及真实下载清理通过。 |
| 4 | 最终集成验证与完整 Gate：整理累计改动和精确身份，补齐真实运行前提。 | 最终集成准备、最小格式收尾及 90 路径候选已通过独立复验；已形成精确 clean commit `842251e48408dcdaa8b3054403be5e2101951f0e`。该提交在专用 rootful 环境的原 Python Gate 4 项、浏览器产品链 5 项全部通过，且已独立验收。该 SHA 的远端 Protected Bulk attempt 2 已通过独立复验，用户接受，详见后文；PR-6 在线 staging 完整链、历史格式债及各历史身份当时未验证的边界继续保留。 |
| 5 | 推送与交付审阅。 | 用户明确确认目标及上传可见性后，上述精确 SHA 已推到原 GitHub 仓库的 `maintenance/20260924-final-gate` 分支；远端引用复核一致，main 仍为父提交。普通 hosted 作业、Lint 和 Lock Check 均成功；Protected 首轮因临时容量监护停止 Docker 而失败。修正并验证临时监护后，仅重试该失败 job 一次，attempt 2 已成功：Python 4项及浏览器5项通过、零失败零跳过；本轮收尾已通过独立复验，用户接受，详见后文。未重跑普通 hosted 作业，未重复 dispatch、合并或发布。本次新增 roadmap 规划和运行记录仅留在源工作区，不属于该已验收、已推送提交。 |

Gate 准备包括受控 runtime root、fixture manifest、测试 Redis、Docker CLI/socket，
以及现有准入契约要求的工具与隔离能力。不得绕过镜像、身份或环境校验。
完整 Gate 绑定已审核、干净的精确 commit；为 Gate 准备该 commit 时仍遵守提交授权。
一次最终版本的完整 Gate 只证明最终身份，G0 与各历史身份继续保留未验证记录，
不以新版本通过追认旧版本。普通 push 不会自动触发 Protected Bulk Gate；
若采用远端 `workflow_dispatch`，须先经授权使候选 SHA 可检出，再显式启用该 Gate。
本地通过与远端 CI 通过分别报告；若 Gate 后实现或执行身份变化，重新评估对应验证。

S3/S9/S10/S11/S14/N2 科学政策、P1/S15/S19 证据不足项、S19/S20 额外展示选择，
以及 bug #3 的历史与完整部署边界保留，不为清零而改代码。不得将部分产物存在
当作平台可下载或完整报告可见的证据。此前临时证据缺失时注明缺口，不重构造历史结果。

### fastp 字段兼容修复状态（2026-09-23）

**实现与定向独立复验已通过，用户接受。** PR-17 三处文案返修已获独立复验通过。
核对 fastp 1.0.1 JsonReporter 及 nf-core/rnaseq 3.26.0 固定源码传递路径后，
共享 `parse_fastp_summary` 先验证 summary 对象，再读取 `summary.fastp_version`；
不增加顶层 fallback。版本、计数、JSON 限额、阈值与聚合约定保持；只修正两份测试
文件的版本布局并补回归，没有改 QC/artifact 调用者、其他解析器或运行时依赖。

固定测试字节：修复前 **21 failed / 34 passed**，同一组修复后 **55 passed**；
QC、artifact 两份完整测试与执行身份回归 **332 passed**，无 skipped。
Ruff、格式及限定文件 diff 检查通过。公开 QC 和 artifact 消费链均实际执行；
输入为上游形状的合成 JSON、临时状态表/占位产物，不是 fastp 二进制或 nf-core 流程实测。
所查 PATH 和已有环境无 fastp；没有安装工具，部署版本及真实报告仍未验证。
NF 在线获取失败，已有两份源码归档重新与固定 source manifest 核对大小/SHA；
fastp 固定 tag 通过官方源码在线核对。详细测试与消费者索引见台账 UC-06/UC-07。

Bulk 闭包仍为 111 项，仅 `adapters/bulk_rnaseq/status_evidence.py` 字节改变；
开工前保留全部受控源码及 manifest/qualification 精确副本与 SHA。
正式生成器更新 manifest/qualification，persistence/schema projection 不变。
旧 aggregate：`1016090c76abc2012b1b7f875d9262851dbf973c0bdad8b50a47d7098a282dff`；
新 aggregate：`12fa7ec0d8c46d8689ad3042636b677468dca1957765eca2425c1193c174d62a`。
本单元不启动已知缺运行坐标的 Gate；完整 Gate 留在第四单元集成收尾。
资格记录生成及定向身份测试不是 Gate 通过；本新身份待验证，G0 与各历史身份欠验独立保留。
未运行完整 ENCODE 身份捕获，未打包、提交、推送。

证据、开工恢复材料、独立 diff、完整命令与原始日志：
`/tmp/helix-fastp-version-s51wlpg9/{report,commands}.md`、`baseline/`、`this-round.diff`。
独立复验 `/tmp/helix-fastp-independent-review-8q60pff8/report.md` 亲跑 332 项通过；
用户已允许进入第三单元。fastp 新身份完整 Gate、真实工具/部署边界继续保留。

### 成功流式下载超时调查状态（2026-09-23）

**调查及定向验证已通过独立复验，用户已接受；分类：环境问题，无需生产或测试代码修改。**
旧 PR-16 超时日志仍可读取；本次先在受限执行环境原样单独运行
`test_artifact_download_is_audited_with_actor_and_target`，35 秒未完成、有界终止（exit124）。
原测试及 fixture 的字节未变。阶段诊断显示认证、prepare、审计及响应头已完成，
工作线程已生成首块，事件循环尚未接收结果；后台清理和 fixture teardown 尚未进入。

不加载产品代码的 AnyIO/Starlette 最小对照同样阻塞：Python 3.12.13 的
`asyncio/selector_events.py:152` 写 self-pipe socket 得到 `EPERM`，工作线程回调
已在事件循环 ready 队列中，无法唤醒 selector。没有修改网络配置或线程池。
在允许本机 socket 通信的执行环境中，用相同解释器、依赖、测试源码、断言及
35 秒上限重跑原测试，约 4 秒正常退出；最小 AnyIO 对照也正常完成。
因此不将本次超时归因于 PR-16 的 500 声明修改，也不从旧栈单独推断所有历史超时的原因。

原审计测试使用替身下载计划，只证明路由与 actor/target 审计。
补充两组原应用、真实临时 SQLite、合成认证及原下载服务的对照：ASGI 和仅绑定
`127.0.0.1` 的临时 Uvicorn HTTP。270004 字节合成产物按原 descriptor revision
与路径身份登记，每组连续下载两次，均按 5 个数据块返回、完整字节/SHA 与输入一致；
状态 200、下载头、审计每请求一条、401/409 拒绝及拒绝不新增下载审计均通过。
每个下载计划的 9 个描述符关闭，连接归还，原 lifespan 关闭队列/持久化，临时源文件可删除；
服务器正常退出，仅余主线程。观察器仅委派原 prepare、记录计划与 ASGI 消息，不替换文件读取。

原测试单独 **1 passed**；与相关认证/下载路由/下载服务组合 **77 passed**，包含原测试，
不重复累计数量；既有迭代中断、读错、路径替换、cleanup 等断言保留。
没有证据提示顺序污染，未扩大到整套 API。真实网络客户端中途断开未新增实测；
既有路由中断用例仍是替身发送异常，服务层停止迭代用例仍按其原边界解释。
本机正常下载不是完整部署、代理、大文件/并发负载或科学工作流验证。

本轮仅修改本 roadmap；限定基线核对代码、测试及 Bulk 111 项、manifest/qualification
字节不变，保留 fastp aggregate
`12fa7ec0d8c46d8689ad3042636b677468dca1957765eca2425c1193c174d62a`。
未重建身份、打包或运行 Gate；完整 Gate 留待第四单元，所有历史身份与 G0 欠验分别保留。
Ruff/格式仅因没有 Python 修改而不适用；对本 roadmap 执行限定 diff 检查。
证据、版本、原始超时栈、命令、资源清理、基线与独立 diff：
`/tmp/helix-download-timeout-m2lzchf2/{report,commands}.md`、`logs/`、`baseline/`、`this-round.diff`。
本单元结束停下验收，不启动最终集成 Gate 或 Hi-TrAC，不提交或推送。
后续独立复验：`/tmp/helix-download-independent-review-76d4zlef/{report,commands}.md`。
“仅余主线程”限于原测试及真实下载完成清理；短诊断的 daemon 观察线程随其独立进程退出，
不是产品泄漏。旧报告的 77 项组合由复验核对日志，本次复验亲跑原测试与机制/真实下载探针；
不将旧日志核验写成独立重跑 77 项。用户已允许进入第四收尾单元。

### 最终集成回归与候选准备状态（2026-09-23）

**本地回归、前提核验和候选范围准备已通过独立复验；不是完整 Gate 通过。**
本轮只更新本 roadmap，未修改生产代码、测试、依赖、客户端、资产或身份文件。
开工 HEAD 为 `13d6c8ed26961150a6ed1bc93e445b3793125f43`，分支 main；
现工作区未提交内容按来源逐项核对，未创建临时 commit、暂存、推送或触发远端工作流。

亲跑当前候选字节：平台 577、API/下载/契约 298、Bulk adapter 792、科学规则及脚本
526 个独立 Python 用例最终通过，合计 **2193 passed**；前端完整单元集合
**392 passed**，桌面/移动原浏览器用例 **6 passed**。原成功下载测试及真实
ASGI/loopback 下载均在允许本机 socket 的环境通过，body、审计和资源释放保持。
最初 bootstrap 参数拒绝、科学测试的既有 namespace 收集顺序问题、canary 临时目录
重定向遗漏单列保留；未改源测试/断言，修正 harness 或采用既有收集顺序后通过，
复跑节点没有重复累计。含保护 profile 的完整身份/执行 fixture 未运行，不声称全仓 CI 通过。

前端基于当前 337 个文件的隔离副本，使用 Node 22.23.1/npm 10.9.8，typecheck、
构建及官方 OpenAPI 导出/客户端生成步骤通过；副本生成的 OpenAPI、客户端、全部
静态资产和 asset-manifest 与源工作区逐字节相同。临时资产输出未写回仓库或部署。
科学证据分层：samtools/preseq 15 项、真实 R 4 项；MNase 14 项保留真实 samtools、
原规则/helpers/summary，但 Picard 是替身、bigWig 预置；consensus 10 项预置合法峰、
执行原消费者。其余 DAG/调度/脚本测试按各自边界解释，不冒充完整测序或科学 bigWig 验证。

Ruff 对 47 个候选 Python 文件的检查和格式检查通过。补齐任务局部锁定 shfmt 后，
按配置的 snakefmt 13 文件通过；显式检查被配置排除的 qc/consensus 仍失败，HEAD
对照也存在格式债。qc 新增 PE params 沿用现有格式，仍落入格式差异区域，不能笼统
声称新增行都满足 snakefmt；未改排除配置或批量格式化。Snakemake lint 的 65 条
正文和定位与现基线一致，原工具退出 1、比较器退出 0；限定文件 diff 检查通过。

候选清单区分 88 项整文件改动、bug-log 的 #3/#5 必要 hunk、混合 roadmap，
以及默认排除的两项前置文档和五项保护内容。bug-log #6、roadmap 的前置 Tier 1
背景另列待审阅，不把所有当前差异自动纳入。用户已有 Hi-TrAC 计划保留，未开始实施。

只读核对 Bulk 111 项、manifest/qualification 及 aggregate 均匹配 fastp 后身份
`12fa7ec0d8c46d8689ad3042636b677468dca1957765eca2425c1193c174d62a`；
未重建身份。ENCODE 完整动态摘要因保护 profile 未捕获。五项 Gate 坐标在当前任务
环境均未设置；实际 runtime/reference/fixture、受控 daemon、真实 Redis 和 protected
runner/offline 缓存仍待核实。本机 unshare 小命令和浏览器成功不替代这些前提。
PR-6 未补跑在线 staging、完整 admission 或科学小样本；没有连接既有 Docker daemon。
完整 Gate 须先审核候选范围，再明确授权形成干净精确 commit、远端可检出 SHA 和
显式 workflow_dispatch；普通 push 不触发 Gate。本轮不越过此授权边界。
G0 与各历史身份欠验分别保留，未来最终身份通过不能追认旧身份。

报告、实际命令、JUnit/原始输出、候选文件/hunk、SHA、Gate 就绪表及本轮独立 diff：
`/tmp/helix-final-integration-9qause_r/` 的 `report.md`、`commands.md`、
`regression-matrix.md`、`candidate/candidate-scope.md`、`identity-inventory.json`、
`gate-readiness.md`、`this-round.diff`。本轮临时服务已结束；保留科学政策、证据不足项
和部署边界，停止等待独立验收，不开始下一单元。

### 当前批次最小收尾（2026-09-23）

**限定格式整理、已选范围的候选包与定向验证完成，待独立验收。**
前一集成单元的独立复验见 `/tmp/helix-final-independent-review-p6zmu82n/report.md`：
主审亲跑 Bulk 792 项通过；2193 个 Python 唯一通过节点、前端 392 项和浏览器 6 项
为核验原日志，不改写为本轮全部重跑。用户接受集成准备，并已裁定本轮范围。

唯一生产规则变化为 `qc.smk::preseq_complexity` 的新增 params 块：
`paired=` 去掉等号两侧空格，将两行移到 conda 后；lambda、layout 判定、输入、
输出、shell 与环境均保持。前后临时副本经 snakefmt 2.0.3 和锁定 shfmt 格式化后
字节完全相同；当前 paired 块与 formatter 输出一致。全 qc 文件显式检查前后仍
退出 1，剩余历史格式债、consensus 和原排除配置未修改。未批量格式化。

亲跑原 complexity DAG/argv **9 passed**，显式 real_execution 层原 samtools/preseq
小测试 **15 passed**，共 **24 passed、0 failed、0 skipped**；默认收集中的15个
real用例 deselected 已由第二命令全部执行。保留 SE/PE、去重策略、MAPQ、数值、
实际原规则执行及合法外推失败边界，不把它们称作完整测序。没有重跑2193项或前端。

候选采用已选的88项整文件改动（qc纳入本次格式调整）、bug-log仅#3/#5，以及
**排除前置Tier 1背景**的roadmap。两份前置文档、bug-log #6、roadmap前置背景和
五项保护内容不纳入；排除只作用于 `/tmp` 候选包，源工作区用户内容未回退。
按当前字节重新生成patch，在仅含许可路径的HEAD副本用原生git apply核验，包含
新增、删除及资产；应用后的所有生产/测试/资产字节与本轮源工作区一致，两个文档
的有意排除另列。没有临时index、commit或暂存。旧minified JS的原生Git重放已由
独立复验通过，前次仅是自制重放器限制；科学测试曾因原收集顺序出现两项导入错误，
catalog先行后通过，没有亲跑旧HEAD同命令，不能声称动态证明旧HEAD也失败。

只读Bulk111项、manifest/qualification与aggregate继续匹配
`12fa7ec0d8c46d8689ad3042636b677468dca1957765eca2425c1193c174d62a`，
未重建身份。qc规则字节改变ENCODE动态身份输入；完整摘要仍因保护profile未捕获。
完整Gate、PR-6完整链、当前身份及G0/各历史身份欠验保留。后续需要单独授权候选
commit、远端可检出和workflow_dispatch，并提供或批准新建隔离runtime、fixture、
Redis及单Docker daemon；本轮未部署或运行Gate，未开始Hi-TrAC。

本轮基线、独立diff、候选文件/hunk及patch SHA、原生命令重放、原始日志和后续
隔离Gate方案：`/tmp/helix-final-minimal-fyq3vqe2/` 的 `report.md`、`commands.md`、
`this-round.diff`、`candidate/manifest.json`、`candidate/candidate.patch`、
`gate-next-steps.md`。限定diff检查通过；不提交、推送，停止等待独立验收。

### Gate 环境空间收尾（2026-09-23）

**限定清理已通过独立验收；该清理单元未部署新环境或运行 Gate。** 最小格式收尾及
90 路径候选范围已通过独立复验。按用户本轮授权，核对指定 demo SQLite 中
`bb1cf001-7f63-4f94-b4ec-17750e7fcda7` 已 succeeded、45 个任务退出 0、无重新
入队请求，以及宿主进程无目录占用后，仅删除该 run 的 `engine/work` 内 102 个
大型暂存/派生中间文件；保留 1,035 个普通文件、符号链接及所指文件、诊断和命令。
results/reports/logs、配置与身份、Nextflow cache 未删除；54 个登记产物和 9 个
外部输入/资源的 SHA 前后一致。既有 QC 索引失败及证据保留；被删中间文件不再
可作 Nextflow resume 缓存命中，不将此次清理称作保留完整 work 恢复副本。

两个明确授权的 packaging 测试副本另删 29,988 个普通文件；保留原失败场景、
当时构建归档、外层 before/XML/日志/哈希、全部链接及保护同名条目。合计分配块
约 154.77 GiB，清理后文件系统可用约 187.13 GiB。只读夹具目录的权限失败及
余项脚本日志名冲突均有记录；经限定权限执行与纠正后，精确清单全部完成。
未清理下载缓存、已安装环境、浏览器、参考或运行资产，未停止既有服务。

当前固定 runtime 的原静态闭包验证、指定 daemon 的 34 个镜像身份/RootFS 与
原只读 Docker 可用性检查通过；不是 canary、完整 admission 或科学/Gate 验收。
可按复用已校验 runtime、另建隔离 daemon 的方案规划新增峰值 135 GiB，并保留
至少 50 GiB 文件系统余量；这是准备预算，不是实测最低需求或通过保证。
新根拟为 `/home/yangzichen/helixweave-gate/20260923`，尚未创建；服务、工具和
下载仍待授权。当前方案暂不需要新盘，完整复制 runtime 的更宽方案仍需额外容量。
完整 Gate、PR-6 完整链、G0 及各历史身份欠验不变；不启动 Hi-TrAC。

本轮清单、清理前后检查、实际命令、权限/探针失败及容量方案见
`/tmp/helix-gate-space-cleanup-927jgdy7/{report,commands,capacity-plan}.md`。
仓库仅本 roadmap 更新；旧候选 patch 未覆盖，该文档新增记录须在后续候选包
重新核对，不能继续以旧 patch SHA 代表含本轮文档的候选版本。

### regulation_of_t 限定清理（2026-09-23）

**116 项限定清理已通过独立验收；该清理单元未准备环境或执行 Gate。** 用户已确认
该项目用于测试学习，只需现有峰与 CPM BigWig，接受以后重算 MACS 信号。
依照独立审查 `combined-candidates.json` 的 candidates 数组，重新核验全部
路径、父目录和文件身份、62 个原始 FASTQ、关联 37 个 BAM 实体，并用现有
samtools 1.23.1 对这 37 个 BAM 运行 quickcheck（退出 0，仅基本结构检查）。
宿主可见进程检查无项目/候选占用，精确执行计划先封存后逐文件 unlink。

实际删除 42 个非空裁剪 FASTQ、74 个 MACS pileup/control bedGraph，116 项
均成功；分配块合计 175,770,890,240 B（163.6994 GiB），文件系统净增约
163.6992 GiB，后置核验可用约 **350.82 GiB**。原始输入、所有 BAM/索引、峰、
BigWig、QC/日志/配置及记录保留；8,776 个保留条目核验通过，只有直接删除所在
目录的预期时间戳变化。缺 BAM 的 20 个裁剪 FASTQ 和 20 个 bedGraph、两个
SE 零字节占位均保留。保护文件未读取/哈希/改动；无目录、缓存或环境清理。
没有验证完整测序内容，也不承诺 MACS 自动重建或与历史信号逐字节一致。

按此次实际余量计算，后续新增 135 GiB 后约余 215.82 GiB，比至少保留 50 GiB
另有约 165.82 GiB 余量；仍只是准备预算，未实测 Gate 峰值或宣布 Gate 通过。
未创建 Docker/Redis、下载、挂盘、重建身份、打包、提交或推送；历史欠验保留。
证据、逐项操作日志、保留核验及独立 diff：
`/tmp/helix-regulation-cleanup-nxde8l41/{report,commands}.md`、`execution-plan.json`、
`logs/`、`this-round.diff`。旧审查证据与候选 patch 未覆盖；本轮新增 roadmap
记录须在后续候选包重新核对，不沿用旧 patch 摘要代表新文档。完成后停止验收。

### 本机隔离 Gate 环境准备（2026-09-24）

**核心环境准备及系统依赖补齐已通过独立复验；完整 Gate 尚未执行。**
用户确认前述两轮限定清理已通过独立验收；regulation_of_t 独立记录为
`/tmp/helix-regulation-cleanup-independent-bc04h2c2/report.md`，首轮容量复核为
`/tmp/helix-cleanup-independent-review-mj303_m9/runtime-capacity-review.md`。
本轮在新根 `/home/yangzichen/helixweave-gate/20260924` 准备隔离环境，
只读原地复用既有固定 runtime，没有复制整个 runtime、修改既有服务或删除用户数据。
开工实测可用约 350.71 GiB；授权新增峰值 135 GiB、至少保留 50 GiB、下载 20 GiB。
20260923 根与 187.13 GiB 余量仍仅是前一日期的提案/记录。

- 新 rootless Docker 29.1.3 / containerd 2.2.0 使用独立 data-root、exec/state、
  Unix socket 和 namespace，确认 `io.containerd.snapshotter.v1` 存储。
  34 个固定归档全部装载，逐项核对精确 archive target image ID 和完整 RootFS diff IDs，
  不用 tag、config digest 或 distribution digest 替代。原 `--phase verify` admission
  与离线 Java/Nextflow canary 在重启前后均通过；这是复用归档的组件验证，不是在线 staging。
- 原 tiny 生成器产生合成输入；新 daemon 中真实 STAR 2.7.11b、Salmon 1.10.3、
  SortMeRNA 4.3.7 分别完成索引。原 finalize 和 `load_acceptance_fixture` 均通过，
  manifest SHA 为 `986afba027f617b3049494e9df939ca6572208544b393372327f621282922462`。
  保留真实 argv、镜像/config 身份、35 个索引输出文件摘要及 provenance；不代表科学流程执行。
- 专用 Redis 7.0.15 仅监听新 Unix socket，凭据/坐标为 0600，原 API/worker client、
  认证、锁与 RQ JSON 入队/读取/清除通过；未执行任务。重启验证通过，测试键清零。
- ci-fast 175 项锁定依赖在真实断网 namespace 创建；用户补充授权后，对明确范围、
  逐字节核对且排除保护文件的 290 文件准备副本正常 editable 安装，pip check 与原
  canonical bootstrap verify-checkout 通过。未安装源工作区、共享环境或用户 site-packages；
  该副本不冒充经审核的 clean commit 或完整 ENCODE 发行副本。
  Node 22.23.1/npm 10.9.8、npm 离线安装、Playwright 1.61.1/Chromium 1228
  离线安装和双视口合成页面启动/关闭通过，不是产品浏览器 Gate。
- 环境准备交付时，原 CI Playwright 系统依赖 dry-run 退出 1；用户已授权精确
  16 项新装与 4 项 GLib 升级，20 个官方 deb 已下载并核验。当时 `sudo -n`
  退出 1、要求密码，因此该次非交互尝试未进入安装。后续交互安装、中断、用户恢复
  与独立复验分开记录如下；sudo 认证不再列为当前阻塞。未修改 sudo 配置或全局代理。

五项运行坐标已生成并完成各自组件检查，完整值保存新根 `private/gate.env` 与
`private/gate-coordinates.json`，报告脱敏。结束前新 Docker/containerd、Redis 和容器
均已停止，环境/镜像/fixture 保留；重启及精确停止方法随交付提供。
系统安装前的观察值为新增根分配块约 36.95 GiB、文件系统可用约 313.70 GiB；下载载荷
257,685,810 B（约 245.75 MiB，预算保守计 0.31 GiB）。逐步容量记录均低于授权范围，
这些是准备阶段观察，不是完整 Gate 峰值保证；本次文档收尾未重测安装后磁盘占用，
不将上述容量写成当前实测值。

Bulk 111 项、manifest、qualification 前后核对一致，aggregate 仍为
`12fa7ec0d8c46d8689ad3042636b677468dca1957765eca2425c1193c174d62a`，未重建身份。
未调用触及保护 profile 的完整 ENCODE 身份捕获。完整 Gate、PR-6 在线 staging →
admission → 科学小样本全链、clean 审核 commit 与远端 protected CI 均未完成；
G0 及各历史身份欠验独立保留，不跨身份追认。

报告、完整命令/原始日志、失败搭建记录、预算、停止核验及 roadmap 独立 diff：
`/tmp/helix-gate-env-p8n2d5y1/{report,commands}.md`、`this-round.diff`。
环境准备单元只改 roadmap；当时两轮清理及环境记录尚未纳入旧候选补丁。
核心组件和服务收尾的独立复验为
`/tmp/helix-gate-env-independent-5Rkg6fCG/report.md`，已获用户接受。

### 系统依赖恢复与最终候选收尾（2026-09-24）

- 首轮交互安装在等待安装子进程时被用户中断，出现 `KeyboardInterrupt`。
  安装器缓冲 stdout/stderr、子进程结束后才写安装日志；首轮缺安装日志不能证明
  未执行安装。随后两次重试均被未完成配置检查拒绝，不能记录为旧安装器完整成功。
- 用户手动配置已解包的 16 项及正常触发器，再按授权精确版本离线安装剩余四项
  `xfonts-cyrillic`、`xfonts-scalable`、`xserver-common`、`xvfb`；未扩大包清单。
  安装器的 hook 检查范围是已登记的 `DPkg::Post-Invoke` 列表，
  不表示它拒绝所有未知 APT hook；本轮不重跑或返修安装器。
- 审阅方随后只读核验：20 包均为授权版本及已安装状态，相对安装前包基线恰好
  20 项变化（16 新装、4 升级），无额外包变化，`dpkg --audit` 退出 0 且无输出。
  原 Playwright 依赖预检退出 0，输出 `All system dependencies are installed.`；
  断网 1280×800 / 390×844 合成页面探针通过，浏览器已退出。这关闭系统依赖缺口，
  不代表真实产品浏览器链或完整 Gate 通过。证据：
  `/tmp/helix-system-deps-postinstall-review-KXMOYSXM/report.md`、
  `packages-review.{md,json}` 与该目录命令、stdout/stderr。本单元核对这些旧记录，
  未重跑系统检查或浏览器实验；Docker/Redis 保持停止，本轮未启动服务。
- 最终候选范围仍为 88 项整文件变更及两份局部文档，共 90 路径；bug-log 仅纳入
  #3/#5，roadmap 只在候选副本中排除前置 Tier 1 背景，保留清理、环境、恢复记录
  和既定 Hi-TrAC 计划。源工作区被排除的用户内容原样保留。新候选补丁及逐项清单在
  `/tmp/helix-final-candidate-2ad2jsed/`，旧补丁未覆盖，旧 SHA 不代表新文档。
  仅在许可路径的临时副本做原生重放及字节/模式/删除状态核对；本轮只读核对
  Bulk 111 项及原 qualification，aggregate 仍为上节值，未重建身份。

**该候选准备单元结束时：待独立验收，尚未形成 commit。** 当轮证据为
`/tmp/helix-final-candidate-2ad2jsed/{report,commands}.md`、`this-round.diff`、
`candidate.patch`、`candidate-inventory.json` 及 `next-steps.md`。
之后须另获暂存/提交授权，在隔离范围形成经审核的精确 commit，并核验包含
untracked 的 clean 状态；完整 Gate、push、workflow_dispatch 分别按授权执行。
本轮未暂存、提交、推送、dispatch 或执行完整 Gate，未开始 Hi-TrAC。
本机环境准备、未来本地完整 Gate、远端 protected CI 分别验收；PR-6 在线 staging
与完整科学链仍未补验，不能用归档装载/admission 代替。G0 及各历史身份欠验保持独立。

### 精确提交与本机 Gate 验收（2026-09-24）

90 路径候选已形成提交 `842251e48408dcdaa8b3054403be5e2101951f0e`，tree 为
`bb075b7e76aea272f223cc390c2b793c493ce843`，父提交为
`13d6c8ed26961150a6ed1bc93e445b3793125f43`。隔离 checkout 含 untracked 的
clean 检查通过，Bulk 111 项 aggregate 仍为
`12fa7ec0d8c46d8689ad3042636b677468dca1957765eca2425c1193c174d62a`。
专用 rootful 环境中，原 Python Gate 4 项与原浏览器产品链 5 项均零失败、零跳过；
该结论已通过独立验收，代码、测试与身份未因本轮验收而修改。

实施证据在 `/tmp/helix-rootful-gate-rp0o022v/`，独立验收在
`/tmp/helix-rootful-gate-review-0yyxwrqc/report.md`。独立验收核对原始记录及身份，
没有重跑科学流程或浏览器。五项浏览器测试不等于五次科学执行；下载验证覆盖
已记录的 Salmon 元信息产物，不外推所有产物或完整 HTML 可见性。测试服务已停止。
此前 rootless/监护失败记录保留；当前成功不追认 G0 或各历史身份，PR-6 在线
staging 全链及远端 protected CI 仍未执行，ENCODE 完整动态摘要未新增捕获。

本次仅补充验收状态与 Hi-TrAC 规划，不改写或追加到上述已测试提交。
推送目标为 `https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline.git`
的 `maintenance/20260924-final-gate`；远端检查时 main 仍为父提交且目标分支不存在。
首次推送被自动审批在执行前拒绝；用户随后明确确认该目标及上传内容可见性。
重试推送退出 0，远端维护分支精确指向上述 SHA，main 仍为父提交。推送单元未合并、
发布或 dispatch。新增 roadmap 规划仅留在源工作区，未混入已推送提交；该单元证据在
`/tmp/helix-push-hitrac-plan-ra8cqokr/`。

### 远端 CI 与 Protected Bulk：接管、执行和结果（2026-09-24）

用户随后明确授权远端 CI。核对维护分支仍指向上述 SHA 后，分别 dispatch 了
CI、Lint、Lock Check；GitHub API 均返回 204。CI 的参数为
`bulk_rnaseq_real_execution=true`、`bulk_rnaseq_expected_sha=842251e48408dcdaa8b3054403be5e2101951f0e`。
三次 run 的 head SHA 均匹配；维护分支普通 push 本身不触发这三份 workflow。

- [CI 35918197084](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/35918197084)：接管复查时普通 hosted 的9个作业均成功，含两个 shard、fast-checks 汇总与 coverage；Protected Bulk job 仍等待 `bulk-rnaseq-real-execution` 环境审批、无 runner 分配，尚未执行。CI 整体为 waiting，不能标为通过。
- [Lint 35918199896](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/35918199896)：成功，job/step 记录确认 Ruff、snakefmt 和原 lint baseline 均通过。
- [Lock Check 35918203241](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/35918203241)：成功；手动 dispatch 只证明 lock 存在，该入口按原规则不检查事件 diff 中 YAML/lock 的同步更新。

运行时观测：仓库 runner API 返回 `total_count=0`；保护环境存在 required reviewers
和 branch policy，待审批记录存在。因此批准环境也不能代替合格的 self-hosted runner。
本轮没有批准环境、改保护规则、注册 runner 或部署服务；这些前提需另列具体方案。
执行 session 接管上述现有 run，不重复 dispatch；归档失败及真实结果，不能将等待
或跳过记为通过。历史各身份保留当时未验证记录，不默认逐个回跑已被当前版本替代的
提交；PR-6 在线 staging 全链仍是独立验证边界，远端原 Gate 也只有 verify 阶段。
本轮命令、API 状态和 job/step 证据在 `/tmp/helix-remote-ci-jou6h_ma/`。

接管复查（2026-09-24，UTC 2026-09-23 21:04）：三个 run 的 head SHA 与维护分支
仍为已审核 SHA，main 仍为父提交；未再次 dispatch。远端原日志显示两个 shard
3626+3658=7284 项通过，汇总 JUnit 零失败/错误/skip；科学真实工具层43项、平台
真实执行层10项、前端392项与普通浏览器13项通过，coverage 检查通过。
这些是本轮查询/归档的远端执行结果，不是本机重跑，也不代替 Protected Bulk。

runner API 仍返回0；五项环境 vars 均存在，但除 Docker executable 外，四项与
已验收本机坐标不同，未连接未知旧端点。环境 custom branch policy 的规则列表为空，
维护分支的允许状态还需管理员核验，不能因 waiting 推断已通过保护条件。
本机离线资源可复用但服务停止，尚无 Actions runner 的 cache/镜像检出接线。
`local-mirror-only` 是 checkout 的 token 值，不会自动创建本地 mirror；须另行准备
受控镜像、runner 专属 Git 配置和缓存，并实际验证固定 checkout action。
最小方案、资源范围和需授权的 runner注册/坐标/服务/环境审批在
`/tmp/helix-ci-h1-xv32gtkc/protected-runner-plan.md`；本轮没有实施这些动作。
API原始记录、全部已完成job日志和命令在该目录；没有失败作业日志可归档。
Lock Check 手动入口仍只证实 lock 存在，不验证事件差异的YAML/lock同步。

后续授权执行（2026-09-24，UTC 2026-09-23 21:54–22:06）：沿用同一 CI run，
没有再次 dispatch 普通 CI。官方 Actions runner 2.337.0 以仓库专用、单 job ephemeral
方式注册（id 66）；固定 checkout action 已在断网条件下实际从只读 local mirror
完整检出。runner 使用用户另行允许的专用 HOME、离线工具缓存和本轮私有 Redis
Unix socket；只更新四项不同环境坐标，Docker executable 不变；只新增精确维护分支
准入规则，不删除已有保护。再次核对唯一待投递任务、分支 SHA 和前提后，仅批准
既有 pending deployment 一次。实际 hook 接受的 run/job/SHA 与授权目标一致。

[Protected job 107375096912](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/35918197084/job/107375096912)
最终为 **failure**，CI run 整体也是 failure。精确 checkout、锁定工具、安装及原 runtime
verify 均通过；原 Python Gate 为 **1 passed、3 failed、0 skipped、77 deselected**。
通过项是实际取消及清理；超时项报 `accepted lifecycle cleanup is incomplete`，
STAR/Salmon 与 rapid-quant 两项在入口因 Docker socket 消失失败。后续四个 Node/浏览器
步骤被跳过，浏览器产品链未执行；缺少产品 evidence 的检查失败不能记为浏览器失败复现。
`Require no managed containers` 失败源于连接不到 socket，不能据此认定仍有残留容器。

直接原因在本轮临时容量监护脚本，未证实产品缺陷：`du` 遍历专用根时，containerd
snapshot 和 Nextflow `engine/tmp` 下的临时 class 文件正被删除，返回 ENOENT/退出 1。
脚本的有限重试只接受 snapshot 路径，混合错误被拒绝，随后主动停止专用 daemon。
用保存的原 stderr 调用原判断函数，已复现“snapshot 单独允许重试、混合/Nextflow 不允许”。
最后完整容量样本约 75.8 GiB、可用约 274.2 GiB，未达到预算；不采信失败 `du` 的不完整计数。
最小后续方案为仅修正该临时监护的并发删除识别与有限完整重测，保留预算、未知错误和
持续失败的拒绝，先验证再重启；本轮未重启或重跑失败 job，未修改产品、测试、workflow 或身份。

收尾核对：ephemeral runner 自动注销（API 返回0），Listener 退出；专用 Docker/containerd
及 Redis 停止，三个 socket 均不存在，停止前专用 daemon 的最后容器清单为空，宿主可见
范围没有本轮匹配的存活进程。最终 checkout 含 untracked 的 clean 检查通过，1147 个
tracked 文件字节/模式及 90 路径范围不变，Bulk 111 项 aggregate 仍为
`12fa7ec0d8c46d8689ad3042636b677468dca1957765eca2425c1193c174d62a`。
job 原日志、带服务端摘要核验的原 artifact、失败因果与清理记录在
`/tmp/helix-protected-runner-f_t4nmkx/`。runner 默认清理了 RUNNER_TEMP，未外置的科学
私有诊断只保留上传的摘要及 pytest 失败链，不声称完整私有日志仍在。
本机已验收 Gate 与本次远端失败分别记录；G0/旧身份仍为当时未验证，PR-6 在线 staging
完整链没有新增验证。Hi-TrAC H1 规划不变，未推进 H2。本轮收尾待独立验收。

### Protected Bulk 定向重试与容量监护修正（2026-09-24）

前次失败归因及收尾已获独立复验确认，见
`/tmp/helix-protected-failure-review-luj0mump/report.md`；上述首轮失败证据仍保留。
本轮仅修改新的临时容量监护及环境接线，没有修改产品、测试、workflow、超时、
断言或执行身份。新监护逐行校验 exit 1/ENOENT/完整路径，只允许专用 snapshot 树
和本次 runner 指定 workspace 的 Nextflow `engine/tmp/nxf-*` 并发删除进入有限重测。
最多6次，失败的部分计数丢弃，每次重测前检查可用空间；未知路径、权限错误、持续失败
及预算超限仍停止。52项匹配/替身/预算检查及18项组装复核通过；有界真实 du 对照
24次中13次受控 ENOENT、11次成功，停止创建/删除后完整测量成功。这些不替代 Gate。

只调用一次原失败 job `107375096912` 的 rerun API，未重新 dispatch；实际新目标为
[attempt 2 / job 107415516205](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/35918197084/job/107415516205)。
九个普通 hosted 作业在新 attempt API 中有新的 alias ID，但时间、runner 和 steps 与
首轮成功记录完全相同，未重新执行。新 runner id67 为仓库专用 single-job ephemeral，
hook绑定实际run/attempt/job key/审核SHA；只更新私有 Redis 坐标，其余四项坐标及保护
规则不变。前提核验后，仅批准本次新 pending deployment 一次。

**远端实际结果：job 与 CI run 均为 success，24个步骤成功。** 原 Python Gate
`4 passed / 0 failed / 0 skipped`（77 deselected），原 JUnit零skip检查通过；保留1条
RQ `job.exc_info` 弃用警告，未压低警告。原浏览器产品入口 `5 passed / 0 failed / 0 skipped`，
覆盖强制可用 registry→创建→真实运行→QC/artifact→下载及桌面/移动展示；产品证据含
56项artifact、93项QC和1234字节下载的SHA，另有9张原截图。原path-free证据、显式endpoint
零容器及上传步骤均通过。归档artifact id `10782430842`（23个文件）的SHA256为
`7dde650047a9ed1cd81941235ff53e74fdd0719cd970ac6945b8a34d986958d4`，与服务端摘要一致。

实际执行仍为 `842251e48408dcdaa8b3054403be5e2101951f0e`；停服后完整核对tree/parent、
90路径（59修改/25新增/6删除）、1147个tracked字节/模式及含untracked的clean状态均通过。
Bulk111及manifest/qualification未变，aggregate仍为
`12fa7ec0d8c46d8689ad3042636b677468dca1957765eca2425c1193c174d62a`，未重建身份。

运行中7次真实snapshot ENOENT分属5轮容量测量，最长一轮第4次取得完整计数，监护无失败；没有将这些记录
冒称真实Nextflow混合竞态（混合stderr和Nextflow路径另有固定输入/有界du验证）。
累计占用峰值采样79.44GiB、最低可用270.48GiB，RSS合计峰值采样3.79GiB；CPU限于0–7，
未突破既定预算。采样值不代表连续瞬时峰值。停止时最后完整占用样本79.15GiB；停服后
可用270.79GiB，未再次以root测量整个store占用。

收尾：runner自动注销，API返回0，Listener退出0；原Gate及停服前容器清单均为空；
专用Docker/containerd、Redis和已知PID/startticks均停止，三个socket不存在，root停止
记录无失败或清理错误。可见进程范围无匹配的存活任务；另有31个进程cwd权限受限，
不声称已看见全机所有cwd或证明所有既有服务不变。私有诊断在RUNNER_TEMP外0700目录，
没有进入公开artifact；预期超时负向用例的诊断已保留。资产、缓存、旧失败记录未删除。

证据：`/tmp/helix-protected-retry-bdt4jw34/report.md`、`commands.md`、原job日志/JUnit/
artifact、`monitor-vs-prior.diff`、`roadmap.diff`及停服后身份核对。本轮已通过独立复验，用户接受；
没有提交、推送、改分支或推进H2。本机此前通过与本次远端通过分别归档；G0/旧身份仍为
当时未验证，不由本次追认；PR-6在线staging完整链没有新增验证。Lock Check手动入口的
YAML/lock变更同步验证边界、Hi-TrAC H1计划及其他政策裁决不变。

独立复验依据：`/tmp/helix-protected-retry-review-D4qGAw6T/report.md`。
复验核对原远端日志、JUnit、artifact 与身份/停止记录，没有重跑科学 Gate 或浏览器。
验收仅对应上述精确提交的 attempt 2 / job `107415516205`；下载证据核对 HTTP、
generation/revision、大小并记录 SHA，未另以源文件独立摘要作等值校验。
此次验收状态及后续 H1/H2 文档更新不属于该已测试提交。

### 后续 Hi-TrAC 预处理 adapter：五个 PR

使用独立 workflow ID `hitrac-preprocess`。第一版范围为双端 FASTQ 的 linker
处理、Bowtie2 比对/MAPQ 过滤、PET BEDPE 与 QC，以及平台配置、执行、取消、结果查看
和下载；不包含 loop/domain/peak calling、bigWig 或差异分析。adapter 拥有科学参数、
输入与参考契约、输出解释；继续复用 SQLite 生命周期、RQ worker、schema 表单和
产物发布机制，不另建平台或在通用 API/UI 硬编码 Hi-TrAC。
用户进一步确认：科学处理范围严格等于 `tracPre2.py` 本身，不增加后续分析或
新的科学产物。下列校验与结果解析用于安全接入该脚本，不另行实现一套预处理算法。

实施暂按 H1–H5 五个 PR 排序，与已完成的 PR-2～PR-17 分开；一次一个，完成后停下
独立验收。本审核会话只负责规划与复验；本次没有开始实现 adapter。优先采用固定
原脚本加薄启动器，初版不引入 Snakemake 或重写科学步骤；若真实基线证明该方式
无法满足失败检测、资源管理或恢复契约，先提交具体方案，不在接线时顺带重构。

| PR | 用户可验收的结果与主要落点 | 验收与停止条件 |
| :--- | :--- | :--- |
| H1：固定上游与科学契约 | 记录固定源码、依赖和参考身份方案；确定输入配对、参数、PET/QC 单位、零结果及中间产物政策。形成简短 adapter 契约和有预期 PET 集合的 tiny fixture 设计。本轮已形成[H1设计](hitrac-preprocess-contract.md)，零 PET/无 cis 政策 A 及中间产物私有保留政策均已获用户确认，设计待独立验收。 | 用固定源码逐项核对；明确哪些是原算法、哪些是 HelixWeave 的校验。默认计划为每样本一对 gzip FASTQ，可含多个样本；多 lane 需预先合为一对，首版不自动合并、不按 lane 分别去重。任一样本 all/noBg 为空或需要 cis 分母的 QC 集合无 cis，整个 attempt 不发布成功，保留原输出与诊断；政策 B 留待后续；BAM/裁剪FASTQ及诊断在所有attempt中私有保留，不自动删除，不提供首版下载。两项政策已收口，进入 H2。 |
| H2：原脚本真实基线与锁定 runtime | 在 adapter 专用 runtime/启动入口中固定 Python、cLoops2、Bowtie2、samtools、bedtools、gzip 及依赖；实现受控 staging、参考准入、全新 attempt 和完成核验。原算法不变。退出码收集与独立配对来源核验方案已获用户批准；H2私有资格入口已独立验收（159快速＋24资格，另threads=3/8双样本CLI）；旧反例原样保留。 | 真实原脚本跑已知 tiny FASTQ/参考：linker、短/长 cis、trans、重复、低 MAPQ、短 read；核对解压后的 PET 集合与 summary。缺 mate、ID/数量不匹配、坏 gzip、错索引、子工具失败、空/无 cis、旧输出等均有明确结果。替身不能代替这项科学基线；无法安全封装的上游缺陷先报告。 |
| H3：adapter 配置与执行接线 | 新增 `adapters/hitrac_preprocess/` 的 schema、validation、workspace、command、input/reference binding、build identity、availability；复用 registry、默认 runner、worker 与取消机制。 | 配置校验、快照身份、受控执行、失败/取消/超时及清理通过；不暴露任意 shell 参数或私有路径。运行时未准入保持不可执行；H4 结果契约完整前，公共 execution 保持 not_configured/unavailable，仅隔离资格验证入口可执行。共享协议仅在存在真实消费者时最小扩展。 |
| H4：结果、QC 与通用界面 | 严格解析最终 15 个指标列及额外样本索引列（物理 TSV 16 列）的 summary，映射内部样本 token；登记 unique/all BEDPE.gz 与原 summary，落实 BAM/裁剪 FASTQ 保留政策。沿用通用表单、QC 和下载入口。 | 样本集合、计数/比例、gzip/BEDPE、产物来源与结果 generation 一致；缺产物或部分结果不发布成功。实际验证列表、QC、鉴权及下载字节，不能以工作区 TSV 存在代替平台可用。无现成通用表达时先提交契约方案，不硬编码 Hi-TrAC UI。 |
| H5：全链验收与交付 | 用同一固定 tiny 输入完成真实 API → RQ/worker → 预处理 → artifact/QC → 下载，覆盖桌面/移动与重复执行；补部署、升级和上游耦合台账。 | 与 H2 原脚本基线比较规范化 PET 和 QC；验证失败、取消/超时与服务清理，记录资源和完整身份。运行相关 ENCODE/Bulk 回归及适用 Gate，在精确 clean commit 上验收；科学 tiny 通过不外推真实实验的生物学有效性。 |

首版建议直接锁定 Git tag `v0.0.5` 对应 commit
`de6cc732fa00b408551b9f4272933640c08447f1`，并在 H2 锁定完整工具构建产物。
该版本 `scripts/tracPre2.py` SHA256 为
`c3c4ef4e6287fa4ade97a6f5980d20c81ea345b8e8e13c88ae3b7c4bd3a67aec`。
PyPI 同名 `0.0.5` 不能作为该 Git commit 的等价身份；依赖不能仅写包版本号。
一手依据为[固定脚本](https://github.com/YaqiangCao/cLoops2/blob/de6cc732fa00b408551b9f4272933640c08447f1/scripts/tracPre2.py)、
[系统调用](https://github.com/YaqiangCao/cLoops2/blob/de6cc732fa00b408551b9f4272933640c08447f1/cLoops2/utils.py#L103)、
[QC 实现](https://github.com/YaqiangCao/cLoops2/blob/de6cc732fa00b408551b9f4272933640c08447f1/cLoops2/qc.py#L66)
及[许可证](https://github.com/YaqiangCao/cLoops2/blob/de6cc732fa00b408551b9f4272933640c08447f1/LICENSE)。

以下为H1规划时核出的源码约束；后续H2真实实验和当前实现状态见文末记录：

- 输入按 `_R1.fastq.gz` / `_R2.fastq.gz` 查找；原配对占位检查不完整，`zip` 不核
  数量与 read ID。使用安全内部 token/staging 路径，运行前核配对与 FASTQ/gzip
  完整性，保留原始数据。不能直接把任意用户样本名/路径送入上游 shell 字符串。
- linker 固定为 `CTGTCTCTTATACACATCT`，保留长度至少 10 bp；MAPQ 默认 10。
  首版保持固定算法，不把其 PET 去重等同 Picard/Samtools 去重。参考先限定完整六件
  `.bt2` 索引并绑定参考 identity；仅 `.bt2l` 的支持需另行验证上游入口。
- 缺工具/索引等分支可能退出 0，部分子命令状态被忽略，旧输出会触发跳过。
  启动器预创建全新输出目录；完成条件同时核退出状态、最终格式、样本和产物集合、
  计数一致性。summary 会被写两次，文件出现不表示最终结果完成。
- 原 QC 在无 cis/无 unique 时存在将计数设为 1 的分支；空结果与无 cis 对照必须
  真跑并制定不可用指标政策，不能直接透传成可信科学数值。未经裁决不改原算法。
- `n=1,p=1` 会计算出转换 `n_jobs=0`；首版可固定 `n=1,p>=2`，但实际 CPU 预算
  还需计 samtools 固定辅助线程。Python/Biopython/distutils 兼容性及内存上限须
  在 H2 实测；上游旧环境文件不是新的运行锁。

现有接线依据为 `platform/adapters.py`、`platform/registry.py` 和
`services/defaults.py:97–157`：默认 runner 当前通过 Bulk 专有部署查询准入，
单独注册 Hi-TrAC 元数据不等于可执行。应在 H3 接入实际 consumer，不复制完整
Bulk 适配器或为一个新 adapter 建通用插件框架。`adapters/bulk_rnaseq/execution_identity.py:79–183`
已将 defaults、registry、worker runtime 等列入 Bulk 闭包；修改这些文件时按
AGENTS.md 更新对应身份并做适用的 Protected Bulk Gate，不能因目录属于平台而
豁免。Hi-TrAC 自有科学/产品验收与 Bulk Gate 分开记录；当前维护提交的通过结果
不覆盖未来新 adapter 身份。每个 PR 按实际最终 diff 判断范围，纯规划文档不重建身份。

输入先沿用现有已合格的外部输入模式，原始 FASTQ 只读。当前
`services/managed_input_verification.py:68–83` 对受管输入运行仍拒绝准入；本次
接入不顺带解锁 managed input 执行或新建上传/存储系统。工作区的安全配对映射由
Hi-TrAC 启动器负责，不能假设 `WorkspacePlan` 已有 FASTQ symlink 接口。

上游源码、URL/SHA、静态接线调查和本次独立文档 diff 保存在
`/tmp/helix-push-hitrac-plan-ra8cqokr/`；这些静态证据不替代 H2/H5 的真实验收。

### H1 契约与基线设计交付（2026-09-24）

[Hi-TrAC H1契约](hitrac-preprocess-contract.md)已完成固定源码核对及tiny设计，
待独立验收。本轮从官方tag与固定commit重新核验脚本SHA，未使用PyPI替代身份。
明确原PET去重键与QC去重键不同，raw/trim按pairs、mapped按PET行计数；最终表为
15指标加样本索引，原mapping ratio依赖Bowtie2日志结构。给出8对读段的目标坐标、
预期5条unique noBg PET及逐项QC，均为源码手推，必须在H2真实工具下验证。
多样本/预合lane沿用既定边界；推荐BAM/裁剪FASTQ私有保留、不自动删除。

设计交付时，零PET/无cis及中间产物保留均为推荐；后续政策裁决记录见下文。
针对exit0但失败、旧输出跳过和中间summary设计fresh attempt与完成核验；
后置条件能否覆盖所有被忽略的子退出码仍须H2故障注入确认，不能宣称已保证。

证据 `/tmp/helix-ci-h1-xv32gtkc/` 含官方源码/URL/SHA、H1分项、CI接管结果、命令与
相对开工文档独立diff。只改本设计文档及roadmap，未构建环境、运行Hi-TrAC科学任务、
实现adapter、重建身份、提交或推送；未开始H2。本机已验收Gate、远端Protected
Bulk等待、PR-6在线staging以及G0/旧身份当时未验证的边界分别保留。

### H1 政策裁决进度（2026-09-24）

用户已确认零PET/无cis政策A：任一样本all/noBg为空，或需要cis分母的QC集合无cis，
则整个attempt不发布成功，不允许部分样本成功发布。保留原始输出和诊断，明确标识
触发样本、集合及原因；不修改上游计数，不据此评价实验没有科学意义。
政策B（合法PET保留、受影响QC不可用）留作后续能力，本轮不扩大结果与发布契约。

用户同时确认BAM/裁剪FASTQ私有保留、不自动删除，首版平台不提供其下载，公开仅
BEDPE/summary；成功、失败及取消的attempt均保留已有中间产物和诊断。原始FASTQ
始终只读，后续清理或扩大公开范围需另行授权。两项政策已收口，H2进展与停止点见下文；
不进入H3。记录位于
`/tmp/helix-hitrac-h2-bx4wf6lj/`；本次文档不属于已验收维护提交，不重建执行身份。

### H2 原脚本基线与失败检测裁决点（2026-09-24）

该轮交付时H2尚未完成，停在独立审核及最小退出码收集方案裁决点；后续批准和实施见下一节。任务环境和全部探针位于
`/tmp/helix-hitrac-h2-bx4wf6lj/`；隔离源码仅从已验收维护commit导出明确许可路径，
不是完整clean Gate checkout。源码/测试/规则/manifest/qualification均未修改；
源工作区本轮只改H1契约和本roadmap，且保留此前内容。

固定原cLoops2脚本SHA再次核对一致。最终专用科学环境为Python3.11.14、Bowtie2
2.5.4固定构建及124项锁定依赖；真实tiny双样本分别all8/noBg5，15指标及样本索引
与原设计条件预期一致。此前2.5.5日志警告不兼容、bootstrap依赖版本不合格及临时
环境组装错误均保留；通过正常任务局部依赖安装纠正，没有修改上游脚本或校验。

子工具故障注入确认：合法但少一对记录的BAM＋退出73，仍可使原脚本exit0，
并让候选后验错误接受all7/noBg4。固定断言1通过1失败，没有弱化或补造转绿。
另有真实单mate用例：两个不同read ID的孤立记录得到末记录自身配对的cis PET；
必须核每条输出PET是否具有同名两端BAM支持，不能只重跑同一转换工具。
已按用户要求停止仅靠后验的实现方向，交付透明子工具退出码收集和配对后验方案，
不擅改科学算法、过滤政策或公共契约。没有发布资格启动器或正式完成标记。

原脚本边界输出与政策A诊断记录已保留；正式输入/参考准入、完整故障矩阵、取消及
CI接线尚未完成。本轮不是H2完成验收，也未推进H3。Bulk111逐文件及原aggregate
核对不变，不重建身份、不重跑Protected Gate；维护精确提交的已验收结果不能覆盖
未来Hi-TrAC代码。G0/旧身份当时未验证及PR-6在线staging边界保持不变。

### H2 透明退出码与独立配对来源核验（2026-09-24）

用户已批准技术方案，独立审核支持旧退出73和singlemate反例；H1政策A不变，
不是NA映射或部分发布。此次正式实现位于
`/tmp/helix-hitrac-h2-impl-gzfw7n63/checkout/`，只从维护commit
`842251e48408dcdaa8b3054403be5e2101951f0e`选择许可源码，并引入开工H1/roadmap文档。
这是有明确来源的准备副本，不是完整clean Git checkout。源工作区保持不变，
本次文档与代码不属于已经通过远端Protected attempt2/job107415516205的维护身份。

新增私有输入/参考/工具准入、完整原子调用记录、产物后验、按qname磁盘关联的BAM双mate
来源核验、取消/超时进程组清理及最后原子完成标记。原tracPre2和cLoops2包字节不变，
沿用runtime-v3固定124项构建，新tiny由原生成器设计及真实bowtie2-build产生六索引。
原rm范围也纳入计划，仅SAM及两个中间QC，不新增清理；BAM/裁剪FASTQ私有保留。
固定脚本Bowtie2取得stat的正确位置为245（不是247）。

最终快速层159通过；显式真实资格层24通过（12项无替身科学、10项工具故障、2项取消/超时），
均零失败零skip，重复执行不累计独立覆盖。详见本轮report/commands；真实双样本每样本all8/noBg5、
全部15指标和PET多重集合不变，MAPQ10/17均验证。退出73的完整部分BAM、真实singlemate、
零/无cis及多样本一例失败均整批拒绝，保存上游原输出。进程/故障替身与真实无故障科学
分开计数；旧候选accepted=true和旧红灯未修改。实现中发现的信号、收割竞态、落盘/取消
窗口及记录结构问题保存同断言红绿证据，不归为上游科学缺陷。

正式快速测试进入现有CI发现规则；专用真实测试使用既有real_execution标记与显式
坐标，不硬编码本机路径，缺环境失败、不skip。本轮未配置或运行远端Hi-TrAC资格job，
未重跑完整平台/Bulk Gate。字段、输入检查和公共发布仍须H3–H5接线验收，不能以私有
完成标记冒充平台执行可用。线程参数/离散资源观察不是CPU/RSS硬限；不同prefix完整重装
可重复性及主动setsid逃逸进程不在已验证边界。

Bulk111项、manifest/qualification及aggregate保持不变，不重建清单。Hi-TrAC新增代码、
调用计划、入口和锁拥有独立实现摘要，不能沿用维护commit的通过结果。完整ENCODE身份
入口未执行；G0/旧身份当时未验证、PR-6在线staging全链边界继续保留。
H2交付时标为待独立验收。随后独立审核159快速＋24资格全部通过，额外threads=3/8双样本CLI通过，用户接受；证据 `/tmp/helix-hitrac-h2-independent-0d3nw0cc/`。H2不包含H3公共接线，原验证边界保留。


### H3 adapter 配置与执行接线（2026-09-24，待独立验收）

本轮只在 `/tmp/helix-hitrac-h3-0uj0yyu1/checkout/` 继承已验收H2精确字节实施，
未覆盖源工作区、H2或审核证据；选择性准备副本不称完整clean Git。新增Hi-TrAC
schema/validation、服务器参考及runtime绑定、输入/参数/实现身份、workspace/command，
注册默认元数据与作者校验。正式命令通过canonical入口消费原H2 qualifier，未改
固定tracPre2/cLoops2、工具锁、H1政策A或中间产物保留。

公共execution始终not_configured（包括runtime已准入），未声明artifact/QC capability，
未开放公开科学执行/成功发布。快照服务的执行可用注入仅在测试中；隔离资格经过正式
adapter/planner/materializer/command/default runner。H4/H5尚未开始。

实际嵌套进程取消红灯促成worker内部最小控制修正：冻结自有PID/starttime树并等待清理，
取消确认不能早于后代终止。保留原RQ horse wait结果、非SIGKILL行为以及SQLite/job/
事件契约，不改默认ProcessRunner。快速、真实tiny、故障和外层取消证据分别统计于
本轮report/commands；不将无Redis的fork horse控制称作公开产品链。

最终快速/相关回归949项通过（560项核心＋389项Bulk adapter/deployment，零skip）；
原H2资格24项与新增adapter/control 8项分层验证，真实科学与故障/控制实验分别记录。
一项既有完整worker-runtime回归受选择性副本/完整ENCODE身份前提限制未完成，
首轮fixture缺件失败与修复后安全定向结果均保留。Ruff、格式和限定文件patch核验见报告。

共享defaults.py和workers/timeouts.py进入Bulk111，已保存旧闭包及正式同步新身份；
新Bulk aggregate为`a568b21c0c6b0e1c39d386b9cf984bc08d2ea29c6cbb3fddb3f2d2ff0894b1ff`，
111路径不变，仅上述两个受控文件改变；Hi-TrAC 52项实现摘要为
`1dfb2dc4a3628da6212678f40f04794171b6f8846baa3f00a598f1ef3e87d866`。
本H3身份的Protected Bulk Gate待补验。本轮无已审核H3完整clean commit、未授权重启
旧Gate服务或远端dispatch，没有以qualification生成代替Gate。完整ENCODE身份入口
因保护profile边界未执行。维护提交842251e的本机/远端通过不覆盖H3；G0和旧身份
“当时未验证”、PR-6在线staging全链继续独立保留。停止待独立验收，不提交或推进H4。


### H3 独立复验返修（2026-09-24，完成，待独立验收）

独立复验 `/tmp/helix-hitrac-h3-review-u9bvnbg9/` 确认三项缺口：已冻结自有子进程
的身份读取失败被当作消失（C1）；根在完整冻结前死亡并改父后代，空扫描仍被当作
清理完成（C2）；共享runner准入依赖Hi-TrAC helper而原Bulk111闭包未覆盖（I1）。
故障注入不证明默认发生频率；I1只观察到白名单变化，没有执行候选程序，不表述为
已证实的任意命令执行漏洞。返修只在新选择性隔离副本实施，保留H3原交付和审核证据。

worker局部控制边界严格区分不存在/身份变化、仍存在与无法确认；只沿自有进程的
thread children发现下一代，逐层冻结。未知状态、失根或发现不完整均拒绝成功清理
确认，不杀无关收养进程；正常路径保留RQ原wait状态、非SIGKILL及pre-setpgrp保护。
这是拒绝错误确认，不承诺自动找回失根后代；测试按已知身份收尾不能冒充产品恢复。
ProcessRunner、SQLite生命周期、持久job identity和事件语义均未重构。

共享runner实际控制依赖进入Bulk正式闭包，外置工具锁由受控部署代码固定摘要；
没有机械纳入全部Hi-TrAC科学模块。旧新精确恢复材料、正式生成器结果和定向测试
记录于 `/tmp/helix-hitrac-h3-fix-jgvb1rzl/`。最终587项快速/相关回归＋32项真实资格
全部通过（零失败/skip），另有真实SQLite25条组合断言；资格层17项科学正负例、
11项子工具故障注入、4项取消/超时分别报告，不重复累计重跑。C1/C2同字节红灯
1通过/3失败→4通过，I1同字节9失败→9通过。getpgid权限错误遗漏在终审另获红灯
后同边界补齐。Ruff/格式通过；限定diff和原生重放结果见交付。

Bulk闭包111→118，新增7项实际准入控制依赖，正式同步后的aggregate为
`947b5757a9bf27748f1e4505e6e4733c1c2c5bc70cb9d7220d817c81584814b4`；
Hi-TrAC52项实现摘要为`b9a20ee8f849c3c2f7314e01d2eab9c4fb59abe6bf65e6986cf6da3bdc3382e2`。
旧H3及本轮中间/最终身份恢复材料均保留，persistence contract字节未变。
公共execution仍为not_configured，H1政策A及原科学实现不变，未进入H4/H5。

**用户本轮裁决：** 完整Protected Bulk Gate延至H5完成后的集成收尾统一执行；
H3/H4仍执行各自定向验证及必要身份同步，不以完整Gate尚未运行单独阻塞交付。
最终Gate只覆盖最终精确版本，不追认历史或中间身份；维护提交842251e的既有通过、
G0/旧身份当时未验证及PR-6在线staging边界继续分别保留。本轮不启动Gate或旧服务。

原planned-run完整worker-runtime测试仍未完成：旧实际失败首先是选择性副本缺tiny
fixture；之后完整ENCODE指纹需要保护profile是源码核对，不是本轮运行结论。真实
SQLite CREATED种子25条组合断言另行验证，不能代替原PLANNED/Redis-RQ持久链。

### H3 早期并发返修与 H4 连续授权（2026-09-24，阶段 A 本轮验证通过）

用户授权本轮先完成 H3 早期停止／等待窗口返修，定向验证通过并冻结后直接继续
H4；两个阶段独立保存基线、diff、身份和证据，最后一起待独立验收，不进入 H5。
独立复验 `/tmp/helix-hitrac-h3-fix-review-XQf8z76R/` 支持上一轮 C1/C2/I1 顺序与身份
修复，但发现真实 RQ stop 已登记、getpgid/pre-setpgrp 查询仍在锁外的并发确认缺口。

本轮仅在 worker 停止边界补同步：早期准备、重试、所有退出与最终清理决策受同一锁
约束；真实 wait4 仍在锁外。等待确认要求匹配 stopped job、horse PID/starttime 的
完成证明，也拒绝 RQ 已写停止标记但 kill 尚未进入的未确认窗口。未改持久生命周期、
ProcessRunner 或科学代码；未知／失根仍拒绝确认，不承诺自动找回全部失根后代。

证据 `/tmp/helix-hitrac-h3-h4-l8ug4hfs/stage-a/`：同一正式9条并发测试字节，旧实现
3 failed/6 passed，新实现连同原50条 worker 回归59 passed；相关 registry、默认runner、
Hi-TrAC/Bulk身份与准入239 passed；分层真实资格9 passed（5项真实科学正负例、2项
子工具退出73故障注入、2项真实嵌套取消/超时）。均零skip。真实双样本仍all=8/noBg=5，
全部15指标与PET基线一致，旧attempt与政策A拒绝保留；没有运行完整Redis/RQ持久链。

正式生成器同步Bulk118项，路径集合不变，仅worker字节改变，aggregate为
`dbe851ece25be24f959e81dbbcf345df5a7f3db528b89aef02c9c2715e8082a5`；Hi-TrAC52项摘要为
`9fb8516d7dcbe56f3dc3d5e03739e0eb242852e251216e6e7e9efd5ac81b0d99`。旧新精确恢复材料
保留；persistence contract未变。阶段A返修通过本轮验证，待最终独立验收。

完整Protected Bulk Gate仍按用户决定留到H5后的集成收尾；资格生成不是Gate通过。
H3/H4各自定向验证及必要身份同步，最终Gate不追认中间／历史身份。维护提交既有通过、
G0/旧身份当时未验证、PR-6在线staging边界和完整ENCODE身份保护限制分别保留。

### H4 结果、原子发布及样本名契约（2026-09-24，本轮完成，待独立验收）

阶段A已冻结后继续H4，证据 `/tmp/helix-hitrac-h3-h4-l8ug4hfs/stage-b/`。
实现Hi-TrAC Results adapter：原16物理列/15指标、动态MAPQ、精确token映射，重新
核验完整调用/输入/参考/实现/完成标记及原输出和配对来源；仅登记原all、unique/noBg
BEDPE.gz及原summary。BAM、裁剪FASTQ、请求和诊断私有保留，不自动删除。
任一样本不满足H1政策A仍整批拒绝发布，不部分发布或修改原科学结果。

真实消费发现既有两步发布会在后续QC拒绝时已留下artifact。用户批准限定内部bundle：
全部候选先验，再以同一事务提交artifact、QC、generation、attempt、事件与publication；
保留原表结构、公共字段及ENCODE/Bulk普通路径。回滚、并发观察、幂等和陈旧attempt
有两repository及原服务回归。Results子类真实能力完成后按原服务器准入开放execution；
未配置或完整性失败仍拒绝，H3基础类仍关闭，没有用户绕过开关。

原summary的完整比例超过公共Decimal最多12位精度时，仅QC显示候选沿用Bulk的
ROUND_HALF_EVEN，原科学核验/计数/summary字节不改。另一实测缺口由用户独立批准：
只为QC sample_id的原字符集增加ASCII空格，1～255字符；拒绝纯空格及路径/控制字符，
原样保留连续和尾空格，experiment_id/assay及数据注册ID不变。既有metric ID算法未改。
同核心断言8失败/49通过→57全通过；一次新HTTP测试读错响应键的搭建错误已独立纠正、
保留原日志，不混作产品红灯。两次公共方案和决定均在阶段B证据中。

实际验证：结果快速层59通过（其中25项与下列相关集合重叠），sample_id专项57通过，
共享集合249通过及worker17通过，registry/runner/身份/准入239通过，H2快速159通过。
逐命令统计不相加冒称独立覆盖总数。最终当前身份的真实分层11通过、零skip：
3组正式快照→默认runner→SQLite bundle→原认证HTTP列表/QC/下载；MAPQ10/30原结果
对照、篡改整批拒绝、真实singlemate/单样本empty拒绝；另含退出73及真实嵌套取消/超时。
每样本all=8/noBg=5，15指标/PET集合不变，发布链MAPQ17也经过原表头核验。

通用界面152项通过，typecheck/build通过。桌面1440×900、移动390×844使用最终
真实结果数据库：各30 QC/5产物，10次原下载SHA/size一致；连续/尾空格带引号保留展示、
复制、键盘和布局检查通过，服务已停止。正式OpenAPI导出与原JSON逐字相同，生成
客户端不变；前端资产由正式入口同步，identity为
`sha256-acfcf9c82fa169c57c504862507c38287f60198011ebcc0012097c5d1e6bddd3`。
没有手填产物/QC或用合成页面替代该浏览器验证；也没有声称H5完整RQ产品链已通过。

正式工具同步Bulk118→119（新增实际Results准入依赖），aggregate为
`066c7c11cfb39d1cf3f6f86ee9480a6f93e5c73359bba83c555c1014b42ba2dc`；Hi-TrAC63项为
`215f805297c93232994785e37721d1ab083646c972e72c843adf91edacf3b259`。
旧/新闭包及qualification精确恢复材料分别保存；没有DDL或persistence contract变化。
完整Gate未运行，仍留H5后的最终精确版本；维护commit既有通过不覆盖H3/H4新身份，
历史/G0当时未验证和PR-6在线staging全链仍分开记录。完整ENCODE摘要保护边界不变。
H3并发返修和H4均为本轮验证完成、待独立验收；不进入H5，不提交或推送。

### H3/H4 F1/F2 限定返修（2026-09-24，返修完成，待独立验收）

独立复验`/tmp/helix-h3-h4-independent-LUTACHet/`确认两处边界：晚到RQ停止标记
可越过wait返回后的确认，原子bundle失败仍进入成功通知点。本轮证据
`/tmp/helix-h3-h4-fixes-iwa4yxgq/`，只改worker两处及正式回归，源工作区不覆盖。

F1：停止登记与原monitor中每次停止标记消费共享清理同步/证明边界，覆盖原
stopped callback及原handle_job_failure二次消费；wait4仍在锁外。只有匹配
job/PID/starttime的证明才可确认。真实RQ stop、wait4和monitor的屏障反例
5失败/1通过→6通过，连同既有59项worker回归65通过；失败拒绝确认不等于
自动清理失根后代。Redis/回调终点替身范围明确，完整持久取消链未验收。

F2：按现有原子发布内部协议，只允许本次精确完整bundle成功后通知。首次/旧
bundle后QC失败、提交异常、不完整发布和正常/非opt-in对照，经真实内存与SQLite
repository、启用通知服务及内存transport验证；相同断言8失败/10通过→18通过，
原通知及worker35项通过。保留旧结果及科学SUCCEEDED，不发送真实邮件。
最终受影响集合476项通过，包含上述F1的65项，不重复累计；覆盖worker、结果服务、
两repository、registry、默认runner与绑定/身份回归。
当前新身份真实资格9项通过、零skip：默认runner超时与原worker嵌套清理、退出73、
MAPQ10/30双样本all=8/noBg=5、漂移拒绝、singlemate及多样本空集合拒绝，
并经过正式快照/runner/SQLite bundle/认证HTTP QC列表与逐字节下载链。
该链使用进程内ASGI与临时数据库，没有Redis服务、浏览器或H5完整链宣称。

正式入口同步Bulk119项与Hi-TrAC63项，路径集合不变；两处worker字节变化，
Bulk aggregate为`965a9e783eaece90a1c6900d742c8d30d2ea12fbf127c72bc9d1bc8c5810a41d`，
Hi-TrAC为`fe775ae826245c0292e9431742b126c3e31c29e51fe197c76214cab1b16d12f9`。
精确旧新恢复材料保留，未捕获会读取保护profile的完整ENCODE摘要。

旧stage-b前端152通过记录应引用`ui/whitespace-ui-c695ffffcd`；本轮不改前端，
不重复浏览器，也不覆盖旧证据。原planned-run worker-runtime仍有选择性副本
缺tiny fixture及后续保护profile边界；未把测试选择失误或沙箱限制记为产品缺陷。
本轮未修只读QC followup观察，不进入H5。完整Protected Bulk Gate按用户决定
留H5后，最终版本不追认历史/中间身份、G0或PR-6在线staging全链。

### H5 全链验收与交付（2026-09-24，本轮交付，待独立验收）

H3/H4 F1/F2 返修已获独立验收（`/tmp/helix-f1-f2-independent-Z3QODmwn/`）。
本轮在 `/tmp/helix-h5-XTvzkxCP/` 隔离副本（783 项许可文件，逐字节核对，
仍非完整 clean Git checkout）完成 Hi-TrAC 真实平台全链验收。只新增测试与
文档：正式资格入口 `test/hitrac_qualification/test_platform_real.py`、
`h5_stack.py`、`test_publication_real.py`，以及
`docs/development/hitrac-preprocess-operations.md`（部署/准入/升级/上游耦合，
含固定 RQ 停止消费点与 tracPre2 输出/工具身份）。无生产代码、schema、
OpenAPI、前端或身份变化；两闭包路径集合与 aggregate 不变。

真实链证据：认证 HTTP API → validated snapshot → create/preflight/start →
任务专用真实 Redis/RQ → 原 worker 默认装配 → 固定 tracPre2 真实工具 →
SQLite 原子 bundle → 原 API 下载。双样本（连续/尾部空格 ID）各 all=8/noBg=5、
15 指标与规范化 PET 多重集合符合 H2 基线；5 产物、30 QC、5 publication 与
generation 一致；下载字节等于登记摘要；无权/私有拒绝；数据库关闭重开后
终态/结果/下载可读；相同输入新 run 新 attempt、同 snapshot 幂等回放；
成功邮件恰一条且含本次 QC（loopback 捕获 transport，非真实投递）。
取消经真实 RQ stop/monitor/回调持久化并记录 job 与 horse PID/starttime；
超时按既有契约清理嵌套科学进程；政策 A 与原脚本失败整批拒绝、无部分发布。
桌面/移动浏览器消费同一次真实 run，下载 SHA 一致、空格 ID 展示复制保留。

本轮未执行：完整 Protected Bulk Gate（留集成收尾，仅覆盖届时最终精确身份）、
远端 CI、真实 SMTP 投递、真实 worker 内的 bundle 提交失败注入（无安全确定性
接缝，由定向集成层覆盖并明确标注）。H5 资格入口为显式 real_execution 层，
缺坐标失败、不静默 skip，未接入远端 job。G0/历史身份与 PR-6 在线 staging
边界原样保留。本轮不提交、不推送，停止待独立验收。

## Roadmap discipline

New work should advance one product outcome and name its exit evidence. Keep
test counts and coverage floors in the quality baseline, implementation history
in Git, and release evidence in `docs/release-checks/` so this roadmap remains a
short statement of delivered scope and future boundaries.

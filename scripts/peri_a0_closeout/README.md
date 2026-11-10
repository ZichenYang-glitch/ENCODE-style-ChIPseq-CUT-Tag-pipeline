# Peri A0 范围受限的收尾候选

归档阅读入口是 [A0 归档](../../docs/development/peri-a0-archive.md)。下面的 `/tmp` 命令是 2026-11-09 这次实验的历史调用记录，不是阅读结论或重新执行的唯一路径。Python 脚本字节与当时实验一致；本说明增加了这段归档入口。`build_candidate.py` 的默认目录仍指向那些历史路径，重新执行必须显式传入任务目录。

2026-11-09。本目录是独立实验，保留 [原版 r1/r2](../peri_a0/evidence/RESULTS.md) 和 [上一候选 q1/q2](../peri_a0_revision/evidence/RESULTS.md) 的负结果。**最后一项ACP单测已在8GiB单编译任务scope补测通过：固定版本、固定补丁和已声明边界内的A0验收通过，交回主会话审阅。** 见[最新裁决及复现命令](evidence/acp-final/RESULTS.md)；[原收尾报告](evidence/RESULTS.md)保留补测前的阻塞事实。不自动进入 A1，不决定长期 fork。

## 相对 q2 的变化

- 可信启动环境 `PERI_PREDICTION_POLICY=off`，进程首次读取后冻结。未设置保持上游默认；其他字符串（包括 `false`）明确报错。本原型要求回执为关闭，不允许默默用默认。
- `spawn_prediction` 在 clone、辅助任务创建、历史上下文读取和模型调用之前返回。窄只读 ACP `peri/host_policy` 返回有效策略、auto-compact 设置和 gate 计数；不暴露为模型工具或用户 ACP 方法。`prediction_started=0` 指未通过启动 gate，不是“已完成任务数”。
- 严格策略中每个完整身份必须能匹配实际工具；fixture 的参数 schema 也必须匹配。初始装配、refresh 发布前、startup replacement 发布前、最终 pin 都核验。两个固定 meta 工具不锁 schema 值，依赖固定 binary、实际 provider 工具面及源码核验；不声称任意版本通用。
- 宿主每回合前后查策略/连接；每次 `tools/list` 返回前比较完整固定目录，包括 description/schema。资格探针在 blocked 内明确标记为合成流量；真实拒绝、实际 provider 目录和回合结束 prediction 抑制均核验后才 ready。A0 全部是受控合成会话，不接用户会话。
- 自动 compact 用现有配置关闭；用户入口拒绝 slash 命令及任意 ACP，禁止 load/resume/fork/config/mode 调整。没有放宽“无 tools 请求”的 stub 规则。

完整 [上游基线补丁](evidence/candidate-testfix.patch) 与 [q2 增量](evidence/q2-to-closeout.patch) 分开保存。8 个 Rust 文件，继续用官方 loop/权限/ACP 装配。初次 Rust 目录测试的桩缺 wire 身份；[补测身份说明](evidence/testfix-identity.json) 证明修正只涉及测试文件，重新构建 binary SHA 不变。不是生产修复后沿用旧 binary。

## 运行边界

继续使用既有 bwrap 私有 mount/user/pid/net/ipc/uts namespace、清空环境、close_fds、唯一模型 UDS、自有 state、只读 runtime/二进制。bridge/stub/broker 是可信宿主域；Peri 只有 `/work` 私有目录及本地 relay 网络。无平台凭据或控制端点。`sandbox.py` 是实际参数契约，`f2/report.json` 的 OS 探针实际执行 open/connect；不把“没注册文件工具”当作 OS 隔离。

运行 scope：512 MiB 总内存、swap=0、1核 CPU 配额、128任务；`resource_guard` 先检查实际 cgroup 值。relay 每方向64KiB、最多2连接、1秒建立/2秒空闲/8秒总期限、1MiB转发；stub 256KiB body、256请求、2 handler；ACP 单写256KiB/1秒、单帧1MiB、2000事件/8MiB。详细继承边界及此前限额触发证据仍见 [既有运行说明](../peri_a0_revision/README.md)。本次实际复验 native FD 故障，不冒充 native cgroup OOM。

这些是实验限额。没有总磁盘 quota、seccomp 白名单、任意自由文本保密或共享内核漏洞防护证明；不把抽样 RSS 当作总资源限制。宿主同 UID 恶意进程不在本实验威胁模型内。

## 复核与复现

先核对现成材料，不必重跑旧负结果：

```bash
cd /tmp/helixweave-peri-design-20261108
(cd scripts/peri_a0_closeout/evidence && sha256sum -c SHA256SUMS)
sha256sum /tmp/hw-peri-a0-close-20261109/build/target/release/peri
python3 -S scripts/peri_a0_closeout/test_contract.py
```

重新执行时选尚不存在的短 evidence 路径和新 scope 名；保留 bwrap 和真实 cgroup 准入。权限不足时停止并报告，不能改成普通子进程。

```bash
systemd-run --user --scope --unit=hw-peri-close-review \
  --property=MemoryMax=512M --property=MemorySwapMax=0 \
  --property=CPUQuota=100% --property=TasksMax=128 \
  /usr/bin/python3 scripts/peri_a0_closeout/run.py \
  --source /tmp/hw-peri-a0-close-20261109/peri \
  --build-manifest /tmp/hw-peri-a0-close-20261109/evidence/build-testfix.json \
  --evidence /tmp/hw-peri-close-review
```

`--targeted` 仅资格与 prediction 定向阶段，输出的 `not_measured` 不能当作完整矩阵通过。完整命令不加此选项。raw capture 含合成会话、固定测试 key 和任务路径，不含真实账号。

Rust 定向测试使用同一私有 target，依次执行，禁止并发写 target：

```bash
systemd-run --user --scope --unit=hw-peri-close-rust-review \
  --property=MemoryMax=4G --property=MemorySwapMax=0 \
  --property=CPUQuota=400% --property=TasksMax=512 \
  /usr/bin/python3 scripts/peri_a0_closeout/verify_patch.py \
  --source /tmp/hw-peri-a0-close-20261109/peri \
  --manifest /tmp/hw-peri-a0-close-20261109/evidence/build-testfix.json \
  --output /tmp/hw-peri-close-rust-review
```

同一 driver 可用 `--crate peri-acp --pattern trusted_prediction_policy_is_strict_and_default_compatible`；ACP测试crate默认优化编译曾触发4GiB上限，可加 `--test-opt-level-zero`，仅用命令级Cargo包设置降低该测试crate优化，不修改运行binary/全局配置。缓存缺失时 `--online-locked` 仅向任务本地补取Cargo.lock锁定依赖；输出目录必须新建。测试数为0或exit非0均失败。transport测试在运行 scope 内执行 `python3 -S scripts/peri_a0_closeout/test_runtime.py`。

需要重建时在新目录克隆固定 c7d154…，应用 `candidate-testfix.patch`，独立复制已核验缓存，不能共享可写 target：

```bash
HW_NEW=/tmp/hw-peri-close-build-review
mkdir -p "$HW_NEW/evidence" "$HW_NEW/build"
git clone --no-hardlinks /tmp/peri-source-20261108 "$HW_NEW/peri"
git -C "$HW_NEW/peri" rev-parse HEAD  # 必须为 c7d154…
git -C "$HW_NEW/peri" apply \
  /tmp/helixweave-peri-design-20261108/scripts/peri_a0_closeout/evidence/candidate-testfix.patch
cp -a --reflink=auto /tmp/hw-peri-a0-close-20261109/cargo-home "$HW_NEW/cargo-home"
cp -a --reflink=auto /tmp/hw-peri-a0-close-20261109/build/target "$HW_NEW/build/target"
mkdir -p "$HW_NEW/peri/target"
cp -a /tmp/hw-peri-a0-close-20261109/peri/target/peri-rmcp-patches \
  /tmp/hw-peri-a0-close-20261109/peri/target/peri-hyper-util-patches \
  /tmp/hw-peri-a0-close-20261109/peri/target/peri-turso-patches "$HW_NEW/peri/target/"
systemd-run --user --scope --unit=hw-peri-close-build-review \
  --property=MemoryMax=4G --property=MemorySwapMax=0 \
  --property=CPUQuota=400% --property=TasksMax=512 \
  /usr/bin/python3 scripts/peri_a0_closeout/build_candidate.py \
  --task-root "$HW_NEW" --source "$HW_NEW/peri"
```

依赖保留的 Rust/Cargo1.99.0只读工具链和锁定缓存，来源及摘要见构建 manifest。未初始化的 peri-cool 文档子模块只记录 gitlink，不声称参与构建。未知版本/未记录diff/binary不匹配均被运行器拒绝。

主 agent 完成集成、构建和验证；Codex 两个只读 subagent 分别复核源码路径和证据，没有共用可写 target。没有相关科学/产品 skill 适用于这个合成原型，未调用不相关流程。若以后获准提交，实际参与者署名 `Co-Authored-By: Codex <noreply@openai.com>`；本轮不提交。

## ACP 最后一项补测（追加，2026-11-09）

真实peri-acp原单测1/0/0通过；Rust相关合计19项，原native23项和宿主11项不重复运行。
仅test driver增加单任务、profile/环境/资源记录；源码、断言、Cargo.lock及运行binary均未变。
仍使用release/fat LTO/codegen-units=1，仅peri-acp测试包opt-level=0；首轮8GiB通过，
未使用非LTO备选。cgroup内核峰值6.2784GiB属于编译scope，不是native运行需求或最低预算。
本轮scope已退出、无遗留进程/socket/target锁。旧4GiB OOM及错误配置记录不改写；
父目录原SHA256SUMS与acp-final的新SHA256SUMS分别校验。A1未授权，不决定长期fork。

# ACP 补测完成：固定受限 A0 验收通过

2026-11-09。**固定版本、固定补丁和已声明边界内的 A0 验收通过，交回主会话审阅。** 本轮仅补齐此前因编译 OOM 未执行的 prediction 策略单测，没有修改 Rust 源码、测试断言、Cargo.lock、依赖版本或 native 运行配置；没有重跑23项native矩阵。

这不代表A1已获准、真实平台集成通过、生产容量已确定或长期fork已决定。原报告、4GiB OOM、错误配置及未执行记录均保持原字节；[旧收尾报告](../RESULTS.md) 描述的是补测前状态，本文件是追加裁决。

## 实际执行与构建条件

测试：`host::prediction::tests::trusted_prediction_policy_is_strict_and_default_compatible`。

| 项目 | 实际结果 |
| --- | --- |
| crate / harness | 真实 `peri-acp` crate 的 `--lib` Rust test harness |
| 用例结果 | **1 passed，0 failed，0 ignored，0 measured，1219 filtered out**；exit 0 |
| 原断言 | None保持上游默认true；off为false；空串、false、on、OFF、typo均报错。测试文件SHA未变 |
| build profile | `--release`，继承fat LTO、codegen-units=1、默认opt-level=z；仅对peri-acp包用命令级`opt-level=0`覆盖 |
| Cargo / 缓存 | 官方`cargo-rmcp-patched.sh`，`--locked`，离线；任务原有独立CARGO_HOME和target；没有并发构建 |
| 环境差异 | 相对归档运行binary构建，只有`CARGO_BUILD_JOBS: 4 → 1`；profile覆盖在命令参数，不写共享配置 |
| scope | MemoryMax=8589934592（8GiB），MemorySwapMax=0，CPUQuota=400%，TasksMax=512；实际值见report |
| 首轮耗时 | driver开始至报告结束94.076秒；Cargo报告编译1m33s，测试本身0.00s |
| 备选非LTO | **未执行**。首轮已通过，不再增加尝试 |

[完整报告](fat/rust/report.json)、[原始Cargo日志](fat/rust/trusted_prediction_policy_is_strict_and_default_compatible.log)、[scope命令](fat/command.json)。日志中的实际ACP rustc命令含`-C lto=fat -C codegen-units=1 --test`；opt-level=0使用rustc默认值，命令级Cargo覆盖见planned_command。不是把降低包opt-level误称关闭LTO。

本轮仅窄改[test driver](test-driver.patch)：可指定构建任务数、详细日志及备用非LTO命令参数，记录profile、环境差异、开始/结束cgroup资源；[执行时源码](verify_patch-executed.py)留档。备用参数没有运行证据，不算额外测试。宿主API、协议、模型stub、工具、OS隔离及512MiB native限制均未改。

## 资源观察与清理

开始前宿主MemAvailable为173767140KiB（约165.7GiB），未见cargo/rustc等并发编译进程，本任务target无锁持有者；仅检查进程名称及是否引用任务target，没有停止其他任务。见[前置检查](preflight.json)。

内核cgroup `memory.peak`为 **6,741,385,216字节（6.2784GiB）**，scope退出前由driver读取；memory.events的max/oom/oom_kill/oom_group_kill均为0，pids.peak=27。该峰值是整个编译scope（driver、Cargo、compiler/linker及计入cgroup的文件缓存）的总内存，不是rustc单进程RSS，也不是native运行需求。外部观察器200ms采样同样捕获该内核峰值，见[观察记录](fat/resources.json)。8GiB是本次授权上限，不是已测最低需求；单次缓存构建不能确定生产容量。

[清理](cleanup.json)：scope Result=success、inactive/dead、ControlGroup为空；6个观测PID均不存在，本轮任务无socket、target无锁持有者。没有启动模型stub、Redis或native服务；旧OOM失败元数据未清理，其他任务未操作。

## 身份与验收适用性

| 身份 | 本轮前后 |
| --- | --- |
| HelixWeave | `622778f64e59f935d822d92a91bade78807dabce`，未迁移worktree或改标main |
| Peri commit | `c7d154c9985cbe08ba25e964cb7468d15014e6d0` |
| 完整patch SHA256 | `06272e1a66423a9e562509a38d044cd8506a1c10e61ef1fd3fe239107680587b`，8文件补丁未变 |
| Cargo.lock SHA256 | `0692683359fb75073eca3d14b364071bc02058be7b11e0335ccaf2b0f314f5dd` |
| 运行binary SHA256 | `60144bbe0919992b8913b9031696b2e30b1e9f82d95048ae3b07c887d5d85517`，未覆盖或重建 |
| prediction_test.rs SHA256 | `5fa87543e25f94b195f28e658f1232364225d25fa7e361b997a217125c777c64`，断言原样执行 |
| 本次test harness SHA256 | `c483dd2556a464b77c4fe0f2f0d54903d512c071179b3e3f77c872a858b702e4`，独立于运行binary |
| 工具链 | Rust/Cargo1.99.0，与build-testfix记录一致；wrapper SHA亦一致 |

[前置身份](identity-before.json)和[完成后身份](identity-after.json)包含完整核对。运行binary构建使用原release opt-level=z/fat LTO/codegen1，本次test harness对peri-acp用opt0；这只改变测试编译方式，不冒充对运行binary的新执行。

A0当前矩阵：

| 证据组 | 当前状态 | 本轮处理 |
| --- | --- | --- |
| 固定受限native A—F场景 | 23/23通过，25实例、128实际provider请求 | [f2](../f2/report.json)不重跑。所有运行代码、配置、sandbox、协议和binary身份不变 |
| Rust策略/filter/目录 | 18个不同用例通过 | 保留既有结果，不重复计数或重复运行 |
| ACP prediction策略 | **新增1项通过** | 本次真实crate原断言；Rust合计19个不同用例通过 |
| 宿主contract/transport | 7+4通过 | 代码不变，沿用原证据 |
| 静态/身份/清理 | 通过 | [本轮检查](checks.json)，没有剩余本轮必需门 |

[既有完整矩阵](../RESULTS.md)中的分层边界和未测范围保持：同进程多会话、被禁用的native恢复、native cgroup OOM恢复、wire目录主动推送全链、任意自由文本保密、大规模长会话和真实平台均没有被本次补测“升级”为通过。目录generation/宿主回调、native FD故障、坏帧与审批超时仍按原实际验证层次解释。基础幂等、确认失效、撤销及丢响应仍是未来B2前置条件。

## 重复核对或补测

优先复核现成证据；本轮不需要再运行native矩阵。新旧摘要分开：父目录SHA256SUMS仍校验原195文件，本目录SHA256SUMS校验追加材料；历史原报告不改写。

```bash
cd /tmp/helixweave-peri-design-20261108
(cd scripts/peri_a0_closeout/evidence && sha256sum -c SHA256SUMS)
(cd scripts/peri_a0_closeout/evidence/acp-final && sha256sum -c SHA256SUMS)
```

如需独立重现，先检查可用内存和target锁，使用尚不存在的输出目录和新scope名，保持一次一个构建：

```bash
systemd-run --user --scope --unit=hw-peri-acp-review \
  --property=MemoryMax=8G --property=MemorySwapMax=0 \
  --property=CPUQuota=400% --property=TasksMax=512 \
  /usr/bin/python3 scripts/peri_a0_closeout/verify_patch.py \
  --source /tmp/hw-peri-a0-close-20261109/peri \
  --manifest /tmp/hw-peri-a0-close-20261109/evidence/build-testfix.json \
  --output /tmp/hw-peri-acp-review \
  --crate peri-acp \
  --pattern trusted_prediction_policy_is_strict_and_default_compatible \
  --build-jobs 1 --test-opt-level-zero --verbose
```

本轮由Codex主agent完成，没有新增架构研究或subagent构建。原参与者记录保留；未提交、push、开PR、合并、部署，未调用科学/产品Gate或真实账号。A1、长期维护路线及一次确认覆盖create/preflight/start仍交回主会话决定。

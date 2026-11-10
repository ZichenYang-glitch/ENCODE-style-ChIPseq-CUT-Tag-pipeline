# Peri A0 接法修订：工具过滤有进展，整体验收仍不通过

2026-11-09。当前候选仅用于本地合成实验，**不能进入 A1**。完整对照、A—F 矩阵及停止原因见 [结果](evidence/RESULTS.md)。原 [A0 r1/r2](../peri_a0/evidence/RESULTS.md) 未改写；q1/q2 是新证据编号，不是原版通过结果。

本目录继续使用 HelixWeave `622778f64e59f935d822d92a91bade78807dabce` worktree；已验收平台 main `b2f15868c6f74ab8764a73ae6fe4394f67304475` 仅更新取消依赖。Peri 候选基于 `c7d154c9985cbe08ba25e964cb7468d15014e6d0` 加三文件实验补丁；当前官方比较对象固定为 `3288b07bf871276026ae0677129821b096e75118`，没有采用其发布包或升级数据库。

## 最小选择

| 方案 | 实际入口、修改范围 | 隔离及维护成本 | 本轮决定 |
| --- | --- | --- | --- |
| A：现成配置 | ACP 顶层仍 allow-all；CLI print 明确拒绝未实现的 allowed/disallowed 参数；当前版本关键 catalog/filter/tool-search 无增量 | 配置不能代替 dispatch 过滤；OS 隔离另需保持 | 不满足。见 [固定增量](evidence/upstream-comparison.json) |
| B：官方组件宿主 | Rust `build_stage_context`、`MiddlewareChainAssembler`、`try_with_filter` 可复用 loop/权限，但 StageBuildInput 约40字段，ACP 构建桥为 pub(crate) | 属于变化较快的内部装配；需维护协议、生命周期及端口接线；同样需要OS边界 | 相比窄补丁改动较大，本轮不实现 |
| C：完整 catalog 策略补丁 | 三个 Rust 文件：接入可信身份策略、严格 pin 再过滤、严格模式 refresh 错误退出 | 保留 upstream loop/权限/ACP；仍与内部 catalog 版本耦合；每次升级重验 | 选择为实验候选，不决定长期 fork |

策略以 `PERI_PROCESS_TOOL_POLICY` 显式注入，首次读取后进程内冻结。未设置保持旧行为；未知 JSON 字段/版本、通配或重复名称拒绝。精确比较 canonical name、MCP server/wire、namespace、builtin、完整 aliases；不是按名称猜测工具是否安全。每实例 broker 仍绑定唯一 serverId/session/generation。仅准许原有两个 meta 工具和两个 fixture 工具。

同一过滤器用于 initial、refresh、startup 工具表；严格策略的最终 pin 再过滤。直接解析、Search/Execute 转发消费该目录。**这不关闭 MCP 自身的后台资源 RPC，也不限制 prediction 等辅助模型调用。** 合法但不存在的策略项目前只是不匹配，并未实现完整“未知策略工具立即报错”；动态目录变化、alias/namespace 等完整行为覆盖仍有缺口，不能称安全接入已完成。

## 运行边界

沿用 bubblewrap 私有 user/mount/pid/net namespace，清空环境、仅 stdio 继承、只读 runtime/CA/二进制和唯一模型 UDS、自有 state。宿主 driver/stub/broker 是可信域，无真实平台凭据。新增策略 JSON 和固定 `TOKIO_WORKER_THREADS=2` 是受控环境项。

- relay：连接期限 1 s，空闲 2 s，总期限 8 s；非阻塞缓冲每方向 64 KiB，总转发 1 MiB/连接，同时最多 2 个连接，不使用无限 `sendall`。
- stub：单请求 body ≤256 KiB，最多 256 个捕获请求；同时最多 2 个 handler；socket I/O 2 s、连接总期限 8 s；失败记录最多256条，脚本动作最多128条。
- ACP：单写 ≤256 KiB、1 s 写期限；读帧 ≤1 MiB；事件 ≤2000条且总量 ≤8 MiB；RPC超时关闭实例；子进程单文件写上限16 MiB。
- 实验 driver、stub、全部 native/relay 后代统一放在任务临时 **user scope**：memory.max=512 MiB、swap=0、CPU=1核配额、pids.max=128。`resource_guard.py` 检查实际 cgroup 值后才启动。不是仅声明 systemd 参数；另有独立64 MiB/16任务/20% CPU的控制器触发证据。

这些是本次实验限制，不是生产预算。没有总磁盘 quota、seccomp syscall 白名单或跨宿主隔离保证。日志单文件限额不等于总磁盘限额；共享内核及宿主同 UID 恶意进程不在证明范围。namespace 内无宿主网络，仅 UDS→stub；未做 syscall/packet 全面外联统计。

## 核对与重复运行

先核对已有证据，无需重新运行旧 r1/r2：

```bash
cd /tmp/helixweave-peri-design-20261108
(cd scripts/peri_a0_revision/evidence && sha256sum -c SHA256SUMS)
sha256sum /tmp/helixweave-peri-a0-revision-20261109/build/target/release/peri
```

以下会重新执行合成实验，应使用尚不存在的短证据路径和新 scope 名。预期 **退出码1**：资格过程中捕获未批准的辅助模型请求并失败，不能解释成成功。默认命令有意只运行资格及启动故障；修正后的未来候选须去掉 `--qualification-only` 才可验收后续场景。

```bash
systemd-run --user --scope --unit=hw-peri-review-next \
  --property=MemoryMax=512M --property=MemorySwapMax=0 \
  --property=CPUQuota=100% --property=TasksMax=128 \
  /usr/bin/python3 scripts/peri_a0_revision/run.py \
  --source /tmp/helixweave-peri-a0-revision-20261109/peri-candidate \
  --build-manifest /tmp/helixweave-peri-a0-revision-20261109/evidence/build-final.json \
  --evidence /tmp/hw-peri-review-next --qualification-only

/usr/bin/python3 -S scripts/peri_a0_revision/test_contract.py
# 同样在上述限额scope中执行本地socket测试：
/usr/bin/python3 -S scripts/peri_a0_revision/test_runtime.py
```

所需 namespace / user bus 能力必须显式可用；失败不能去掉 bubblewrap/cgroup 后声称同一实验通过。未使用 systemd 常驻服务，也未修改用户 manager/global 配置。资源探针 OOM 会结束整个 scope，检查 `Result=oom-kill`，不要将没有探针 JSON 误当成功或自动重跑。

需要重建时创建新任务根，复制原任务缓存而非在原目录写入；使用固定源码和下列补丁。示例路径必须尚不存在：

```bash
HW_NEW=/tmp/hw-peri-build-next
mkdir "$HW_NEW"
mkdir "$HW_NEW/evidence" "$HW_NEW/build"
git clone --no-hardlinks /tmp/peri-source-20261108 "$HW_NEW/peri-candidate"
git -C "$HW_NEW/peri-candidate" rev-parse HEAD  # 必须为 c7d154…
git -C "$HW_NEW/peri-candidate" apply \
  /tmp/helixweave-peri-design-20261108/scripts/peri_a0_revision/evidence/candidate-final.patch
cp -a --reflink=auto /tmp/helixweave-peri-a0-20261109/cargo-home "$HW_NEW/cargo-home"
cp -a --reflink=auto /tmp/helixweave-peri-a0-20261109/build/target "$HW_NEW/build/target"
mkdir -p "$HW_NEW/peri-candidate/target"
# 复制已校验的上游依赖patch缓存；wrapper仍按固定patch摘要校验。
cp -a /tmp/peri-source-20261108/target/peri-rmcp-patches \
  /tmp/peri-source-20261108/target/peri-hyper-util-patches \
  /tmp/peri-source-20261108/target/peri-turso-patches "$HW_NEW/peri-candidate/target/"
systemd-run --user --scope --unit=hw-peri-build-next \
  --property=MemoryMax=4G --property=MemorySwapMax=0 \
  --property=CPUQuota=400% --property=TasksMax=512 \
  /usr/bin/python3 scripts/peri_a0_revision/build_candidate.py \
  --task-root "$HW_NEW" --source "$HW_NEW/peri-candidate"
```

此复现依赖保留的1.99.0工具链和缓存，工具链来源摘要见构建身份。manifest记录完整diff、锁文件、命令/环境、binary SHA；换源码/补丁/binary即拒绝。`verify_patch.py --source … --manifest … --output <新目录>` 运行三个相关Rust测试过滤器，不能与同一target构建并行。

原版与候选构建、原始错误和结果分别保留；不接真实账号、平台 API/SQLite/Redis，不发送上游消息，不提交/push。B2 的幂等、确认变更/过期、撤销、丢响应仍是未来写操作前置条件；一次确认覆盖三段仍待产品决定。

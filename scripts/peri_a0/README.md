# Peri A0：隔离与工具准入原型

**本轮结论：A0 不通过，不进入 A1。** 固定 native 的实际 provider 请求暴露额外工具；完整身份、逐项结果及原始失败说明见 [RESULTS](evidence/RESULTS.md)。

2026-11-09 接续交付：现有调查已完成，未重跑实验或构建。原型及实验仍在
HelixWeave `622778f64e59f935d822d92a91bade78807dabce`；最新 main
`b2f15868c6f74ab8764a73ae6fe4394f67304475` 的 #210 取消修复只更新设计依赖，
不改变 Peri `c7d154c9985cbe08ba25e964cb7468d15014e6d0` 或旧实验身份。
接续静态核验见 [post210-review.json](evidence/post210-review.json)。

本目录仅用于安全调查，不是生产 bridge。设计见 [ADR](../../docs/architecture/peri-agent-integration-decision.md) 和 [契约附录](../../docs/architecture/peri-agent-integration-contracts.md)。不接真实模型、样本、平台 API/SQLite/Redis；两个工具只有合成查询和内存计数。没有 A1、业务确认表、任务协调器或聊天 UI。

## 运行边界与复现

需要 Linux user/mount/pid/net namespaces、bubblewrap 0.9.0、`prlimit`、本机 Python 3.12 与相应只读运行库。当前启动器是本机 Linux 的固定实验装配，不是跨发行版安装器。默认工具执行沙箱禁止 namespace 所需调用时，应给这条实验命令单独执行授权；不得移除 bubblewrap 或修改系统配置后冒充同一测试。

从包含本目录的设计 worktree 运行。源码须是官方 `KonghaYao/peri` 的固定 commit，不使用已安装 Peri，不读取 `~/.peri`。源码取得和子模块初始化属于独立任务目录中的工作；确认 `git rev-parse HEAD` 后才构建。

已有交付先检查身份和证据即可，不需要重复运行：

```bash
git -C /tmp/peri-source-20261108 rev-parse HEAD
sha256sum /tmp/helixweave-peri-a0-20261109/build/target/release/peri
(cd scripts/peri_a0/evidence && sha256sum -c SHA256SUMS)
```

以下仅为以后获准重验时的命令。保留 r1/r2 和原 build.json；需要重建时使用新的
task-root（示例 `...-rebuild` 必须尚不存在），不要覆盖现有构建证据：

```bash
# 已有源码时无需重复 clone，也不要自动更新 main。
git -C /tmp/peri-source-20261108 submodule update --init --depth 1 peri-cool

/usr/bin/python3 scripts/peri_a0/build.py \
  --source /tmp/peri-source-20261108 \
  --task-root /tmp/helixweave-peri-a0-20261109-rebuild

# 单测和 native 复现也不属于本次接续交付的执行步骤。
# 下面复用已核实的原始 binary；重建版必须改用其独立 build manifest。
/usr/bin/python3 -S scripts/peri_a0/test_contract.py
/usr/bin/python3 scripts/peri_a0/run.py \
  --source /tmp/peri-source-20261108 \
  --build-manifest /tmp/helixweave-peri-a0-20261109/evidence/build.json \
  --evidence /tmp/helixweave-peri-a0-20261109/evidence/r-next
```

`--evidence` 必须是新目录；使用短路径以满足 Unix socket 路径长度。失败返回非零，不能解释为“没有执行科学任务所以通过”。本环境 Cargo/libgit2 HTTPS 下载失败后使用了 `build.py --git-ssh`：仅该构建子进程的 Git CLI URL 映射和 BatchMode 生效，不写全局 Git/SSH 配置，不改变依赖 pin。保留首次失败和重试日志。构建会从官方 Rust manifest 校验 1.99.0 工具链归档，安装在 task-root；执行固定源码自带 `cargo-rmcp-patched.sh build --locked -p peri-tui --release --bin peri`，包含上游仓库已有依赖补丁，不能把它们称为本地新补丁。

`build.json` 记录源码、锁文件、工具链、子模块、完整命令/环境、二进制 SHA256 和本地 diff；`run.py` 拒绝源码变化、补丁或二进制摘要不符。每次运行记录所有原型 Python 文件摘要，身份不依赖未授权的 Git commit。重跑会使用已有任务目录缓存；需复核源码及构建记录。产物、缓存与原始日志保留在 task-root，不纳入仓库。

## 最小装配

- `sandbox.py`：每实例独立 bubblewrap namespaces、只读 binary/runtime、私有 state、空 `/work`，环境重建且仅 stdio 继承。`relay.py` 在该 netns 内固定转发 loopback HTTP 到唯一挂载的模型 UDS；不支持动态目标、DNS或转发宿主 FD。
- `protocol.py`：ACP newline JSON-RPC，配置经 length-prefixed stdin 注入。成功路径为 `blocked → ready → closed`，准入失败为 `blocked → closed`（本轮实际路径）；closed 不可重新开放。只暴露带服务端 instance/session/generation 绑定的 prompt，不暴露任意 ACP、模式修改、load/resume/fork。
- `stub.py`：宿主可信 HTTP 假模型，仅 UDS、固定路径、合成 credential；从 spawn 前开始计数，捕获真实请求。首个实际 provider 目录出现额外工具立即关闭实例，不返回工具调用。正常资格请求只用固定合成文本，基础 mode/MCP 目录检查前无模型请求；用户 prompt 始终等全部资格完成。
- 假 MCP broker 仅声明 `lookup` 与 `increment`。每次工具调用检查实例/代次、当前 session、单活动轮次、参数精确内容及一次性审批 grant。Peri 到达 broker 的请求数与 broker 副作用分开统计，避免 broker 拒绝掩盖 Peri 审批失效。它不是 HelixWeave 的项目权限或业务批准实现。
- `fixtures.py` 仅演示封闭字段投影：固定数值/枚举，不传 description/default/examples/message/hint/ID/path。它不能处理任意科学 schema，更不是通用文本匿名化器。别名化、正则过滤均不能保证自由文本不泄漏。

宿主 driver、stub、假工具属于同一可信测试域，仍具有宿主用户权限；只有 Peri/relay 受上述 OS 边界约束。运行库只读目录可读，不声称完全无文件访问。没有 seccomp syscall 白名单、cgroup 总内存/CPU/PID 上限；nofile/core/单文件大小由 prlimit 限制，不能当作总资源预算。共享内核漏洞及宿主同 UID 恶意进程不在本次证明范围。

## 证据的解释

`evidence/` 保存可审阅的结果、固定身份、原始失败摘要和脱敏协议/provider 记录；详见其中结果说明。运行目录保留完整 stderr、消息、OS 探针与采样，仓库副本不含二进制、会话 DB 或下载缓存。

`os_probe.py` 是真正执行 open/connect 的攻击性进程，使用与 Peri 相同的 namespace/mount/env 策略，只多挂探针脚本；其结果证明指定策略的访问边界，不意味着 Peri 曾调用 shell/file 工具。探针只访问任务生成的假秘密、另一假会话、假平台 socket 和本任务监听端点，不探测真实服务。

`test_contract.py` 只验证宿主 broker、绑定和字段投影。真实 native 验收在 `cases.py`；完整资格失败时后续场景记 `not_measured`，不移除场景也不扩大授权。未知 mode、断开的 mode 设置、缺失/额外/替换工具分别是独立故障注入。`mode_timeout` 是 driver 丢掉该 outbound frame，不声称 Peri 自发超时。坏消息测试直接注入同一解析函数，不声称真实 Peri 曾发坏 JSON。

运行时间从 spawn 前计到 ready 或首次 closed，报告两者区别。RSS 每 100 ms 采样本实验 bwrap/relay/Peri 后代树，峰值是采样值、可能重复共享页且不含宿主 stub，不能推出容量或内存上限。关闭时停止本任务实例/stub、检查观测过的后代，不清理其他任务资源；本地文件与缓存保留供审阅。

B2 的重复提交、丢响应、批准过期/变更、权限撤销及执行中断基础验收仍是未来开放写操作的前置条件，本轮未实施。一次 submit 确认是否覆盖 create/preflight/start 仍待产品决定。

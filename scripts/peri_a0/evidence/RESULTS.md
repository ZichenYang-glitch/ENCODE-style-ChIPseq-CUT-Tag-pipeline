# A0 调查结果：不通过，调整接法后重验

日期：2026-11-09（Asia/Shanghai）。**当前固定 Peri native ACP + 反向 MCP 的纯配置接法不满足工具闭合门槛，不进入 A1。** OS 策略和若干失败前置检查通过，不等于接入验收通过。没有改 Peri 源码，也没有补丁版结果。

## 固定身份与决定性证据

- HelixWeave worktree、初次研究及 r1/r2 实验基线：`622778f64e59f935d822d92a91bade78807dabce`，也是实验当时核实的远端 main。后续 #210 已进入 main `b2f15868c6f74ab8764a73ae6fe4394f67304475`；接续交付只读核实该增量，不迁移 worktree、不倒改实验身份、未 commit。
- Peri：`c7d154c9985cbe08ba25e964cb7468d15014e6d0`；远端 main 已漂移至 `63d012d0750ccbb8d9ff4711ef1209a6809e10c8`，没有切换。`agent-v4.6.0` 与本机其他安装版本不在结论范围。
- Rust/Cargo **1.99.0**，x86_64-unknown-linux-gnu；官方 wrapper、release profile、`--locked`、4 build jobs。构建耗时 **353.90 s**，不含工具链下载/首次失败；完整身份见 [build.json](build.json)。源码 diff 为空，Cargo.lock 前后摘要相同。上游已有依赖补丁由固定脚本应用，未新增本地补丁。
- 实际 binary SHA256：`ecb4eb1394c67db2c57bc74373f194d52c6342e3be425ed0980b665ec29d626f`。完整命令、环境、Rust manifest 归档摘要和 submodule gitlink 在 build.json；原型逐文件摘要在 [r2/report.json](r2/report.json)。运行前再次比对源码和 binary。资格实例 native 日志因原始末尾空行以 JSON 的 raw_text 无损保存，raw_sha256 可验证还原字节；其余原始日志直接保存。
- [真实 native 资格请求](r2/A_B_permission_tools/evidence.json)：`session/new` 回报初始 `bypass`；设置 Default 后返回 `currentValue=default`。`mcp/list` 只报两个 fixture 工具，但实际 provider `tools` 为 **ExecuteExtraTool、ReadContextArchive、SearchExtraTools**。其中 SearchExtraTools/ExecuteExtraTool 是原型预期的发现/执行元工具，**多出的 direct 工具是 ReadContextArchive**。同一请求 system 内容还列出非白名单 deferred **DiscoverMCP、mcp_read_resource**。这是实际捕获，不仅是源码推断；不把两个预期元工具误记为额外能力。
- 宿主 stub 在第一个请求处拒绝额外工具并关闭实例，未返回脚本化模型输出。`counter=0`、`lookups=0`、`tool_call_arrivals=0`、`approvals=[]`、`ready_seconds=null`；不能把这些零计数称为 Peri 实际审批拒绝通过。native 还尝试 `resources/templates/list`，broker 按未授权方法拒绝（`denied_callbacks=["inner_method"]`）。不声称已穷尽所有 deferred 路径。

## 原始失败及分类

1. 初次 build driver 错把 peri-cool 当 Rust crate，已修为校验 gitlink/实际 SHA；[原始错误](build-driver-attempt1.log) 属于原型接线错误。
2. Cargo/libgit2 HTTPS 拉取锁定 mio 失败；[首次构建日志](native-build-attempt1.log)、[失败身份](build-attempt1.json) 保留。随后只在构建进程中采用 Git CLI/SSH URL 映射取得相同 pin，构建成功；无全局配置变化，无版本替换。
3. [r1/report.json](r1/report.json) 的 session/new 在隔离环境缺 CA bundle 时 panic；见 [native 原始日志](r1/A_mode_error/native.log)。**r1 的 mode_error/mode_timeout/missing/extra/replacement “passed” 无效**：当时断言接受了任意启动失败，未命中预期故障边界。原始文件保留不改写；r1 不能作为这些验收的通过证据。补挂唯一公开 `/etc/ssl/certs/ca-certificates.crt`、设置固定 SSL_CERT_FILE，并加强精确错误码/注入标记后，r2 重跑。没有增加外部网络、个人配置或秘密访问权。
4. r2 最终阻断是 **Peri 固定版本/配置能力限制**：MCP middleware 无条件注册上述额外 intrinsic 工具，ACP 没有核实到可关闭单项的配置 seam。独立进程和 OS 沙箱不能改变工具目录。当前接法需调整；尚无证据说明所有 native ACP 或整个推荐分层都不适用。

## A—F 验收逐项归属

| 条目 | 本轮结论 | 证据/局限 |
| --- | --- | --- |
| A：spawn 开始计数、blocked 拒绝用户 prompt | 通过该片段 | r2 各 case 从 spawn 前创建 stub；资格请求只有固定合成文本。六个启动故障模型请求/工具到达均 0；完整准入未成功 |
| A：Default 设置失败/超时、未知模式 | 通过 | [mode error](r2/A_mode_error/evidence.json)、[mode timeout](r2/A_mode_timeout/evidence.json)、[unknown](r2/A_unknown_mode/evidence.json)；mode error 实际改错 session ID 并获得 native -32602；timeout 是 driver 丢 outbound frame，1 s 截止；未知模式在宿主拒绝 |
| A/B：缺失、额外、替换工具 | 通过服务器目录拒绝片段 | [missing](r2/A_missing_tool/evidence.json)、[extra](r2/A_extra_tool/evidence.json)、[replacement](r2/A_replacement_tool/evidence.json)，精确 TOOL_CATALOG_MISMATCH。不能以这张服务器表替代 provider 工具面 |
| A/B：实际 provider 与 deferred 工具闭合 | **失败，A0 阻断** | 多出的 direct ReadContextArchive 触发首请求拒绝；同次请求还列出 deferred DiscoverMCP、mcp_read_resource。未主动调用它们，未穷尽后续发现路径 |
| A：实际审批拒绝、审批连接中断/超时无 Bypass | 未测 | 工具闭合先失败，未抵达新的审批边界；Default 名称不是行为证明 |
| B：允许工具成功；伪造 shell/read/SQL/MCP/配置/扩权；参数越界/未知字段 | native 未测；宿主部分通过 | 场景保留在 cases.py。六项 [宿主单元测试](host-unit.log) 验证 grant/范围/参数/任意 ACP 拒绝，不冒充 native 执行 |
| C：OS 策略探针 | 12 个断言通过 | r2/report.json.os_probe：两次真实文件读取、三次真实 connect，分别针对假 bridge 文件、另一假会话、宿主 TCP、pathname/abstract 控制 socket；其余检查 proc/FD/环境、cap/no-new-privs、路径和私有目录写入。`no_host_mounts` 只检查 `/home`、`/run` 不存在，不是完整挂载审计或外网连接穷举。指定访问尝试被 namespace/mount 边界阻止，没有依赖未注册工具来通过 |
| C：两个真实实例交错运行的配置/权限/回调隔离 | **未测，硬门未满足** | 完整资格无 ready 会话，未运行该场景。独立 namespace 理论及单个假其他目录探针不能代替 A/B native 交错证据。同 native 进程多会话研究为未来可选项 |
| D：kill/restart 后旧连接/代次失效、新实例重新准入 | native 未测；宿主绑定片段通过 | host-unit 检查旧代次/closed/跨实例拒绝；未运行完整重启序列。native load/resume/fork 不暴露；自动 continuation 仅源码证据，未运行复现 |
| E：字段白名单/假秘密与注入文本投影 | 宿主投影片段通过，真实 provider 完整验证未测 | host-unit 使用合成 description/default/examples/hint/ID/path/instructions，固定数值/枚举投影；被前置失败阻断的 E_minimum_egress 未发 prompt。不能宣称任意自由文本已匿名化 |
| E：禁止继承个人配置/遥测环境、限制网络出口 | OS 片段通过 | 真实 env keys 记录在 os_probe；私有 HOME/XDG/work/state，无宿主个人配置挂载。Peri 尝试读取的是**本任务私有且不存在**的 `.claude/settings.json`，没有读取个人配置。1 次 native provider 请求证实允许的 UDS→stub 链路；没有 packet/syscall 级全面外联尝试统计，不把未观察到调用当作遥测代码不存在 |
| F：坏协议/模型超时/工具回调失联/进程主动退出 | 对应 native 场景未测 | 六项启动失败和资格关闭路径已运行；其他场景被准入阻断，不能用它们替代审批中断或运行中断 |
| F：失败不伪造成功、进程清理 | 已运行路径通过；整体非零退出 | r2 `passed=false`；所有 root process 已退出，所有观测过的后代 `still_present=[]`、cleanup_errors=[]；stub listener 关闭并删除，runner 已退出 |

没有运行科学流程、完整 CI、真实平台端到端、Bulk Gate，亦没有用真实模型帐号。没有扩大白名单、去掉测试、改 Peri 源码或以普通子进程替代隔离。

## 时间、资源与清理边界

r2 从 spawn 到目录拒绝为 **0.0984 s**，不是到 ready 的启动时间；本轮没有 ready 测量。其余基础错误约 **0.086–0.095 s**（宿主直接拒绝未知模式 0.0025 s），注入 mode timeout 约 **1.0783 s**。这些是一次本机热缓存样本，不能推出用户体验或并发容量。

采样间隔 100 ms，短生命周期实例只观测到 0–1280 KiB，明显漏过 Peri 实际峰值，**不可作为内存需求**。较长的 mode-timeout 样本观测进程树峰值 **44,544 KiB（43.5 MiB）**，仍是下界性质的抽样，不包含宿主 driver/stub，不是强制上限。原生编译与运行资源观测不能混用。

没有设置 cgroup memory/CPU/PID/总磁盘限制。relay 的 select timeout 是轮询而非 idle deadline，连接/线程无独立并发上限，sendall 可阻塞；nofile128不能补足这些限制。仅用于本轮受控合成调查，不能称为生产资源隔离。候选并发、内存、TTL预算仍待重新测量。

所有本任务验证实例/driver/stub 已停止；原始下载、构建产物、合成 session DB 和日志保留在 `/tmp/helixweave-peri-a0-20261109`，未删除其他任务资源。进程退出检查的范围是该运行创建及观测过的进程，不宣称扫描整个宿主证明无任意未知进程。

## #210 合入后的接续交付核验

现有调查已足以判定**该固定版本和配置接法不通过**，本次不重跑构建、六项宿主测试
或 r1/r2 实验，不通过放宽白名单取得 ready。原始 r1 的无效通过记录、r2 的失败和
未测项均保留。核对 binary、Cargo.lock、构建记录及全部 12 个 Python 文件摘要与
原运行一致；机器可读复核记录见 [post210-review.json](post210-review.json)。

#210 只解除平台取消修复依赖：请求与 acknowledgement 分离、原截止时间内等待
SIGSTOP 可见、证明失败拒绝确认。历史失败的具体 unconfirmed_reason 仍未知。
该平台改动不触及本原型，不能将本报告标为在新 main 执行，也不授权 A1 或写操作。
未测项是下一种接法进入 A1 前必须补齐的门槛，不是本轮可冒充通过的结果。

## 推荐下一步（交主会话决定新实验范围）

先调整接法并重跑 A0，不进入 A1：优先研究上游小范围**完整工具目录过滤**接口，须同时作用于 direct、deferred、资源发现和实际 dispatch；不能只隐藏 provider schema。其次比较使用官方 Rust 组件的最小自定义宿主装配，其维护成本/版本耦合更高。继续采用每活动会话独立 OS 隔离有意义，但纯进程分离不能修复目录问题；TS/WASM 接法也不能未经研究就视作安全替代。当前不决定 fork、不提出平台数据库替换。

若未来补丁获准，原版与补丁版必须分别固定源码/binary/diff，重验全部 A—F，包括本轮未测硬门；不能把补丁通过归给本次原版。产品方面，一次 submit 批准覆盖 create/preflight/start 仍待选择，A1与B阶段授权均未取得。

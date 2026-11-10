# A0 收尾修订：native 矩阵通过，完整验收暂不通过

日期：2026-11-09（Asia/Shanghai）。**调查交付完成；A0 完整验收暂不通过，不建议现在进入 A1。** 固定受限 native 行为矩阵23/23通过、目录库回归通过；新增 ACP Rust 配置单测两次编译被4GiB上限终止，尚未执行。该剩余阻塞属于测试构建环境，不是已发现的运行权限绕过或需要广泛改造的新架构问题。没有接入 A1、真实平台或模型，没有决定长期 fork。

## 身份与原证据

| 对象 | 身份/结论 |
| --- | --- |
| HelixWeave 实验 | `622778f64e59f935d822d92a91bade78807dabce`，原设计 worktree 未迁移；平台已验收 main `b2f15868c6f74ab8764a73ae6fe4394f67304475` 仅更新 #210 依赖，不能改标本实验 |
| 原版 r1/r2 | Peri `c7d154c9985cbe08ba25e964cb7468d15014e6d0`；r1无效通过、r2有效阻断，见 [原报告](../../peri_a0/evidence/RESULTS.md)，均保留 |
| q2 | 同 commit 加三文件 patch；binary `f715e1548bd167ed93b8e761f771d495d57fbc63aaf2d9b9d039e37c2b0cbbff`；工具面首次闭合，但合法 prediction 请求未纳入契约，A0失败；不是泄漏证据 |
| 本轮初次构建 | 同 commit；[build-final](build-final.json)、[8文件补丁](candidate-final.patch)，patch `57d6e37054bfc187134524866840ae890349628c3d822001ecf971e2228fb712` |
| 最终测试桩修正候选 | [build-testfix](build-testfix.json)、[完整patch](candidate-testfix.patch)，patch `06272e1a66423a9e562509a38d044cd8506a1c10e61ef1fd3fe239107680587b`；只补 `#[cfg(test)]` 模块的 wire 身份，生产 diff 不变，见 [同一性记录](testfix-identity.json) |
| 两次构建 binary | 均为 `60144bbe0919992b8913b9031696b2e30b1e9f82d95048ae3b07c887d5d85517`；原产物另存 `/tmp/hw-peri-a0-close-20261109/artifacts/initial/peri` |

Rust/Cargo1.99.0；Cargo.lock `0692683359fb75073eca3d14b364071bc02058be7b11e0335ccaf2b0f314f5dd`。官方 wrapper、`--locked`、离线构建，私有复制的 cargo/target，无共享可写 target。完整 q2→本轮增量见 [patch](q2-to-closeout.patch)。原 r1/r2/q1/q2 的133个文件由 [原始清单](prior-inventory.json) 和收尾检查逐一验证不变。

`f2` 是实际运行代码快照，记录当时 `build-final.json` 与所有 .py 摘要；之后只有构建 driver 支持 testfix 标签、测试 driver 可选择下载锁定缺失依赖和命令级单包测试优化设置，以及 Rust 测试桩修正。没有重标 f2 的运行身份。最终工具链重建得到同一 binary，核心运行代码摘要与 f2一致，因此不重复完整 native 矩阵；最终文件摘要另列收尾检查。

## 改动与实际调用路径

继续选择窄实验补丁 C，没有重新设计平台。上一轮 A/B/C 的接口、隔离及维护比较仍适用（[对照](../../peri_a0_revision/README.md)）。本轮新增启动策略、只读策略回执、完整身份验证和测试；没有重新实现 loop 或权限系统。8个Rust文件不是长期fork决定。维护成本仍是固定内部 catalog/ACP 接口耦合，每次上游升级必须重新核验装配/出口。

- `PERI_PREDICTION_POLICY=off` 在 stdio 入口先校验并冻结。`spawn_prediction` 先统计尝试，再在 disabled 分支返回，位置早于 clone、spawn、上下文投影和模型调用。未知字符串启动失败。没有开放 ACP setter；用户/模型只可走宿主窄 prompt 面。
- `peri/host_policy` 是只读配置/计数回执，不是已验证实际目录的替代品。完整匹配发生 `stage_builder/tools.rs → enforce_strict_pinning → verify_complete`；资格 prompt 成功才证明装配已走通。
- 只有成功普通 turn 的 `prompt_dispatch` 才调用 prediction；固定资格回合、普通成功、审批拒绝、工具错误后成功结束均实际触发 suppressed，未通过 gate 的启动尝试数一直为0。q2 的第四请求与本轮完整 capture 构成关闭前后对照；不是靠一段静默等待。
- 其他可达辅助出口：自动 compact 由既有 `compact.auto_compact_enabled=false` 关闭，手动 `/compact` 被宿主 slash 入口拒绝；Auto classifier 只在 Auto 模式，实验始终 Default；goal、hook、subagent、workflow 受关闭的 middleware、窄 ACP 和工具目录约束。具体固定调用路径见 [源码审阅](source-review.md)。没有承诺所有未来上游辅助出口天然受此开关控制。
- 完整策略包含 canonical/server/wire/namespace/builtin/aliases；fixture schema 精确匹配。合法但不存在的项阻止装配，refresh/replace 不发布不完整目录，pin 后再核验。meta 工具 schema 不通过该字段锁定，仍受固定 binary 和实际模型目录验证。
- 宿主保持 blocked，先核对模式、配置、连接、工具；执行标记的内部合成资格回合并观察真正拒绝，等回合结束的 suppressed 计数满足后才 ready。运行期每回合前后及每次 tools/list 发布前重查，目录描述/参数变化也拒绝。

## 完整验收矩阵

以下分清 native、宿主、OS 与库测试；不要求每层重复同一断言。f2 的23项 case 是分组场景数，OS断言、Rust或unit数量不重复加成“端到端用例数”。

| 门 | 状态 | 实际证据与边界 |
| --- | --- | --- |
| 身份、旧证据不覆盖 | 通过 | 两次build、patch、binary、prototype摘要；[收尾核验](checks.json)，133旧文件未变 |
| A：初始化前禁止用户 prompt；未知模式、设置错误/超时、缺失/额外/替换工具 | 通过 | [f2](f2/report.json) A六项；失败时0模型请求、0工具到达。未知mode在宿主拒绝；mode错误由native拒绝；丢帧由宿主deadline关闭 |
| A：未知合法工具、错误schema、未知字段、prediction错误配置 | 通过 | f2 A_policy四项；均未ready、0 provider、0工具 |
| A：完整ready及prediction关闭 | 通过 | A_prediction_paths；成功/工具错误正常结束，3次suppressed、started0；资格与每回合检查；全25实例128请求无stub_failure |
| A/B：Default拒绝和允许成功 | 通过 | A_B_permission_tools：新增拒绝审批后tools/call0；允许increment1、lookup1。不是broker最后一层拒绝伪装Peri审批成功 |
| B：伪造工具/资源/元工具/扩权 | 通过受测闭合面 | 真正模型输出直接或Execute转发；native返回`Tool not found`；搜索仅返回两个fixture，total_available2；无未批准工具效果。后台资源RPC仍由broker拒绝，不称Peri不发任何resources请求 |
| B：越界/未知参数 | 通过服务端边界 | 错参数可能到达fixture tools/call，broker拒绝且counter/lookup不变；不宣称native在dispatch前校验JSON schema |
| B：目录动态变化/身份变更 | 通过分层契约 | B_changed_catalog回合前关闭；[host unit](host-unit-final.log) 活动tools/list返回前拦截description变更；[native内部generation补测14/0](rust-agent-testfix/report.json)，不是wire推送刷新全链重演 |
| A：缺broker、审批断连/超时 | 通过当前接法 | A_missing_broker是真实ACP回调返回method-not-found；A_approval_disconnect/timeout由宿主关闭，均新审批后tools/call0、效果0，不回退Bypass；不代表Rust None broker构造或native自主超时 |
| C：OS边界 | 通过受测范围 | f2/report 12断言；实际2次open、3次connect攻击假秘密/其他目录/TCP/pathname及abstract socket，其他env/FD/proc/caps/private write检查。探针使用相同sandbox启动参数，非仅未注册工具 |
| C：两个native独立实例 | 通过 | C_independent_instances交错；双方Default，A允许/B拒绝，最终各1计数；session、标记上下文、回调/代次不串用；A目录变化关闭不影响B。未修改成Bypass，也不证明同进程多会话 |
| D：重启及旧代次 | 通过 | D_restart kill后旧连接closed，新实例先blocked、重新资格，旧session/instance/callback拒绝，新代次正常；load/resume/fork始终不暴露，未调用激活历史方法 |
| E：完整合成外发 | 通过封闭字段投影 | E_minimum_egress两次真实lookup：成功和isError结果经投影进入后续provider请求；假秘密/路径/标识符/指令标记均不存在，稳定INPUT_INVALID存在。检查范围为整个E会话capture；B场景中stub主动产生的攻击参数可进入其自身历史，不混算为宿主数据泄漏。不是正则万能过滤，不代表任意自由文本保密 |
| E：个人配置/遥测/非授权网络 | 通过配置与OS限定 | 环境/挂载白名单；固定禁用hooks/skills及其他MCP；唯一模型UDS。未作全面syscall/packet流量统计，不声称遥测代码不存在 |
| F：协议坏帧、模型超时、工具回调失联 | 通过宿主关闭契约 | F_malformed向真实使用的解析函数注入坏帧并关闭活native（不是native吐坏JSON）；F_timeout实际模型无响应；F_callback_disconnect真实工具回调失联，均closed无副作用 |
| F：native资源故障 | 通过受限FD故障，非OOM验收 | F_native_resource定位实际native PID后注入soft nofile0。f1走RPC_TIMEOUT；f2实际native -32010 `store is unavailable: session database rejected the write`，provider_delta0，宿主关闭清理。没有errno/EMFILE证据，不推断具体open/connect失败 |
| F：传输与总资源限制 | 通过既定边界 | [4项真实socket测试](runtime-tests.json)；运行scope实际512MiB/swap0/1核/128任务。之前控制器OOM/PID/CPU触发证据保留，未当成native资源结果 |
| F：清理 | 通过 | f2 all_processes_exited=true；25实例cleanup_errors均空、观测PID均退出；[最终清理](cleanup.json)：244个观测PID均不存在，无任务socket、无活scope；OOM失败元数据保留 |

明确未测且不外推：同native进程多会话隔离；native load/resume/fork行为（本接法禁用）；Native cgroup OOM恢复；wire主动目录推送刷新全链（内部generation和宿主发布边界分别测试）；任意用户自由文本/大规模长会话/跨宿主攻击；真实平台与全部A1/B阶段。这些不作为已完成能力。基础幂等、确认变更/过期、撤销、重复/丢响应仍是未来B2开放写操作前置。

## 中间失败没有覆盖

- [host-unit.log](host-unit.log)：新增fixture字段未在测试桩初始化，1 error；修正后7项通过。没有弱化行为断言。
- [t1](t1/report.json)：三项目录负例已经在更早的LIVE_CATALOG_MISMATCH处安全拒绝，但测试仍预期旧码，所以整体failed；更新准确预期、加强native Tool not found和deferred结果断言后，[t2](t2/report.json)定向通过。
- [f1](f1/report.json)完整23项通过；源码独立审阅指出活动tools/list还要发布前核验，补齐回调及unit后，[f2](f2/report.json)完整23项再通过。
- 初次Rust目录测试13通过/1失败：BridgeTool桩server有值但wire缺失，严格parser拒绝，根本没执行新增generation断言；不能归为生产filter失败或把f2文字当补测已通过。原测试源码、[失败log](rust-agent/session::tool_catalog::tests.log)及[补测14/0](rust-agent-testfix/report.json)分别保存。
- [ACP Rust首次](rust-acp/report.json)因离线缓存缺少assert-json-diff锁定测试依赖而未运行；只补取本任务缓存后重试，默认优化编译触发4GiB scope OOM，未执行测试且driver未及写report，见 [退出记录](rust-acp-online/termination.json)。之后两次低优化命令因Cargo --config参数层次错误使依赖patch未生效，在解析依赖阶段退出（[离线](rust-acp-lowopt/report.json)、[联网](rust-acp-lowopt-online/report.json)）；已修正为同一全局参数层次，不改源码、lock或限额。不能把这些记录算单测通过。修正参数后单包opt-level=0仍触发4GiB OOM，见 [最终终止记录](rust-acp-lowopt-fixed/termination.json)，测试尚未开始，停止进一步重试；运行限额及补丁未放宽。
- q2的合法prediction辅助请求仍是原负结果；未声称其是泄漏，未放行任意无tools请求。

## 测量、维护与交接

[汇总](summary.json)：f2有25个实例，其中15个进入ready；单次热缓存合成启动范围0.1402—0.2649s（包含固定资格探针，不是生产SLA）。100ms后代树RSS抽样最高52480KiB（51.25MiB），不含host stub，可能漏峰或重复共享页；不等于scope总内存。没有长会话容量或多用户并发结论。

初次增量构建543.96s，不含缓存复制/下载；测试桩修正重建约0.85s，binary不变。构建scope4GiB/4核/512任务与运行scope分开。没有新增总磁盘quota或seccomp保证。保留任务本地缓存、产物、合成会话和原始证据；交付目录只复制必要JSON/log/patch，不提交session数据库或凭据。

定向测试归档：Rust策略2项+既有filter2项+目录补测14项，共18个不同用例通过；初次目录13/1不算通过，重复13个不重复计数。prediction配置单测未运行，不能计入19项通过。宿主contract7项、transport4项通过；Ruff、Rustfmt、相对链接、完整新增文件空白和git diff --check通过，见[检查](checks.json)。未跑完整CI、科学或产品浏览器/Bulk Gate。

下一步最小选择：单独评估并授权ACP测试编译预算，保持native隔离和512MiB运行限额不变；或暂缓接入。没有证据要求现在扩大业务补丁、改用通用模型网关或决定长期fork。完成这一剩余验证后再审阅是否进入A1，本轮不自动实施。

复现命令见 [README](../README.md)。实际参与者为Codex主agent及两个只读审阅agent，主agent唯一写入/构建；没有commit/push/PR。尚待决定：是否授权独立A1工作、是否接受固定内部接口实验补丁维护或另开官方组件装配评估。一次submit确认覆盖create/preflight/start仍是待用户确认的产品建议，不在本轮实现。

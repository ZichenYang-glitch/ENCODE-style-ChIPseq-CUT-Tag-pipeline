# 接法修订结果：候选局部有效，A0 仍不通过

日期：2026-11-09（Asia/Shanghai）。本轮选择 **c7d154… + 三文件完整目录策略实验补丁**。首次工具面和一次 Default 审批拒绝已验证，但成功 turn 后的自动 prediction 模型调用不受工具目录策略约束，严格 stub 拒绝后关闭实例。**不扩大白名单或补丁范围，不进入 A1。** 这是接法仍需调整的调查结论，不是宣称 Peri 存在越权漏洞。

## 身份及原版/候选对照

| 对象 | 固定身份与结果 |
| --- | --- |
| HelixWeave | 实验 worktree 保持 `622778f64e59f935d822d92a91bade78807dabce`；平台已验收 main 为 `b2f15868c6f74ab8764a73ae6fe4394f67304475`，没有迁移或将旧结果标为新main运行 |
| 原版 | Peri `c7d154c9985cbe08ba25e964cb7468d15014e6d0`，binary `ecb4eb1394c67db2c57bc74373f194d52c6342e3be425ed0980b665ec29d626f`；[原r2](../../peri_a0/evidence/RESULTS.md) 工具闭合失败；r1无效通过仍保留，不重写 |
| 官方当前比较对象 | `3288b07bf871276026ae0677129821b096e75118`，557文件增量；关键catalog、顶层allow-all、Search/Execute未解决本问题。仅只读比较，未构建此版本。见 [比较记录](upstream-comparison.json) |
| 首个补丁候选 q1 | base同原版；binary `cf6644fe2a67ad9ad9c377230cd997855fc34246e22ce420a14946f97c43a6a6`；[构建](build.json)、[diff](candidate-q1.patch)、[运行](q1/report.json)。归档binary在任务根 `artifacts/candidate-q1/peri`，构建manifest中的临时target随后由最终构建替换，不能再用该target冒充q1 |
| 最终候选 q2 | base同原版；binary `f715e1548bd167ed93b8e761f771d495d57fbc63aaf2d9b9d039e37c2b0cbbff`；[构建](build-final.json)、[完整diff](candidate-final.patch)、[运行](q2/report.json)。patch SHA256 `af1bad255c0644c9b0bdf377294f9a5d7c83829b3cde9b9065e1092792afbfba` |

Rust/Cargo 1.99.0，锁文件摘要保持 `0692683359fb75073eca3d14b364071bc02058be7b11e0335ccaf2b0f314f5dd`。复用原工具链只读文件和独立复制的依赖/target缓存；官方 wrapper、`--locked`、离线构建，未改全局配置。peri-cool 是原源码中的文档站 gitlink，不参与候选 Rust 构建；新候选 checkout 未初始化文档站子模块，manifest 继承的 submodule_gitlink 仅记录 Git 对象身份，不声称进行了其工作树构建。

q1→q2 只收紧补丁兼容性（strict pin 仅在显式策略开启时应用）和 Rust 格式；未增加能力。Python stub 把 q1 的误称 `malformed_provider_body` 修正为准确的 `UNAPPROVED_MODEL_REQUEST_CLASS`；原q1记录不改。运行记录包含各自原型 SHA。q2之后仅修正 Rust 验证driver的零用例检查；其旧源码另存 [driver原稿](verify_patch-before-zero-test-fix.py)，实际运行核心未变。对应核验见 [checks.json](checks.json)。

## 能力与停止原因

[方案比较及策略契约](../README.md) 给出 A/B/C 接口、成本和选择依据。三项原额外能力的精确边界：

| 能力 | 固定源码事实 | 本轮限制 |
| --- | --- | --- |
| ReadContextArchive | [wrapper](https://github.com/KonghaYao/peri/blob/c7d154c9985cbe08ba25e964cb7468d15014e6d0/peri-middlewares/src/mcp/context_archive.rs#L27) 从可信session_id读取canonical或登记的inherited archive；[reader](https://github.com/KonghaYao/peri/blob/c7d154c9985cbe08ba25e964cb7468d15014e6d0/peri-agent/src/context_archive.rs#L137)校验anchor/digest | 非任意文件工具，但历史原文可能绕过平台字段投影重新进入模型；不批准 |
| DiscoverMCP | [实现](https://github.com/KonghaYao/peri/blob/c7d154c9985cbe08ba25e964cb7468d15014e6d0/peri-middlewares/src/mcp/discover_tool.rs#L152)读目录、pool及registry缓存 | 不执行工具；可返回URL/schema/资源URI/描述，不能视作只有批准字段；不批准 |
| mcp_read_resource | [实现](https://github.com/KonghaYao/peri/blob/c7d154c9985cbe08ba25e964cb7468d15014e6d0/peri-middlewares/src/mcp/resource_tool.rs#L194)绑定session可见server、generation并调用resources/read；可缓存或刷新registry | 非任意本地open，但可触发MCP I/O；资源URI不是HelixWeave白名单。description构造还可能取pool-wide摘要；无本轮跨实例泄漏证明；不批准 |

三者 intrinsic 身份不等于 external MCP；Default 不是普遍审批所有只读能力。原两个fixture则因绑定的external MCP来源进入审批。候选在工具对象目录过滤，不删除工具实现，不改prompt来隐藏能力；后台resources/templates/list仍会尝试，broker拒绝并计数，不能称过滤补丁消除了该RPC。

q2 [资格证据](q2/A_B_permission_tools/permissions/evidence.json)：前三个provider请求只有 SearchExtraTools/ExecuteExtraTool；deferred搜索只返回fixture，Default发出一次increment审批，被宿主reject_once，**tools/call到达0、副作用0**。这次实际拒绝不依赖broker最终拦截。

第四个provider请求不带tools，system为 `<prediction_directive>`，携带刚才的合成对话及工具结果，准备预测输入并维护标题/标签/摘要。[prompt_dispatch](https://github.com/KonghaYao/peri/blob/c7d154c9985cbe08ba25e964cb7468d15014e6d0/peri-acp/src/host/prompt_dispatch.rs#L471)在成功普通turn后启动；[prediction](https://github.com/KonghaYao/peri/blob/c7d154c9985cbe08ba25e964cb7468d15014e6d0/peri-acp/src/host/prediction.rs#L74)直接调用模型，`caps.prediction`只控制后续通知，并非模型调用开关。未核实到本候选的现成关闭入口。stub没有为该请求返回脚本响应，随后native被关闭。

原型曾在0.1154s设置ready，约0.1166s因该请求关闭；因此“主模型资格通过即可ready”也不构成完整准入。该辅助路径是 **上游native装配/模型出口契约限制**，同时暴露原stub把全部请求当作工具推理请求的接线假设。请求本身合法，不能称协议坏消息；也没有真实秘密或真实模型外发。若下一轮继续，应先明确关闭或批准辅助模型路径的契约，不能直接把无tools请求加入允许列表以获得绿灯。

## A—F 完整矩阵

| 验收 | 状态 | 证据及边界 |
| --- | --- | --- |
| 身份、原始证据不覆盖 | 通过 | 构建/patch/原型SHA；checks核对原目录53文件；原r1/r2没有变化 |
| A：blocked下拒prompt；模式错误/超时/未知；工具缺失/额外/替换 | 六场景通过 | q2/report.json；错session导致native拒绝、driver丢mode帧的1s超时、宿主未知mode拒绝分别记录；无模型请求或工具到达 |
| A：完整ready准入 | **失败** | 成功主turn后额外prediction请求触发关闭；一次短暂ready不是全部模型路径已准入 |
| B：首次provider目录及一次deferred搜索 | 部分通过 | 两个预期meta；fixture检索结果；没有原三项intrinsic广告。未穷尽目录刷新/所有后续发现 |
| A/B：Default审批拒绝 | 一次实际拒绝通过 | permission request 1，broker tools/call 0，counter/lookup 0；非最后一层拒绝掩盖 |
| B：允许工具成功、伪造工具/扩权、参数越界 | native未测；宿主片段通过 | q2资格整体失败，尚未到允许成功和伪造调用循环；host-unit 6测试不冒充native；场景保留在cases.py |
| B：资源、alias/namespace、元工具及实际dispatch闭合 | 部分源码/单元证据；完整未测 | 2个新policy测试、2个现有filter测试、13个catalog测试；禁止目标同名异来源/严格pin解析拒绝已测。动态更新和refresh错误新分支、完整alias/namespace矩阵尚未实测 |
| B：未知配置/未知能力 | 部分通过，有缺口 | JSON未知字段/版本/通配拒绝；合法但不存在的allow身份尚不立即报错。未知配置不能泛称全部已闭合 |
| A：审批断连/超时/缺broker不Bypass | native未测 | 因候选资格未整体通过而停止；源码保证不替代实验 |
| C：OS实际访问探针 | 12断言通过 | q2/report.json：2个open与3个connect针对假秘密/另一假目录/宿主TCP及pathname/abstract socket；其余env、FD、proc、cap、私有写入。不是12个访问攻击，no_host_mounts也非全面挂载审计 |
| C：两个native实例交错配置/权限/目录/回调隔离 | 未测 | 单实例OS隔离不能代替双实例交错；同进程多会话仍为未来研究 |
| D：重启/旧代次/native恢复拒绝 | native未测；宿主片段通过 | 仅暴露session/new；跨实例、旧代次、任意ACP拒绝的host-unit保留，不宣称完整restart通过；未调用load/resume/fork |
| E：最小外发/提示注入投影 | 宿主片段通过，完整provider未测 | host-unit封闭字段投影；没有执行E_minimum_egress；捕获辅助请求本身表明还需覆盖额外出口。别名不等于匿名化 |
| E：个人配置、hooks、skills、遥测与网络 | OS/config片段通过 | env白名单、私有mount、唯一UDS出口；未读个人配置；无packet/syscall全面外联统计，不宣称遥测代码不存在 |
| F：relay/stub期限、断连、连接/体积上限 | 4个宿主行为测试通过 | [runtime-tests.log](runtime-tests.log)；真实socket和背压，不靠观察RSS判断资源限额 |
| F：总内存/CPU/PID限制 | 控制器触发证据通过 | [OOM scope](resource-oom-systemd.txt)：64MiB触发oom-kill/stop；[PID/CPU](resource-pids-cpu.json)：pids.max拒绝、CPU限流递增、子进程回收。触发对象是受控探针，非native运行中故障验收 |
| F：坏协议/模型超时/工具回调失联/native资源超限 | native未测 | 保留场景，资格失败后不放宽权限继续；宿主断连测试不替代这些native门 |
| F：已运行路径失败与清理 | 通过该片段 | q2 passed=false/all_processes_exited=true；[cleanup](cleanup.json)核对38个观测PID均不存在、无模型socket，所有scope无ControlGroup；OOM failed元数据保留，无存活服务 |

## 测试、测量与证据限制

- Rust：新增policy 2 passed、既有filter 2 passed、catalog 13 passed，合计17个不同用例，0失败。初次 `tool_catalog_test` 过滤器匹配0，不算通过；首次 `patch-tests/report.json` 的总passed因此无效。修订验证driver时有一次NameError，保留该目录，不算验收；最终 [catalog补测](patch-tests-catalog2/report.json) 为13/0/0。重复执行的2+2不重复计数。
- 宿主contract 6 passed、transport 4 passed。Rustfmt、Ruff check/format、链接、新增源码/文档空白及git diff --check见checks。未运行科学任务、完整CI、真实平台端到端或Bulk Gate。
- q1构建314.74s；最终增量构建223.03s，包含复制缓存后的重编译，不含下载/缓存复制。编译使用4GiB/4核/512任务scope，不能混作运行容量。
- q2资格首次ready约0.1154s、拒绝关闭约0.1166s；单次热缓存样本，**不是可用会话启动SLA**。100ms进程树RSS抽样峰值46080KiB（45MiB），不含宿主stub，可能重复共享页/漏峰；与实际512MiB cgroup硬限额分开报告。没有成功长会话或并发容量测量。
- OOM探针因scope stop未写出自身JSON；systemd记录提供退出证据。PID/CPU探针独立执行，不能把二者合写为一次完整probe通过。总磁盘quota和seccomp未实现；没有借nofile替代总内存/PID/CPU限制。

进程已回收；仅保留本任务缓存、二进制、合成session DB、原始日志与证据，未清理其他任务资源。可重复命令见 [README](../README.md)。q1/q2负结果、初次format工具启动缺库、测试driver错误均保留，不以最终通过的局部测试覆盖历史失败。

## 交回审阅的选择

本轮选型、最小补丁和复验已交付，**goal完成不等于A0通过**。后续如授权继续，先研究显式关闭辅助prediction的窄入口，或重新比较能控制全部模型调用的官方最小宿主；同时补齐策略未知身份拒绝和所有未测硬门。不能默认放行prediction，也不能认定TS/WASM天然更安全。暂不决定fork，不向上游发消息。

若产品确需ReadContextArchive、预测输入或标题摘要，应单独决定消费者、外发投影、会话绑定及中断/升级契约，再验收；本轮不批准这些能力。A1和B阶段仍未授权，B2基础幂等/确认失效/撤销/丢响应门不后移，一次确认覆盖create/preflight/start仍待用户决定。

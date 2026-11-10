# 固定候选调用路径审阅

2026-11-09；Codex只读审阅agent，主agent核对并保存。基线c7d154…及本目录candidate-testfix.patch。没有构建、修改或运行实验；测试由主agent执行。

| 路径（补丁后的源码行） | 裁定 |
| --- | --- |
| peri-acp/src/host/stdio/mod.rs `run_acp_stdio` | 开始时读取/验证prediction和工具策略；进程OnceLock冻结。配置错误先失败 |
| peri-acp/src/host/prompt_dispatch.rs:472 | 仅`!is_continuation && result.is_ok()`调用spawn_prediction；错误/取消路径不应被误计为成功回合 |
| peri-acp/src/host/prediction.rs:61–65 → :74 → :86 → :93 | disabled先return，后面才clone、task spawn、读取有效历史投影；started在spawn前递增，0表示未越过gate |
| peri-acp/src/host/requests.rs `peri/host_policy` | 只读空参数回执，无setter；宿主用户方法面不暴露此ACP方法，也不让模型直接发ACP |
| peri-agent/src/session/exec/stage_builder/tools.rs:84 | 初始enforce_strict_pinning才验证实际目录；不能用策略回执替代该调用成功 |
| peri-agent/src/session/tool_catalog.rs:447、530、564 | refresh发布前、startup静态替换发布前、pin后完整identity核验；530/532重复验证无害冗余，本轮不做顺手清理 |
| peri-agent/src/agent/stages/reason.rs:37 | strict refresh失败返回，不继续使用旧目录发模型请求 |

其他模型入口：

- `peri-acp/src/host/prompt.rs:493`从配置读取compact，:498应用环境覆盖；sandbox清空环境，没有compact覆盖变量。`peri-agent/src/agent/stages/compact.rs:84`检查disabled，`compact_v2/mod.rs:278`非force同样检查。
- 手动compact可强制：`peri-agent/src/session/exec/executor_helpers/intercept.rs`处理slash命令，因此宿主拒绝所有leading slash，不能仅关闭auto后声称手动也关了。
- `peri-middlewares/src/permission/mod.rs:603`的Auto classifier只在AutoMode；本接法固定Default。
- Goal/Hook/SubAgent/Workflow middleware关闭，工具及ACP入口亦收窄。未来重新开启任何一项需重新审模型出口；没有建立通用模型网关。

发现并修正的宿主缺口：最初仅回合前后比较`declared_tools`，活动回合tools/list可能先发布description变化。最终回调先比较完整固定TOOLS再返回；host-unit活动回合测试证明拒绝时不返回恶意描述、实例closed、无grant/效果。Rust策略本身不锁description，因此这层宿主约束是必需的。

未发现必须扩大本轮内部架构补丁的其他缺口。该结论限定固定版本、受控宿主和当前关闭能力，不是未来Peri版本或所有协议面的证明。

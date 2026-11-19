# Peri 活动模型绑定（IRENE-44）

基线为主干的 A1 受限启动器和 [P2 离线授权期](peri-p2-model-transport.md)。
外发门按 [#247 的设计](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/blob/1b055fd33183c79d101055cb2afcf39d1e6d8e6d/docs/architecture/peri-p2-vendor-outbound.md)。
HelixWeave 负责业务授权、投影、额度和采用，Peri 负责供应商绑定及协议。

## 固定版本的观察

用产品 pin `5cc82caf6e04578f19ff5790fbc4eb8c9e115793ff6950a699c9e2ad28c50f47`
（Peri `c7d154c9985cbe08ba25e964cb7468d15014e6d0` 及已资格验证的补丁）经现有
Bubblewrap、私有 Unix socket 和 systemd scope 启动。输入只来自临时合成 Peri
设置。没有访问用户的配置目录。

- 协商 `peri.sessionWorkspaceV1` 后，`session/metadata` 的模型相关字段为
  `modelAlias`、`modelName`、`providerName` 和 `effort`。读到的是显示值。
  `session/new.configOptions` 和 `config_option_update` 也只投影显示选项。
  `peri/host_policy` 回报工具、预测和压缩策略，未提供模型绑定代次。
- 使用受限启动器的 `--settings-stdin` 时，Peri 对有效模型别名的持久化切换
  返回 `-32603`（Failed to persist config）。探针用新的 Peri 进程激活新设置，
  A/B 可以保持相同显示值。使用不存在的活动别名时，进程在 ACP 初始化前退出，
  同一配置里的其他已配置档位没有收到请求。
- `session/prompt` 没有 expected binding 参数。该实现为执行捕获 provider
  快照，显示值的前后两次读取无法把请求和快照绑定。ACP 的 `session/config_update`
  也不能作为确认回执。它是写通知，无法提供一次交换的预期代次条件。

参考固定源码的
[模型选项](https://github.com/KonghaYao/peri/blob/c7d154c9985cbe08ba25e964cb7468d15014e6d0/peri-acp/src/session/state_builders.rs)、
[会话元数据](https://github.com/KonghaYao/peri/blob/c7d154c9985cbe08ba25e964cb7468d15014e6d0/peri-acp/src/host/requests/session_restore.rs)、
[注入设置权限](https://github.com/KonghaYao/peri/blob/c7d154c9985cbe08ba25e964cb7468d15014e6d0/peri-config/src/settings.rs)和
[执行快照](https://github.com/KonghaYao/peri/blob/c7d154c9985cbe08ba25e964cb7468d15014e6d0/peri-acp/src/host/prompt.rs)。

## 接入所缺的最小契约

1. 受信任 Peri 部署提供只读的可用状态与不透明绑定代次。有效供应商、模型、
   端点、凭据或影响授权的配置变化时，代次推进，同名换绑及 A→B→A 也适用。
   展示、授权和回执仅需要该代次，无需读取 Peri 设置。
2. 发送接受预期代次。Peri 在建立外发前用该代次核对并冻结同一个 provider，
   通过有界的单次交换返回所用代次或稳定失败码。HelixWeave 将许可冻结在该代次
   上，并在响应和采用时重查授权期与绑定。启动子进程的 PID 或平台会话号仅标识
   生命周期。当前 ACP 缺少上述两个操作，真实供应商入口继续返回 `vendor_unselected`。
3. 实际外发须证明业务体仅有 `call_id` 与获准投影，并能证明与扣减既有费用上限。
   常规 ACP prompt 会套用聊天上下文及工具描述。本地验收用合成会话，当前没有
   证明常规 ACP prompt 可作为生产外发通道。

## 本地验收接口

复用 `POST /api/v1/agent/authoring-sessions/{sid}/turns` 及同一路径下的 turn
读取入口。在临时 SQLite 上，通过测试注入的本地 Peri 启动装配运行真实固定进程。
初始化使用现有资格 stub 行为。业务验证固定投影至合成字段 `syn.threads`，
响应经现有 P2 deliver/adopt，再通过真实 authoring 校验生成建议。

- A/B 各在独立的私有 socket 上。通过临时 Peri 设置选中不同端点，保留同一
  供应商 ID、模型名和档位。相同 HTTP 建议请求读取不同的响应值。
- 缺活动绑定时模拟端点收到零请求。换绑后用新的 Peri 进程和新的平台授权期。
  在模拟端点阻塞时换绑或撤权，放行迟到响应后，旧请求仍没有提案或可用回执。
- 测试的允许集固定在本地授权期，费用检查只用合成额度，发送以后不退预留。
  不合格的字段或额度在业务外发前被拒绝。回执仅含合成调用标识、字段名和结果状态。
- 检查私有端点标记及合成凭据未出现在 HTTP、平台 SQLite 和日志。验收仅生成
  本地非敏感证据。生产接入须由 Tibo 核对新的 Peri 契约与既有隔离界面。

## 本机复现

先按[运行时资格记录](../development/peri-runtime-qualification/README.md)准备锁定的
`ci-fast` 环境、用户 systemd 会话和固定 Peri 二进制。用本检出运行：

```bash
HELIXWEAVE_PERI_BINARY="$PINNED_PERI" \
  ./.local/envs/ci-fast/bin/python -I -S scripts/checkout_bootstrap.py \
  --repository-root . pytest -m platform_real_execution \
  test/api/test_peri_binding_local.py test/services/test_peri_authoring_bridge.py
```

两套模拟端点分别回读 `config.threads=7` 和 `11`，而请求里的偏好是 `4`。
原生绑定代码与 API 路由参与这条执行路径。测试装配持有合成绑定代次、额度和收据，
并提供模拟端点。换绑后使用新的原生进程和平台连接。保留同一拟议字段与业务入口，
晚到响应无建议，已经发送的合成预留仍消耗，下一次请求额度不足时零发送。

回执在本机测试进程内生成。浏览器可见的 HTTP 响应、临时平台 SQLite 和日志均用
合成端点与凭据标记检查。普通 ACP prompt 的业务消息只含合成投影和关联号，外层
聊天协议仍带有原生协议字段及合成的资格上下文。该路径的运行装配只存在于测试中。
真实供应商、生产绑定证明、费用计量、个人设置、研究数据和真实浏览器未在本轮验证。

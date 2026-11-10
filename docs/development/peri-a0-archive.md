# Peri A0 归档

当前状态（2026-11-09，ACP 补测之后）：固定 HelixWeave `622778f64e59f935d822d92a91bade78807dabce`、固定 Peri `c7d154c9985cbe08ba25e964cb7468d15014e6d0`、固定补丁和已声明边界内，A0 验收通过。A1 没有实施。没有真实平台集成、真实模型外发或生产容量结论，也没有决定长期 fork。

这是阅读入口。历史报告里的“暂不通过”保留为当时结论，不在这里改写成通过。

## 读哪一份材料

| 需要的内容 | 位置 | 它证明什么 |
| --- | --- | --- |
| 当前通过范围 | [ACP 补测裁决](../../scripts/peri_a0_closeout/evidence/acp-final/RESULTS.md) | 补测后的验收结论 |
| 补测前的阻塞 | [收尾报告](../../scripts/peri_a0_closeout/evidence/RESULTS.md) | native 23/23 已过，ACP 单测当时因 4GiB OOM 未执行 |
| 原版 r1/r2 | [原报告](../../scripts/peri_a0/evidence/RESULTS.md) | 纯配置接法不通过，保留 |
| q1/q2 | [修订报告](../../scripts/peri_a0_revision/evidence/RESULTS.md) | 工具面有进展，合法 prediction 请求未纳入契约，候选失败 |
| 已接受补丁 | [candidate-testfix.patch](../../scripts/peri_a0_closeout/evidence/candidate-testfix.patch) | 相对 Peri `c7d154…` 的实验补丁，不是官方发布 |
| 补丁许可 | [PATCH-NOTICE](../../scripts/peri_a0_closeout/evidence/PATCH-NOTICE.md) | Apache-2.0 上游，本仓库不内置 Peri 树 |

设计背景在 [ADR](../architecture/peri-agent-integration-decision.md) 和 [契约附录](../architecture/peri-agent-integration-contracts.md)。两份文档里带日期的失败叙述是当时状态。

## 固定身份

| 对象 | 值 |
| --- | --- |
| 实验时的 HelixWeave | `622778f64e59f935d822d92a91bade78807dabce` |
| 归档时的 main | `b2f15868c6f74ab8764a73ae6fe4394f67304475`，只包含 #210，不是这次实验的运行基线 |
| Peri | `c7d154c9985cbe08ba25e964cb7468d15014e6d0` |
| 完整补丁 SHA256 | `06272e1a66423a9e562509a38d044cd8506a1c10e61ef1fd3fe239107680587b` |
| 运行 binary SHA256 | `60144bbe0919992b8913b9031696b2e30b1e9f82d95048ae3b07c887d5d85517` |
| 工具链 | Rust/Cargo 1.99.0，`cargo --locked` |

`candidate-final.patch` 是测试桩修正前的 8 文件补丁。`testfix-identity.json` 记录修正只在 `#[cfg(test)]` 模块，生产 diff 和运行 binary 不变。不要把两次补丁摘要混成一个。

## 已接受的证据范围

- native 23/23，25 个实例、128 次模型请求，见 [f2/report.json](../../scripts/peri_a0_closeout/evidence/f2/report.json)。本归档没有重跑。
- Rust 19 个不同用例通过。其中 ACP 原策略单测 1 passed、0 failed、0 ignored，见 [acp-final](../../scripts/peri_a0_closeout/evidence/acp-final/RESULTS.md)。
- 宿主 contract/transport 共 11 项通过，沿用原证据。
- ACP 编译 scope 为 8GiB、单编译任务，内核峰值 6.2784GiB。native 运行限制仍是 512MiB。编译峰值不是运行需求。

仍未验证：同进程多会话、被禁用的 native 恢复、native cgroup OOM 恢复、wire 目录主动推送全链、任意自由文本保密、大规模长会话、真实平台和 A1。4GiB OOM、错误 Cargo 参数和未执行记录保留在 `rust-acp*` 目录，不能改记为通过。

## 复核现成证据

不需要原来的 `/tmp` 目录。在仓库根目录执行：

```bash
(cd scripts/peri_a0_closeout/evidence && sha256sum -c SHA256SUMS)
(cd scripts/peri_a0_closeout/evidence/acp-final && sha256sum -c SHA256SUMS)
python3 -S scripts/peri_a0_closeout/test_contract.py
```

`test_contract.py` 只检查宿主 broker 单元行为，不代替 native 矩阵。证据 JSON 里的 `/tmp` 和 cgroup 路径是当时运行事实。

本归档没有包含 Peri 源码树、Cargo 缓存、`target`、运行 binary 或合成会话数据库。因此上面的命令不能证明从零构建。

## 重新执行

重新执行使用调用者自己的任务目录，例如 `TASK=/path/to/new-a0-task`。目录必须尚不存在于本次输出用途，并且不要和已有实验共用可写 `target`。

前置条件：

- 只读的 Rust/Cargo 1.99.0。
- 将 Peri `c7d154c9985cbe08ba25e964cb7468d15014e6d0` 克隆到 `"$TASK/peri"`，再应用仓库内的 `scripts/peri_a0_closeout/evidence/candidate-testfix.patch`。
- bwrap，以及能强制 `MemoryMax`、`MemorySwapMax=0`、`CPUQuota`、`TasksMax` 的用户 scope。权限不足时停止，不能改成普通子进程。
- 原验收构建复制了已核验的 Cargo 缓存和 `target`。那些缓存不在本仓库。没有它们时，`--online-locked` 只允许为任务本地补取 `Cargo.lock` 锁定依赖；这不是已验收 binary 的从零构建证明。

native 运行 scope 仍是 512MiB、swap 0、1 核、128 任务：

```bash
systemd-run --user --scope --unit=hw-peri-a0-review \
  --property=MemoryMax=512M --property=MemorySwapMax=0 \
  --property=CPUQuota=100% --property=TasksMax=128 \
  python3 scripts/peri_a0_closeout/run.py \
  --source "$TASK/peri" \
  --build-manifest "$TASK/evidence/build-testfix.json" \
  --evidence "$TASK/run-evidence"
```

ACP 单测若要重跑，使用新的 scope 名和 8GiB、单编译任务、`peri-acp` 测试包 `opt-level=0`。首轮通过后原实验没有再跑非 LTO 备选。编译峰值不是 512MiB 运行限额。

`build_candidate.py` 不传参时仍默认历史 `/tmp` 路径。归档后的复现必须传入 `--task-root`、`--source` 和 `--reuse-root`。Python 原型脚本没有为这次归档改运行逻辑。新增的是本页、补丁声明，以及设计文档和收尾说明里的当前状态指针。`pyproject.toml` 把这三个原型目录排除在 Ruff 之外，以保留与 2026-11-09 证据摘要一致的脚本字节。

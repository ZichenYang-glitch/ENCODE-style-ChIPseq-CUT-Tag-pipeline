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
| 搬迁预检 | [plan.py](../../scripts/peri_a0_reproduce/plan.py) | 用当前工具链、缓存和任务目录构造计划，不启动构建 |

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

本归档没有包含 Peri 源码树、Cargo 缓存、`target`、运行 binary 或合成会话数据库。因此上面的命令不能证明从零构建。下面的搬迁预检同样不启动构建。

## 重新执行

历史脚本、补丁和证据保持实验当时的字节。参数不能跨入口混用。

`scripts/peri_a0_closeout/build_candidate.py` 接受 `--task-root`、`--source`、`--reuse-root`、`--label`。它读取历史 manifest 的 `environment`，只替换 `HOME`、`CARGO_HOME`、`CARGO_TARGET_DIR` 中能对上 `--reuse-root` 的前缀，不重定位 `PATH`。它不接受 `--online-locked`。不传参时默认目录仍是 2026-11-09 的 `/tmp` 路径。它不是搬迁后的复现入口。

`--online-locked` 只属于 `scripts/peri_a0_closeout/verify_patch.py`，以及证据目录里冻结的同名 driver 副本。打开后，该测试 driver 把 manifest 环境中的 `CARGO_NET_OFFLINE` 改成 `false` 并启动测试。它不构造新的工具链或缓存路径，也不是从零构建开关。`run.py` 只接受 `--build-manifest`、`--source`、`--evidence`、`--qualification-only`、`--targeted`，同样没有这个参数。

搬迁后的预检是 [plan.py](../../scripts/peri_a0_reproduce/plan.py)。它不启动 Cargo，不复制缓存，不应用补丁。调用者显式传入当前工具链、缓存来源和任务输出目录：

```bash
python3 scripts/peri_a0_reproduce/plan.py \
  --source "$TASK/peri" \
  --patch scripts/peri_a0_closeout/evidence/candidate-testfix.patch \
  --toolchain "$TOOLCHAIN" \
  --cache "$CACHE" \
  --task-root "$TASK" \
  --historical-manifest "$MANIFEST"
```

`$TASK` 的父目录必须已经存在。`$TASK` 本身应为空，或者只含作为直接子目录的 Peri 检出。`$TOOLCHAIN` 是只读 Rust/Cargo 1.99.0，其中有可执行的 `bin/rustc` 和 `bin/cargo`。`$CACHE` 位于任务目录之外，并含 `registry` 或 `git`。`$MANIFEST` 用来核验固定身份；其中的绝对路径不会写入新的命令或环境。

预检核验 Peri 提交、仓库内补丁摘要、`Cargo.lock` 摘要、`rust-toolchain.toml` 的 channel，以及 `rustc -V` / `cargo -V`。通过后只写 `$TASK/reproduce-plan.json`。计划中的 `PATH` 来自 `--toolchain`，`HOME`、`CARGO_HOME` 和 `CARGO_TARGET_DIR` 来自 `--task-root`。`CARGO_NET_OFFLINE` 固定为 `true`。Cargo 命令是相对路径的 `./scripts/cargo-rmcp-patched.sh build --locked ...`，工作目录是 `--source`。

预检会拒绝把任务目录放进原实验目录 `/tmp/hw-peri-a0-close-20261109`、`/tmp/helixweave-peri-a0-20261109`、`/tmp/helixweave-peri-design-20261108`，放进历史 manifest 记录的本地绝对路径、`--cache` 或 `--toolchain`，以及放进含 `conda-meta` 或 `pyvenv.cfg` 的共享环境。源码树落在原实验目录或那些历史路径里也会拒绝。与历史可写 target 重叠时同样拒绝。拒绝时不写计划文件。

这个预检不证明从零构建。原验收构建复制过已核验的 Cargo 缓存和 `target`，那些缓存不在本仓库。缓存缺失时不能改用 `--online-locked`，因为本入口不接受该参数；应停止并另行准备任务本地缓存。真正的构建由调用者在阅读计划后自行执行。native 运行限制仍是 512MiB、swap 0、1 核、128 任务。已有 build manifest 之后，冻结的 `run.py` 才使用 `--source`、`--build-manifest` 和 `--evidence`。ACP 编译若要重跑，仍是 8GiB、单编译任务；编译峰值不是运行需求。

`pyproject.toml` 的 Ruff `extend-exclude` 只列出已经冻结的历史 Python 文件。`scripts/peri_a0_reproduce/` 和新增测试不在其中。以后新增的维护代码不会因为放进三个历史目录名就自动豁免。

## 空白检查

冻结补丁里的上下文空白和原始日志末尾空行是实验字节，不裁剪。`.gitattributes` 只对出现这些问题的具体证据文件关闭 `blank-at-eol` 和 `blank-at-eof`，其他空白规则仍然有效。本页、复现入口和测试不在该列表中。

审阅整个归档范围时比较两个提交，不用提交后空工作区的 `git diff --check` 代替：

```bash
git diff --check b2f15868c6f74ab8764a73ae6fe4394f67304475..<archive-tip>
```

`<archive-tip>` 是包含上述属性的归档分支尖端。

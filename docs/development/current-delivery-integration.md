# 当前 PR 收口与集成清单

## 2026-11-16T08:22:10Z 当前状态

排期见 [roadmap](workflow-platform-agent-roadmap.md)，实施分工见
[session 计划](next-delivery-session-plan.md)。12:18:06Z 的快照及更早记录不倒写。
本页在 Draft #213，仍未合并。本次更新不写回自身新 SHA。

受限配置建议与校验 A1 已交付并进入 main
`1b1db421ba3172bc39cf31ce0d03976511379fda`。
这是 [#219](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/219)
的普通 merge，父提交为 `101983f9480bb72306da69875e79dc6ad65bca5a`
与批准 head `e112223fd03795b43cce6c8dce4566dcb4739f4b`。
merge 与批准 head 的文件树都是 `1f9ac53f6b34625fcbb6c09965427a3d6867e1ef`。
C4、launcher 生命周期和最终 Bulk Gate 不再写成当前阻塞。

证据身份分开记录：

- Bulk Gate 跑在 `e112223…`。[CI 43792685991](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43792685991)
  为 workflow_dispatch attempt 2 success，head 就是该批准提交。
  job `133061548271` success，记录为 4 项执行测试和 5 项产品浏览器测试通过。
  Artifact `14195706467` 的 ZIP SHA256 为
  `54ffa1c4aa24454ffc0dc63e8d6a0728004345c977c2633b0e05b5985adb88a6`。
  本次查询确认了 run、attempt、head 和该 job 的 success，没有重新下载 ZIP。
- C4 跑在 `82d091278a5cf46f5c33b40f4cbb26d3d12fe2a8`：Vite 页面、真实平台、
  固定 Peri、模型 API stub，5 passed。
- 文件树相同，因此这两项证据可以复用到 merge。不能写成 `1b1db4…` 上新跑过 Gate 或 C4。
- 合入后 push attempt 1 全部 success，head 为 merge SHA：
  [CI 43816767856](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43816767856)、
  [Lint 43816767900](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43816767900)、
  [Lock Check 43816767895](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43816767895)。
  两个 shard、frontend、browser-e2e、fast-checks、coverage 成功。
  real-execution、platform-real-execution、container-smoke、bulk-rnaseq-real-execution
  为 skipped，不是执行通过。

2026-11-16T08:22:10Z，GitHub 把 #212、#214、#218、#219 标为 merged。
它们的已审阅内容包含在上述 main 中。#216 的科学内容此前已在 `101983f…`；
同一次查询它也被标为 merged，记录的 merge 提交仍是 `101983f…`，不是一次新的合并。

当前进行中只有两件：

- #217 仍是 open Draft，head 仍为 `c4c49098d3c05b4d6913a684b459508c866c39b7`。
  它独立推进，与新 main 的组合仍需核验。12:18:06Z 对该 head 的取消恢复复审不改写。
- #213 仍是文档 Draft，未合并。

后续阶段都还没实施。ACP 兼容不是自动通过接入。已有 A1 草稿采用权限保持不变。

- 真实模型 API：设计和实施尚未开始，不要求本地模型。
- Agent runtime 适配：平台拥有授权、项目范围、草稿、确认、回执和生命周期；
  Peri 特有协议与版本留在适配层。这是设计要求。
- Agent 提交与取消运行：独立操作契约、确认和幂等仍待设计。
- 数据目录：按课题、细胞类型、实验类型分类与查询仍是设计。
- AI 整理：成员查询和建议，管理员在网页或终端确认。不把已有 A1 草稿采用改成管理员专属。
- managed execution、生产容量和科学确认另列，不属于这次 A1 交付。

非阻塞注记：relay 第三连接测试 `test_relay_refuses_a_third_connection_and_closes`
仍只证明某个 handler 启动。断言补强不改生产 relay。若另交测试提交，
`82d0912…` 的 C4 证据保留原 SHA。

## 2026-11-15T12:18:06Z 快照

排期见 [roadmap](workflow-platform-agent-roadmap.md)，实施分工见
[session 计划](next-delivery-session-plan.md)。07:21:49Z 的 `23681d9…` 快照及更早记录不倒写。
本页在 Draft #213，仍未合并。本次更新不写回自身新 SHA。
同一次查询：main 仍是 `101983f…`；#219 head 为 `6d570822620ccb873c26a1dabad1a8baebe5593b`；
#217 为 `c4c49098d3c05b4d6913a684b459508c866c39b7`；#214、#218、#212、#216 的 head 未变。

- [#219](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/219)
  已是 open Draft。已审阅 head 是 `6d5708…`。短命创建线程已由常驻 launcher 修复。
  新的待返修：创建超时后，晚到进程可能无人接管。当前不再把 `d2a0cc…` 加未提交接线当作现状。
- C4 本地记录有 5 passed，完整验收未完成。没有采用后再校验；丢响应测试在请求到达 API 前中止；
  未执行 Reject；“过期保存”实际是另一次写入造成的 revision 冲突。
  原先 2 passed、1 failed、2 skipped 保留为更早历史证据。这次本地 JUnit 早于最终干净 SHA，
  不是该 SHA 或生产静态包上的运行。
- #219 的 pull_request [CI 43693526422](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43693526422)
  为 attempt 1 success。两个 shard、frontend、browser-e2e、fast-checks、coverage 成功；
  四个条件执行 job skipped。
  同一 head 的 workflow_dispatch [43693536449](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43693536449)
  为 attempt 1 failure。real-execution、platform-real-execution、container-smoke、frontend、
  browser-e2e、两个 shard 和 fast-checks 成功。
  Bulk job 步骤数为 0，在 runner 执行前被分支保护拒绝，不是科学执行失败。
  coverage job `132755072368` 失败：services 全套覆盖率 85.37%，低于 87.28%。
  该步骤里更后面的门槛没有执行，不能记为通过。不能把这次 dispatch 写成只剩环境准入。
  本次查询核对了 run 与 job 结论，没有重新下载 coverage 日志，也没有重跑这两次 CI。
- #217 `c4c4909…` 的取消恢复增量已通过独立定点复审。`db3d8a5…` 的 P2 不再列为未修复。
  原回归在 `db3d8a5…` 失败、在 `c4c4909…` 通过。准确提交的生产浏览器证据仍待收口；
  作者浏览器记录仍标在 `db3d8a5…`，这次复审没有重跑浏览器证明。#217 未纳入 #219。
- #214 `a566943…` 与 #218 `ca40bce…` 已各作为祖先纳入 #219 一次。Frontend 待命。
  Science 没有新增任务。169 passed 仍属于 `0fef7a2…`。
- #212 `ff4bdd5…` 仍是已复审的文案与署名增量。C4、整个 A1、最终组合 Gate、
  真实 provider 和 Agent 运行操作分别记账。
- `23681d9…` 的 shard 失败里，PID 123 的 SIGKILL 是测试 mock 捕获的 `os.kill` 调用，
  不是已证明向宿主 PID 123 实际发信号。该 CI 结论仍是 failure，不改写成通过。

## 2026-11-15T07:21:49Z 快照

排期见 [roadmap](workflow-platform-agent-roadmap.md)，实施分工见
[session 计划](next-delivery-session-plan.md)。19:18:06Z 及更早快照不改写。
本页在 Draft #213，仍未合并。本次更新不写回自身新 SHA。
同一次查询：main 仍是 `101983f…`；#212 head 为 `ff4bdd5…`，#214 为 `a566943…`，
#217 为 `db3d8a5…`，#218 为 `ca40bce…`，#216 仍 open。没有为这些文档增量另查一套新 CI。

- C4 已尝试，未完成。手工保存、刷新回读和真实 snapshot 成功。
  native 创建线程生命周期缺陷阻断提案及后续步骤。原 JUnit 为 2 passed、1 failed、2 skipped。
  集成 HEAD 是 `d2a0cc683aefd3a05fe507dce0fba60e79ecc320`，另有未提交接线；交接时没有集成 Draft。
  受控复现确认 SIGKILL 来自短命 handler 线程退出后的 `--die-with-parent`，不是 OOM 或资源不足。
  两个确认性测试通过只说明缺陷已复现，不表示 C4 通过。
- #217 `db3d8a51af10c51ba77c5100b4f20e6a58664bf6` 的生产外发 P1 已通过复审。
  新 P2：取消进行中的预检后再选择，主动 AbortError 被留在失败集合，成功重选不能恢复。
  等待窄修复。该复审没有跑完整 CI、Gate 或 native 测试，不是 Viewer 全部通过。
- #212 `ff4bdd5dbb206d653e73fed6627bfc0979cea568` 的文案与署名增量已通过复审，生产和测试代码未变。
  `a08aa62…` 的摘要绑定结论仍只属于该 SHA。
- #218 `ca40bced3e0e954b87e7cde82017ca8dfd2d6d53` 的单行文档增量已通过复审，生产和测试代码未变。
  169 passed 的真实执行证据仍属于 `0fef7a2…`，不改写该测试身份，也不重跑科学层。
- #214 `a566943…` 继续可消费，没有新的 C3 任务。
- C4、整个 A1、最终组合 Gate、真实 provider 和 Agent 运行操作分别记账。
  不因某个 session goal complete 写成全部完成。原 binary `60144bbe…` 的缺失和历史 A0 记录保留。

## 2026-11-14T19:18:06Z 快照

排期见 [roadmap](workflow-platform-agent-roadmap.md)，实施分工见
[session 计划](next-delivery-session-plan.md)。下面是各 owner 本轮最终交接。
17:40:08Z 及更早的快照保留在后文，不改写成这次才知道的结论。
本页已在 Draft #213，该 PR 仍未合并。本次更新不写回自身新 SHA。

| 项 | 准确 head | 本轮结论 | 下一步 |
| --- | --- | --- | --- |
| [#214](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/214) | `a56694394d6c02d8375ec921b0f3870514c74f71` | 保存返修复审通过。普通 [CI 43622972873](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43622972873)、[Lint 43622972896](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43622972896)、[Lock Check 43622972843](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43622972843) 均为 pull_request attempt 1 success。四个条件执行 job skipped | 进入 C4 准备。不是整个 A1 批准 |
| [#218](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/218) | `0fef7a2100755f69778771f1a37e6f9632aa9045` | 真实执行 169 passed、0 failures/errors/skips。整次 [CI 43628729299](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43628729299) 为 workflow_dispatch attempt 1 success | 保留该 head。sample-sheet 的 `final.bam` 文案是非阻塞项，不因此重跑科学层 |
| [#212](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/212) | `a08aa62a3df6a8253b96bd63b1abf2e9b1ffa22e` | 摘要绑定的窄增量复审通过。[CI 43629047387](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43629047387)、[Lint 43629047391](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43629047391)、[Lock Check 43629047437](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43629047437) 均为 pull_request attempt 1 success | 供 C4 使用该 head。不是整个 A1，也不是最终 Gate |
| [#217](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/217) | `aa0ef9fe86caf3aee99e35bae1687afdbf48f2d8` | gzip 与空白返修通过。新 P1：生产读取失败会尝试外部 CORS 代理 | 独立返修该回退。不是 Viewer 全部通过 |
| #216 | `85ac7a2…` 已在 main `101983f…` | PR 对象仍 open。禁止再次合并 | 不把 PR 异常当成未合入代码 |

#218 的 workflow success 来自整次 run，不是从 real-execution 一个 job 推断。
该 run 中 real-execution、两个 shard、frontend、browser-e2e、platform-real-execution、
container-smoke、fast-checks、coverage 为 success；bulk-rnaseq-real-execution 为 skipped。
复审记录 real-execution job `132550360035` 的 checkout 就是 `0fef7a2…`。
本次查询核对了 run 与 job 结论，没有重新下载 checkout 日志。skipped 不是执行通过，也没有派发 Gate。

#212 同一 CI 的四个条件执行 job 为 skipped。窄增量没有重跑 23 项 native 矩阵或 ACP 大构建。
119 个 Bulk 文件相对 `e00130e…` 未变，aggregate 仍是
`979bf1f410471d2ba3d78fcb473c17c688ba24afe0ecaf1af96ed08528e4d2ae`。
这次摘要绑定不增加 Gate；整个 A1 最终组合仍要适用的 Gate。
原 binary `60144bbe…` 的缺失和历史 A0 记录保留。`5cc82caf…` 只在固定受限范围内被接受为替代运行身份。

#214 复审记录其浏览器 job checkout 为 `e86f1ffea90ca7a7df8a2d208def5c7c42ce00d8`，
PR head 为 `a566943…`。本次没有重下该日志。schema smoke 与浏览器通过数以那次复审为准。

C4 由 Backend 在独立 worktree 主导，Frontend 配合，组合 #212、#214 和 #218 的上述 head。
Viewer 返修独立推进。#217 未修好时不阻塞 A1 准备，但不能作为已验收 Viewer 纳入交付。
C4、整个 A1、真实 provider、Agent 运行操作和最终组合 Gate 都还没有执行。

非阻塞文案：#218 `docs/sample-sheet.md` 的对照示例仍写 `final.bam`，应为 `analysis.bam`。
#212 资格说明把历史 ACP 单测写成旧 binary 上的执行，应改为独立 peri-acp harness 的源码层证据。
两处都不因此重跑科学层或重新编译。

## 2026-11-14T17:40:08Z 快照

排期见 [roadmap](workflow-platform-agent-roadmap.md)，实施分工见
[session 计划](next-delivery-session-plan.md)。本页保存身份及证据边界，不替代源码审阅。
下面先给这一时刻的下一步。15:57:43Z 的检查记录留在后文，不改写成新的查询结果。

- main 仍是 `101983f…`。#216 代码已进入 main，PR 对象仍 open、未标记 merged。禁止再次合并。
  dispatch CI `43612614528` 已 failure。源码核对见测试 PATH 遗漏 Picard，
  且测试准备仍保留旧 profile 分流和复制 BAM。下一步由 Science 用独立小 PR 修复，
  并验证完整真实执行层。不以 `--maxfail=1` 已经跑过的数量当作整个套件通过，
  也不写没有剩余阻塞。`85ac7a2…` 的 schema 审阅不因此被取消。
- #214 `399c3b3…` 的独立复审要求返修。PE→SE 保存确认后草稿仍 dirty；
  样本复选框与实际保存成员不一致。原定向测试 100 passed，独立新增的 3 项回归断言失败。
  这不是 C4 或整个 A1 验收。下一步是返修这两项后再复审。
- #217 `cbe5a0e…` 的独立审阅已完成：多成员 gzip 在第一个 member 后漏掉后续数据，
  PR 范围空白检查不干净。`414872e…` 的父提交是该审阅 head；17:40:08Z 的 head 是
  `aa0ef9fe86caf3aee99e35bae1687afdbf48f2d8`，其父提交为 `414872e…`。
  新 head 等正式交付和复审，不转移旧裁决，也不提前宣布修复通过。
- Peri 候选 `5cc82caf…` 仍未通过资格。宿主支持 systemd/cgroup v2，用户服务未启动。
  启动 `user@1000.service` 并验证任务资源隔离的限定指令已经交给 Backend；
  没有实际执行证据之前，不更新为环境合格或候选通过。
- #212 head 仍是 `e00130e…`，本轮没有新审阅。#213 仍是 open Draft，本更新不写回自身新 SHA。

## 2026-11-14T15:57:43Z 查询快照

远端 PR/main 的上一份查询时间是 2026-11-14T15:57:43Z。head 未移动的检查和审阅继续只适用于原 SHA；
新 head 的 CI、源码审阅和产品验收分开记录，不继承旧 head 的批准。
新提交出现后，本表只作为本次快照，不自动沿用 CI 或审阅结论。
本页已在 Draft #213 中，该 PR 仍未合并。本次更新不把随后生成的提交 SHA 写回本文件。

## 已合入基线

2026-11-14T15:57:43Z 的 main 为 `101983f9480bb72306da69875e79dc6ad65bca5a`。
这是 GitHub 生成的 merge 提交，父提交为 `5bc21a040951ba28333bcc04dfb151735089fcd6`
与 `85ac7a26b9958ead4a9ddde6e56dc38c07adeca2`，说明为 Merge pull request #216。
该 SHA 没有 push 检查。workflow_dispatch、attempt 1 的
[Lint 43612607625](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43612607625)
与 [Lock Check 43612610860](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43612610860)
在 2026-11-14T15:57:43Z 已是 success。
[CI 43612614528](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43612614528)
当时 in progress，于 2026-11-14T16:12:17Z 结案为 failure。
real-execution 失败，步骤名为 Scientific real-execution suite 与
Require passing, zero-skip scientific outcomes。frontend、browser-e2e、
platform-real-execution、两个 shard、container-smoke、fast-checks、coverage 为 success；
bulk-rnaseq-real-execution 为 skipped。本次没有下载失败日志，
不把该失败写成对 `85ac7a2…` 审阅的否定，也不是新的 Gate。
commit 在 15:57:43Z 没有关联 PR。
[#216](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/216)
的 Pulls API 同时仍是 open、merged=false、非 Draft。
不把这个 PR 对象写成已关闭，也不把已审阅 head 的 pull_request 检查改记成 merge SHA 的检查。

`5bc21a0…` 是
[#215](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/215)
的 merge，修复 browser-e2e 的 apt 镜像步骤和契约测试。
该 SHA 的 push 检查均为 attempt 1、success，不是 `101983f…` 的检查：

- [CI 43210921001](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43210921001)。
- [Lint 43210921261](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43210921261)。
- [Lock Check 43210921108](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43210921108)。

上述 job 结论属于已记录的 `5bc21a0…` push 检查。本次没有重新下载这些 job 日志或 artifact。
两个 shard、frontend、browser-e2e、fast-checks、coverage 成功；四项条件真实执行 job skipped。
另有同 SHA 的 schedule CI 43336374574 success，不混作上述 push attempt 的证据。
#210 merge `b2f15868c6f74ab8764a73ae6fe4394f67304475`、
#211 merge `ff30329607c4b6c9849bfc7369906e81b428937d` 已在此 main 历史中。

## 查询时的 PR

下表是查询时的 PR head 与分支基线，不等于 PR CI 的合成 checkout。
#212、#213、#214、#217 为 open Draft。#216 为 open、非 Draft；main 已有包含其 head 的 merge 提交，
PR 对象当时未标记 merged。#194/#195 为 open、非 Draft。本次不合并或关闭任何 PR。

| PR | head | 记录基线 | 当前动作 |
| --- | --- | --- | --- |
| [#212](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/212) | `e00130e2d5e55b3197624ef6bff9c25fd687cf87` | `b2f15868c6f74ab8764a73ae6fe4394f67304475` | head 未变。2026-11-13 边界复审仍只适用于该 SHA，可交 C3；本轮无新审阅 |
| [#214](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/214) | `399c3b34b42810f6a042a8aa3d86f92ae005eb6f` | `ff30329607c4b6c9849bfc7369906e81b428937d` | 17:40:08Z head 未变。独立复审要求返修；C4/A1 未完成 |
| [#216](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/216) | `85ac7a26b9958ead4a9ddde6e56dc38c07adeca2` | `ff30329607c4b6c9849bfc7369906e81b428937d` | 代码已在 main。PR 仍 open，禁止重复合并。真实执行测试仍有剩余阻塞 |
| [#217](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/217) | `aa0ef9fe86caf3aee99e35bae1687afdbf48f2d8` | `ff30329607c4b6c9849bfc7369906e81b428937d` | `cbe5a0e…` 审阅要求返修。当前 head 比 `414872e…` 再前进一提交，等正式复审 |
| [#213](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/213) | `ec7779b4a0771e98f8253f3b1ed8d72c09740f05` | `5bc21a040951ba28333bcc04dfb151735089fcd6` | 当前入口已在该 Draft；查询时未合并。本状态更新不写回自身新 SHA |
| [#194](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/194) | `4825ed62bf94b86bdb9afe5b3896e97ea46c2cea` | `fa0d8401cbc3fecb765fed6d7e01e18bb8fbc96d` | head 未变。open、非 Draft；旧 Cloud 工具链说明待内容复核 |
| [#195](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/195) | `a64434962c3f497a83188b46f8b2f488a353cd8e` | `fa0d8401cbc3fecb765fed6d7e01e18bb8fbc96d` | head 未变。open、非 Draft；审计建议不直接转成新开发批次 |

记录基线不承诺新 session 仍从旧 main 开始；实施前核实目标分支、merge-base 和检查 checkout。

## 已记录 head 的普通检查

以下均为 attempt 1、pull_request，并只覆盖表中写出的 SHA。Lint/Lock 的准确 run 和 head 也已核对。

| PR | CI | Lint / Lock Check | 结论 |
| --- | --- | --- | --- |
| #212 | [43383020734](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43383020734) | [43383020771](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43383020771) / [43383020680](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43383020680) | head 未变，三项 success；A1/C4 未验收 |
| #214 | [43610191318](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43610191318) | [43610191450](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43610191450) / [43610191314](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43610191314) | `399c3b3…` 三项 success。普通 CI 不取消返修结论 |
| #216 | [43409346164](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43409346164) | [43409346644](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43409346644) / [43409346283](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43409346283) | 只覆盖 `85ac7a2…` 的 pull_request，三项 success |
| #217 `cbe5a0e…` | [43409373596](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43409373596) | [43409373630](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43409373630) / [43409373577](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43409373577) | 只覆盖该旧 head 的普通 CI。不覆盖 `aa0ef9f…`，本轮未查询新 head 的 CI |
| #213 | [43409481140](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43409481140) | [43409481214](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43409481214) / [43409481220](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43409481220) | 只覆盖查询时 head `ec7779b…`，不覆盖本次状态更新 |

#214 旧 head `2c3c439…` 的 [CI 43222397611](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43222397611)
为 cancelled，只属于该 SHA。新 head `399c3b3…` 的 CI 已 success：两个 shard、frontend、
browser-e2e、fast-checks、coverage success；四个条件执行 job skipped。
这不是该 head 的返修验收或 C4/A1 验收。
#212 与 `85ac7a2…`、`cbe5a0e…` 的上表 pull_request 普通 jobs 为 success。
#217 当前 head `aa0ef9f…` 的 CI 本轮没有查询。#216 merge SHA 的 workflow_dispatch CI 已失败。
这些已完成 PR CI 的 real-execution、platform-real-execution、container-smoke、
bulk-rnaseq-real-execution 为 skipped，不计实际执行通过。
#216 的 PR 合成提交元数据为 `5865e552c9fb03858f4dae5868c5bddf88c876da`，
父提交为 `5bc21a0…` 与 `85ac7a2…`。#217 的对应元数据为
`5a903d12cd78d5b332c39927ad1a837e41824f6f`，父提交为 `5bc21a0…` 与 `cbe5a0e…`。
本次没有重新下载这些 job 的 checkout 日志。PR 记录基线仍是 `ff303296…`。

## 已有审阅的适用范围

#212 的下列证据属于未移动的 `e00130e…`，不是 2026-11-14 的新审阅。
该 head 重跑过原独立边界复现、文件访问和目录到真实 snapshot：
41 passed。祖先目录移走再恢复仍拒绝生成 snapshot；归档池返回空目录；
合法输入正例仍通过。F1/F2 复审通过，C3 可消费此版本的目录及样本计划契约。
`control_sample_revision_id` 位于 `workflow_samples` 的计划行内，不是草稿顶层字段。
目录不发给模型。监视不覆盖 SQLite 插入，后续执行仍需独立资格核对。

#216 `85ac7a2…` 的 2026-11-14 整体独立审阅当时没有发现继续阻止合入的执行层缺陷。
该句只描述那次审阅，不是当前“没有剩余阻塞”。审阅核对了 base 到 head 的范围，并复用 `c44c4d6…` 的 Picard 定点证据，
没有把那次定点复审重跑，也没有把它写成唯一结论。该独立审阅运行 87 项契约和 7 项 DAG，
是合成 fixture 与 dry-run，不是大规模科学矩阵、真实 Picard 重跑或生物学确认。
两项非阻塞文档整理建议不改变“可作为 A1 已审阅科学 schema 目标”的结论。
main 上的 `101983f…` 包含这个 head。17:40:08Z 的 PR 对象仍 open、未标记 merged，禁止再次合并。
merge SHA 没有 push 检查；其 workflow_dispatch CI 已于 16:12:17Z failure。
随后的源码核对把失败定位到测试 PATH 遗漏 Picard，以及仍保留旧 profile 分流和复制 BAM 的测试准备。
这是剩余阻塞。Science 用独立小 PR 修复并验证完整真实执行层。
`--maxfail=1` 下已经运行的数量不是整个套件通过。这也不等于 C4、A1 或科学默认冻结。

`c44c4d6…` 上的 Picard 定向回归仍是 27 passed（24 项配置，3 项真实规则测试）：
直接删除开关被拒绝，默认路径只标记且保留输入 records，缺 Picard 不回退。
arguments_file 的四种既有实测在解析阶段失败，不写成已证实的删模板漏洞，
也不外推到未来工具版本。正式支持时须另定契约。

#217 `cbe5a0e…` 的独立审阅要求返修。多成员 gzip 在解压到第一个 member 结束时返回，
没有处理或拒绝后续 member，合法双 member BEDPE 的派生结果少了一行。
PR 范围空白检查退出码为 2，问题在生成的前端资源尾随空白，不手改 bundle，
也不因此全局关闭空白检查。`bcdbd32…` 的静态抽查留在更早的 head。
`414872e…` 和当前 `aa0ef9f…` 都还没有这次审阅之后的正式复审。
上述状态不替代 C4、最终组合 CI 或 Protected Bulk Gate。

## Peri runtime 候选

历史 A0 接受的运行 binary 是 `60144bbe0919992b8913b9031696b2e30b1e9f82d95048ae3b07c887d5d85517`。
该产物现已缺失。固定 Peri `c7d154c…`、固定补丁和当时验收结论保留在
[A0 归档](peri-a0-archive.md)，不改写成历史失败。

新构建候选 SHA256 为
`5cc82caf6e04578f19ff5790fbc4eb8c9e115793ff6950a699c9e2ad28c50f47`。
2026-11-13 的只读复核确认产物字节、权限和摘要，以及源码、补丁和锁文件身份；
它只允许该候选进入隔离的资格验证，没有运行候选或重新编译。
同一源码和补丁不能单独证明新二进制行为，也不能恢复 `60144bbe…` 的身份或继承其 native 验收。
宿主支持 systemd 和 cgroup v2，但用户服务未启动，所以还不能把宿主支持写成用户会话合格。
Backend 已获得启动 `user@1000.service` 并验证任务资源隔离的限定操作指令。
本轮没有执行该指令，也没有操作服务。只有收到实际执行证据后，才能更新为环境合格或候选通过。
完成前，该候选不是 A0、A1 或 C4 runtime。

## Gate 与科学边界

- 已有 Gate [42553090512](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/42553090512)
  对应 `6978349ca515a37a55a93abc3ec65519fd0be2f8`。2026-11-13 查询时 jobs 均 success；
  本次没有重新下载或核算 artifact。#210 merge 的同树证据见历史交接；这不是后续 SHA 的新 Gate。
- #212 正文登记 aggregate 为 `979bf1f410471d2ba3d78fcb473c17c688ba24afe0ecaf1af96ed08528e4d2ae`，
  与旧闭包不同，最终候选需独立 Gate。该值来自交接，本轮未重新生成。
- #216 `85ac7a2…` 的整体审阅核对 119 个 Bulk 实现文件与已提交 manifest 一致，
  aggregate 仍为 `94b577b62f31ed53264ccb77b9f4d92d61246bca280398e79d36ff47eb481f36`。
  该核对不是新 Gate，也不覆盖最终 A1 组合。#217 新 head 本轮没有新的 Bulk 结论。
  最终组合仍检查文件集合、字节、runtime/参考身份和接线。
- SCI-7 本批 no conclusion；JUND 跨靶点检查完成。其他未执行 confirmation、
  MNase 未评估指标及默认未冻结的边界见 roadmap，不由 CI 代替科学裁决。

## 集成安排

1. 原 session 交出准确 head、范围、生成物、测试和未完成项；Codex 独立复审。
2. A1 已随 #219 进入 main。#216 不需要再次合并。169 项真实执行证据仍是 `0fef7a2…`。
   创建超时和当时的覆盖率失败留在 12:18:06Z 快照，不再是当前阻塞。
   如以后公共契约变化，后端正式生成 OpenAPI/client，再交前端适配，不手改生成物。
3. #214 已含部分 #212 历史，合成候选前核对祖先关系。只纳入已审阅提交，
   不将整个活动分支盲目拉进候选。候选使用独立 worktree/PR。
4. #217 `c4c4909…` 继续独立推进。它没有进入 #219。与 main `1b1db4…` 的组合仍需核验，
   因此不作为已经随 A1 交付的 Viewer。
   它与 #214 的公共导航、前端打包资产由集成 owner 协调重生成。
   最终测试按实际组合执行，不能简单相加来源分支的通过记录。
5. 普通 CI、C4 和适用 Gate 完成后交回主会话；合并/发布仍单独执行授权。
   本次只更新 Draft #213 的当前状态，不创建、关闭或合并远端 PR，不启动其他 session。

## 历史证据

[整理前的完整集成记录](current-delivery-integration-history-through-20261112.md)逐字保留。
`f83bd7d…` 与查询时的 `ec7779b…` 都留在 #213 历史中。当前入口以本文件和 roadmap 为准；
本次状态更新提交后 head 会前进，该新 SHA 不写入本文件。

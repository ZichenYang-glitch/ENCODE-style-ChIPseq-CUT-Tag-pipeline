# 当前 PR 收口与集成清单

## 2026-11-14T17:40:08Z 当前状态

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
2. #216 不再合并。A1 继续使用 `85ac7a2…` 这个已审阅 schema，同时等待真实执行测试小 PR。
   A1 C4 使用该 schema，
   并等待完成资格验证的 runtime。如公共契约变化，后端正式生成 OpenAPI/client，
   再交前端适配，不手改生成物。
3. #214 已含部分 #212 历史，合成候选前核对祖先关系。只纳入已审阅提交，
   不将整个活动分支盲目拉进候选。候选使用独立 worktree/PR。
4. #217 先完成新 head 的正式复审，再谈独立验收；它与 #214 的公共导航、前端打包资产由集成 owner 协调重生成。
   最终测试按实际组合执行，不能简单相加来源分支的通过记录。
5. 普通 CI、C4 和适用 Gate 完成后交回主会话；合并/发布仍单独执行授权。
   本次只更新 Draft #213 的当前状态，不创建、关闭或合并远端 PR，不启动其他 session。

## 历史证据

[整理前的完整集成记录](current-delivery-integration-history-through-20261112.md)逐字保留。
`f83bd7d…` 与查询时的 `ec7779b…` 都留在 #213 历史中。当前入口以本文件和 roadmap 为准；
本次状态更新提交后 head 会前进，该新 SHA 不写入本文件。

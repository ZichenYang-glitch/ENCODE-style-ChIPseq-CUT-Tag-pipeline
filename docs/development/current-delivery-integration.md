# 当前 PR 收口与集成清单

## 2026-11-13 查询快照

排期见 [roadmap](workflow-platform-agent-roadmap.md)，实施分工见
[session 计划](next-delivery-session-plan.md)。本页保存身份及证据边界，不替代源码审阅。
远端 PR/main 查询时间 2026-11-13 09:31 UTC；随后按准确 head 查询 Actions。
2026-11-13T10:26:21Z 再次查询，main 与 #212、#214、#216、#217、#213 的 head
没有变化，因此下表审阅只适用于这些 SHA。
新提交出现后，本表只作为本次快照，不自动沿用 CI 或审阅结论。

## 已合入基线

远端 main 为 `5bc21a040951ba28333bcc04dfb151735089fcd6`，即
[#215](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/215)
的 merge。该 PR 修复 browser-e2e 的 apt 镜像步骤和契约测试。
main 的 push 检查均为 attempt 1、success：

- [CI 43210921001](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43210921001)。
- [Lint 43210921261](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43210921261)。
- [Lock Check 43210921108](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43210921108)。

本轮核对 run 的 head/event/attempt 与 job 结论，未重新下载全部 checkout 日志或 artifact。
两个 shard、frontend、browser-e2e、fast-checks、coverage 成功；四项条件真实执行 job skipped。
另有同 SHA 的 schedule CI 43336374574 success，不混作上述 push attempt 的证据。
#210 merge `b2f15868c6f74ab8764a73ae6fe4394f67304475`、
#211 merge `ff30329607c4b6c9849bfc7369906e81b428937d` 已在此 main 历史中。

## 未合入 PR

下表是 PR head 与分支基线记录，不等于 PR CI 的合成 checkout。
除 #194/#195 外均为 open Draft。本轮没有修改 PR、分支或合并状态。

| PR | head | 记录基线 | 当前动作 |
| --- | --- | --- | --- |
| [#212](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/212) | `e00130e2d5e55b3197624ef6bff9c25fd687cf87` | `b2f15868c6f74ab8764a73ae6fe4394f67304475` | 最新 F1/F2 定点复审通过，目录与样本计划契约可交 C3 |
| [#214](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/214) | `2c3c439f0c1bf3905caf4888486e16df0cd86134` | `ff30329607c4b6c9849bfc7369906e81b428937d` | 正在把目录接到真实 snapshot；该 head CI cancelled，C4/A1 未完成 |
| [#216](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/216) | `c44c4d6f43912b8dc26317c33a35d3a77ba885bd` | `ff30329607c4b6c9849bfc7369906e81b428937d` | Picard 返修定点复审通过，整体合入审阅待完成；尚未成为 main 行为 |
| [#217](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/217) | `bcdbd3256afb2eb61be955003fc318985cfca9d1` | `ff30329607c4b6c9849bfc7369906e81b428937d` | 最新 head 普通 CI 通过，待独立审阅 |
| [#213](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/213) | `f83bd7dff34f4bbfb1bda85879634d11491f4c72` | `ff30329607c4b6c9849bfc7369906e81b428937d` | 与本轮文档协调一次纳入，旧快照不作当前状态 |
| [#194](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/194) | `4825ed62bf94b86bdb9afe5b3896e97ea46c2cea` | `fa0d8401cbc3fecb765fed6d7e01e18bb8fbc96d` | open、非 Draft；旧 Cloud 工具链说明待内容复核 |
| [#195](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/195) | `a64434962c3f497a83188b46f8b2f488a353cd8e` | `fa0d8401cbc3fecb765fed6d7e01e18bb8fbc96d` | open、非 Draft；审计建议不直接转成新开发批次 |

记录基线不承诺新 session 仍从旧 main 开始；实施前核实目标分支、merge-base 和检查 checkout。

## 当前 head 的普通检查

以下均为 attempt 1、pull_request。Lint/Lock 的准确 run 和 head 也已核对。

| PR | CI | Lint / Lock Check | 结论 |
| --- | --- | --- | --- |
| #212 | [43383020734](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43383020734) | [43383020771](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43383020771) / [43383020680](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43383020680) | 三项 success；A1/C4 未验收 |
| #214 | [43222397611](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43222397611) | [43222397634](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43222397634) / [43222397628](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43222397628) | CI cancelled；Lint/Lock success |
| #216 | [43385234431](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43385234431) | [43385234478](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43385234478) / [43385234417](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43385234417) | 三项 success；Picard 返修定点复审通过 |
| #217 | [43394445393](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43394445393) | [43394445120](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43394445120) / [43394445154](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/43394445154) | 三项 success；本轮未重新执行绘图验证 |
| #213 | [42775173026](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/42775173026) | [42775173062](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/42775173062) / [42775173268](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/42775173268) | 三项 success；内容仍是较早快照 |

#214 的 browser-e2e、frontend、shard 1 success；shard 2 cancelled，fast-checks 和 coverage failure。
取消原因本轮未取证，不能称为已确认的产品缺陷，也不能称必过检查全部通过。
#212/#216/#217 的普通 jobs 均 success。
这些 PR CI 的 real-execution、platform-real-execution、container-smoke、
bulk-rnaseq-real-execution 为 skipped，不计实际执行通过。

## 同日定点复审：C3 可以开始接线

#212 在上述准确 head 上重跑原独立边界复现、文件访问和目录到真实 snapshot：
41 passed。祖先目录移走再恢复仍拒绝生成 snapshot；归档池返回空目录；
合法输入正例仍通过。F1/F2 复审通过，C3 可消费此版本的目录及样本计划契约。
`control_sample_revision_id` 位于 `workflow_samples` 的计划行内，不是草稿顶层字段。
目录不发给模型。监视不覆盖 SQLite 插入，后续执行仍需独立资格核对。

#216 的 Picard 定向回归 27 passed（24 项配置，3 项真实规则测试）。
直接删除开关被拒绝，默认路径只标记且保留输入 records，缺 Picard 不回退。
额外 arguments_file 调查未复现成功删除模板：固定工具与生产参数语法组合均解析失败。
这不是任意参数文件或未来工具版本的安全结论；正式支持时须另定契约。
本轮不重跑整个单路径科学矩阵，不将定点修复通过写成生物学确认。

#217 本轮仅静态抽查本地证明脚本与许可打包修复，未独立重跑真实绘图，仍待验收。
上述复审不改变 PR Draft 状态，也不替代 C4、最终组合 CI 或 Protected Bulk Gate。

## Gate 与科学边界

- 已有 Gate [42553090512](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/42553090512)
  对应 `6978349ca515a37a55a93abc3ec65519fd0be2f8`。本轮再次查询 jobs 均 success，
  未重新下载或核算 artifact。#210 merge 的同树证据见历史交接；这不是后续 SHA 的新 Gate。
- #212 正文登记 aggregate 为 `979bf1f410471d2ba3d78fcb473c17c688ba24afe0ecaf1af96ed08528e4d2ae`，
  与旧闭包不同，最终候选需独立 Gate。该值来自交接，本轮未重新生成。
- #216/#217 正文登记其 Bulk 闭包未变。最终组合仍检查文件集合、字节、runtime/参考身份
  和接线，不要求逐个文档/UI 提交机械重跑 Gate。
- SCI-7 本批 no conclusion；JUND 跨靶点检查完成。其他未执行 confirmation、
  MNase 未评估指标及默认未冻结的边界见 roadmap，不由 CI 代替科学裁决。

## 集成安排

1. 原 session 交出准确 head、范围、生成物、测试和未完成项；Codex 独立复审。
2. #216 可先独立进入合入审阅。A1 C4 面对最终科学 schema；如契约变化，
   后端正式生成 OpenAPI/client，再交前端适配，不手改生成物。
3. #214 已含部分 #212 历史，合成候选前核对祖先关系。只纳入已审阅提交，
   不将整个活动分支盲目拉进候选。候选使用独立 worktree/PR。
4. #217 可独立验收；它与 #214 的公共导航、前端打包资产由集成 owner 协调重生成。
   最终测试按实际组合执行，不能简单相加来源分支的通过记录。
5. 普通 CI、C4 和适用 Gate 完成后交回主会话；合并/发布仍单独执行授权。
   本轮交付文档与定点复审，不创建/关闭/合并远端 PR，不启动其他 session。

## 历史证据

[整理前的完整集成记录](current-delivery-integration-history-through-20261112.md)逐字保留。
上表 #213 head `f83bd7d…` 是协调前的查询值，保留在该 PR 历史中；
本提交纳入后，当前入口以本文件和 roadmap 为准。

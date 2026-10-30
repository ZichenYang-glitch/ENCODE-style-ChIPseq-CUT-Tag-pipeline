# 当前 PR 收口与集成清单

核对日期：2026-10-30。当前远端 main 为
`fa0d8401cbc3fecb765fed6d7e01e18bb8fbc96d`。本清单记录准确 head 的
远端检查、依赖和拟议交付顺序；本次未合并、重定基线、推送既有 PR，
也未执行最终组合版本。状态发生变化后必须重新查询，不能沿用本表批准新 head。

## 当前状态

17 个 PR 均 open、未合并；#194、#195 为非 Draft，其余为 Draft。
其中 15 个准确 head 的 CI、Lint、Lock Check 成功，#194、#195 的 CI 失败。
本次读取了各 CI 的 job 列表，成功状态不等同于已完成代码审阅、获准合并或科学验收。

| PR | 完整 head | 当前 base | CI | Lint / Lock Check |
| --- | --- | --- | --- | --- |
| [#184](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/184) | `483a56ed76bda57c7ce50bd61772c33d333e3639` | main | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/39736844720) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/39736844652) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/39736844896) |
| [#185](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/185) | `5c86076c8f2db1facbe0a6244c1360c8d4b09c69` | #184 | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/39790305946) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/39790305954) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/39790305945) |
| [#186](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/186) | `faa3d1dc714d4bf6105dcb1b004419a0c6c09016` | #185 | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/39915397365) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/39915397329) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/39915397346) |
| [#187](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/187) | `659e849a07f3129065fced1f44e2e2c95972f704` | #186 | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/39925518724) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/39925518747) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/39925518788) |
| [#188](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/188) | `0d7f7fd5c623ae559ff2490356a36fbceba5c29c` | #187 | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40180744887) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40180744833) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40180744834) |
| [#189](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/189) | `66d0ce1004f74fa43321bf5938e5ca61c6ff5eb9` | #188 | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40315543051) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40315542998) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40315543096) |
| [#190](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/190) | `b6f3a49f45a6b8093f8c42db6684bce5f1425475` | main | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40492125802) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40492125706) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40492125866) |
| [#191](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/191) | `dbc66d6a66cbc6d2ef95783e90203d9b353135a0` | #189 | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40502544641) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40502544644) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40502544692) |
| [#192](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/192) | `f93bf0c110f8c0d859e959be27448e780a79e4c6` | #191 | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40707037910) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40707038007) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40707037939) |
| [#193](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/193) | `39a5484119d688a3658c7e5f7d359704ee274866` | #191 | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40811304880) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40811303761) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40811304064) |
| [#194](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/194) | `4825ed62bf94b86bdb9afe5b3896e97ea46c2cea` | main | [failure](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40812112759) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40812112762) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40812112844) |
| [#195](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/195) | `a64434962c3f497a83188b46f8b2f488a353cd8e` | main | [failure](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40819941637) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40819941640) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40819941635) |
| [#196](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/196) | `28bf02ac11a86a52afe4fce029c3961d32dc9f8d` | #192 | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40841961937) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40841962068) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40841961973) |
| [#197](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/197) | `57d8097eb96f431a543ac4d2fba94c513244fa7e` | #196 | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40855845842) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40855845649) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40855845874) |
| [#198](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/198) | `e9abd5803f60ff07b8606f3dfc5900116d44c022` | #192 | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40992858966) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40992858951) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40992858997) |
| [#199](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/199) | `3950c7e77834c28ebde2a09b7b079809a0eacce6` | #192 | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40885233340) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40885233337) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40885233342) |
| [#200](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/200) | `0ce4ec1ef46eb445a4ae0c2f5b06fb66b0b32be2` | #197 | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40872357034) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40872357192) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40872357031) |

上述普通 PR CI 中，real-execution、platform-real-execution、container-smoke、
bulk-rnaseq-real-execution 均为 skipped。它们不记为通过；人工执行的分层证据
另列于下文，不能把某条文档 PR 的绿色 CI 当作重新执行过科学运行。

所有 PR 当前对各自 base 的 GitHub `mergeable` 为 true；#194、#195 的
`mergeable_state` 为 unstable，其余为 clean。这不是所有分支组合后的
冲突或兼容性证明，本次没有执行合并演练。

## 拟议顺序与依赖

| 阶段 | 范围与依赖 | 下一步 | 退出条件 |
| --- | --- | --- | --- |
| 1 | #190，共享 CI 与 cleaner 测试维护，base main | 独立审阅后优先纳入集成基线 | 人工批准；最终候选包含对应修复；不放宽生产 Python 准入 |
| 2 | #184 → #185 → #186 → #187 → #188 → #189 → #191 | 按父到子处理科学计划与 SCI-1～6，实现审阅与验收证据随提交保留 | 每条 PR 的准确差异、依赖、CI 和适用定向证据复核；获得合并授权 |
| 3 | #193 与 #192 都依赖 #191；前者提供执行 Gate 的维护修复，后者提供 SCI-7 准备材料 | 先确认两者均进入候选，再消费后续验证分支 | 不遗漏 #193 的工具链、profile 与容器身份接线；保留 #192 的 manifest/参考边界 |
| 4 | #196 → #197 → #200；该链起点依赖 #192 | 收拢 CUT&Tag 对照、扩展与组蛋白结果及比较脚本 | 保留输入、运行源码和比较口径；无结论仍可交付，不将峰数变化等同于优劣 |
| 5 | #198、#199 各自依赖 #192；与阶段 4 没有互为父分支的关系 | 纳入 MNase 结果与 skills；交叉检查 README、技能引用及最终默认行为 | 保留 MNase 原登记指标未评估；skills 在组合后的实际源码上仍准确 |
| 6 | #194、#195 各自 base main，现有 CI 失败 | 在具备 #190 修复的基线上另行更新并复验；分别审阅操作指导与审计建议 | 新准确 head 的 CI 成功；审计建议不冒充已确认缺陷或已实现功能 |
| 7 | 当前独立文档分支 | 更新状态入口与本清单，准备独立文档 PR | 链接、范围和历史记录完整；获得推送/创建 PR 授权 |
| 8 | 以上获准纳入的最终组合版本 | 形成一个明确候选，执行下节集成验收 | 精确身份、必要检查和证据齐全后，另行决定发布或默认冻结 |

这是依赖顺序，不是要求一次合并全部 17 条 PR。#194/#195 的内容或 CI 未收口时，
可以继续独立处理，不以红灯分支阻塞其他已获准的交付，也不忽略它们的失败。

堆叠 PR 的父提交需要保持可追溯。具体合并方式由维护者决定；若选择 squash
或重写历史，必须另行授权并安排子分支重建与复验。父 PR 进入 main 后，逐个确认
子 PR 的 base、祖先关系和剩余 diff，避免直接将子 PR 合回父分支后就声称已交付 main。
本次不修改任何 PR 的 base、Draft 状态、历史或保护规则。

## #194/#195 的已知失败

[#194 CI](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40812112759) 与
[#195 CI](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40819941637) 均显示：

- shard 1 的 24 项 bamCoverage fixture 在 `conda info --json` 处失败，
  日志为 `conda: command not found`；
- shard 2 的 4 项 deployment 测试找不到 `/usr/bin/python3.12`；
- browser 的 schema/ENCODE cores 前置 smoke 因 conda/profile 失败；
- fast-checks 与 coverage 的红灯由其生产者失败向上传递，不能作为独立根因。

两条 PR 的差异分别仅为 AGENTS.md 和 9 个审计文档，没有包含 #190 的修复。
#190 已增加受支持 runner/Python 检查与 tiny Snakemake profile 隔离。
因此优先复用共享基线，避免在文档 PR 再复制一套维护修复；能否消除全部失败，
仍以更新后的实际检查为准。本轮没有修改这两条 PR 或替它们重新发起 CI。

#194 在 AGENTS.md 中新增 Cloud 环境说明，其环境假设和可复现性仍需内容审阅；
#195 的建议进入候选清单后逐项核实，不直接启动测试删除或架构改造。

## 已执行证据与未完成的科学范围

| 范围 | 已有证据 | 当前边界 |
| --- | --- | --- |
| SCI-6 合成平台与真实工具集成 | 已报告 platform 10 passed、scientific 170 passed、零失败和 skip；原始交接记录保留 | 历史运行记录，不是本次重新执行；不自动覆盖最终组合版本 |
| #193 Protected Bulk Gate | [run 40832642755](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40832642755)，attempt 1；准确 head `39a5484119d688a3658c7e5f7d359704ee274866`，本轮再次查询结论 success | 正式 Gate 已通过；无需继续报告旧的镜像/runner 欠验，也不追认其他身份 |
| #196/#197 CUT&Tag | 对照与扩展记录、真实运行及比较脚本；当前普通 CI 成功 | inconclusive；缺少足以冻结默认的独立确认，不以普通 CI 代替科学裁决 |
| #200 组蛋白 H3K36me3 | retain/exclude 两臂及共同 universe 指标已记录；当前普通 CI 成功 | 单组输入；无默认冻结 |
| #198 MNase Rep1 | 三条 arm 退出 0；额外诊断 `6,883/4,083,944` 已正确标注为 200-arm 基准 | 原登记的 500-arm 短模板指标未评估；no conclusion；没有 confirmation |
| #199 workflow skills | 真实 validate/dry-run、拒绝未标记输出目录；Codex 目录发现和三个显式调用、Cursor 三个显式调用已核对 | 仅显式调用范围通过；Cursor 交互式菜单与隐式自动触发未测试，不作为本轮阻塞 |

以上科学批次的实际运行源码、输入、参考和工具身份以各自结果记录为准，
文档提交的 head 不能替换历史运行 SHA。默认冻结需要单独证据与裁决；
允许 inconclusive/no conclusion，不通过新增阈值或重复挑选输入强求结论。

剩余 TF、ATAC discovery、预先定义的 confirmation、必要的重复一致性与最终科学汇总，
继续沿 [比较设计](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/blob/f93bf0c110f8c0d859e959be27448e780a79e4c6/docs/development/sci7-public-data/comparison-design.md)
分批推进。本清单没有授权新增数据下载、重跑全部科学矩阵或默认冻结。

## 最终候选的集成验收

1. **锁定组合身份。** 记录最终 commit、干净工作树、准确导入位置、锁定工具和参考身份。
   在任务本地环境通过 canonical bootstrap 验证；保护已有运行 outdir，不复制或手写
   marker 绕过身份。文档、科学结果和 skills 的同名文件按内容集成，不能整文件覆盖
   其他已验收分支的结果。
2. **检查真实差异。** 比较候选与已验收实现，特别检查 #190/#191 中共享维护补丁、
   #193 的工具链与容器 UID/GID 接线、比较脚本及 QC/artifact 消费链。
   重新核对堆叠 PR 去掉祖先后剩余的 diff，保留实际参与者署名。
3. **执行对应层验证。** 最终候选的 CI、Lint、Lock Check，以及生产代码影响所需的
   validator、DAG/tiny execution、适用 runtime bundle 和 platform/scientific 集成。
   比较脚本修改执行针对性测试；skills 若指令依赖的行为改变，补相应显式调用验证。
   不因合并了一份文档而机械重跑全部真实文库。
4. **裁定 Gate。** 对比已通过的 #193 身份与候选的 pinned runtime、工具、参考、
   execution closure 和科学 artifact/QC 语义。只改文档/skills 等且闭包未变时，
   不自动重跑 Bulk；相关闭包变化则按既有准入与保护入口执行候选的准确身份。
   历史 Gate 保留原身份，不能改写成候选已通过。
5. **保存与交接。** 保存准确 SHA、run/attempt、JUnit/命令与退出码，分开记录成功、
   失败、取消、skip 和未执行。清理本任务进程后再验收；代码集成、科学结论、
   默认冻结、合并和发布分别裁决。

本轮交付为可审阅的顺序与验收清单。最终组合版本尚未构建或运行，历史证据
未重写；本清单不代表合并批准、发行验收或新增计算授权。

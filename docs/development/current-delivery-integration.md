# 当前 PR 收口与集成清单

## 集成候选（2026-11-06）

本节是当前入口。下方 09:50 快照和更早正文是历史证据，不描述本候选。
main 仍是 `fa0d8401cbc3fecb765fed6d7e01e18bb8fbc96d`。本候选没有合入 main，
没有发布，也没有冻结科学默认。

来源合并树是 `526e26d1ff8dc6e39807fae9c9406f9bd0a19e7e`。
技能启动说明对齐是 `b436087176f2c82d7c01f30dca1434eaecab3349`。
本节首次写入于 `954014fd01e8ed840c0279ffec279718ea42831b`。
评审身份是 Draft PR head 的完整 SHA，不把任一来源 PR head 写成候选身份。
候选 CI 是该 head 上新跑的普通 CI、Lint 和 Lock Check，
不沿用下表里来源 head 的结果。

### 纳入

祖先核对后，子提交已经包含的父提交没有再次叠补丁。#190 的两个 CI 提交
（`4e7c5a46945f69a703d0d3ec0d90ac8ebc7d599c`、
`b6f3a49f45a6b8093f8c42db6684bce5f1425475`）与科学链上的
`978651a45974216066488449f49e43dced5bb40e`、
`50c28743dfa3693ec90de3d2bd6106cbdbf4ed8a` patch-id 相同，只保留一次。
#190 独有的清理预算 `32570adf50725206e0bfe40a6af1d3d81391478d` 已进入候选。

| PR | 准确 head | 处理 |
| --- | --- | --- |
| #190 | `b6f3a49f45a6b8093f8c42db6684bce5f1425475` | 快进纳入 |
| #184–#192、#196、#197、#200、#202 | 见 09:50 表，本轮复核未移动 | 由 #203 祖先包含 |
| #203 | `65760720c4609de633ef9f636ab9d4eea4f3e8c0` | 合并；普通 CI [42154760414](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/42154760414)、Lint [42154760469](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/42154760469)、Lock Check [42154760405](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/42154760405) 均为 success |
| #193 | `39a5484119d688a3658c7e5f7d359704ee274866` | 合并。PR CI [40811304880](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40811304880) attempt 1 success。Protected Bulk Gate 仍是 dispatch [40832642755](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40832642755) attempt 1 success，只属于该 head |
| #198 | `e9abd5803f60ff07b8606f3dfc5900116d44c022` | 合并 |
| #199 | `3950c7e77834c28ebde2a09b7b079809a0eacce6` | 合并 |
| #204 | `ae08bd2544ae47a3691cd30ce299ebaaf96c5a68` | 合并 |
| #205 | `be4a2db56811adc1e9702382234cdddb5819b998` | 合并 |
| #206 | `ce86c52cc489825e04e52cef6b268e0a678ceba7` | 合并 |
| #201 | `ce10b73717af4c602b3210a0ebd6451ed52b50db` | 合并。该 head 的 CI [42146755914](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/42146755914) 仍是 failure；Lint 与 Lock Check 成功。失败留在来源 head，不记成候选通过 |

### 排除

#194 `4825ed62bf94b86bdb9afe5b3896e97ea46c2cea` 只改 `AGENTS.md` 的 Cloud Agent
工具链说明。#195 `a64434962c3f497a83188b46f8b2f488a353cd8e` 只新增九份工程审计笔记，
这些文件只存在于该 PR，不在本候选中。二者都不改变科学行为、
工作流、CI、测试、锁文件或 API。它们各自的 CI 仍失败。本候选不纳入，避免把
Cloud 运行手册和一套并行的代理流程路线图写成当前交付。后续应在候选之后单独
审阅，并先核对已过时的环境、Alembic 和 Actions 陈述。

### 冲突处理

- `docs/development/sci7-public-data/README.md`：保留各 assay 已交付记录，并写入
  MNase Rep1 的 no conclusion。不再写“MNase 记录留在另一分支”。
- `workflow/Snakefile`：使用 `import_checkout_package`。包已定位但初始化失败时
  停止，不退回另一个 checkout。SCI-1～6 的 insert bound、policy revision、
  run identity 和 effective config 导入全部保留。
- `test/config/test_validation.py`：删除 `test_mnase_se_still_rejected` 和
  `test_mnase_narrow_still_rejected`。保留同断言的 `test_mnase_se_rejected`、
  `test_mnase_narrow_peak_mode_rejected`，以及 `test_scientific_policy_revision_key_rejected`。

#193 的 conda 工具链、profile 隔离和容器 UID/GID 接线在合并后与 gate head 的
`.github/workflows/ci.yml`、`workflow/envs/ci-conda.lock`、
`workflow/envs/ci-conda.yml`、`containers/Apptainer.ci-conda.def` 字节相同。

### 闭包与 Gate

相对 `39a5484119d688a3658c7e5f7d359704ee274866`，119 个 bulk 实现路径和
`execution-implementation-manifest-1.0.0.json` 字节相同，aggregate 仍为
`15a16a23f8a97afa006b475a9a9be741ba84750ddedd81f96ca10326383c5b86`。
`default-execution-qualification-1.1.0.json` 也相同。
`scripts/source_provenance.py` 和 `workflow/Snakefile` 有变化，且不在这 119 个路径里。
因此本候选不重跑 Protected Bulk Gate，也不把历史 Gate 写成候选已通过。

### 科学边界

SCI-7 本批测量和记录已交付，结论是 no conclusion。JUND 是跨靶点检查，不是
NR2F2 的生物学重复。CUT&Tag、histone、MNase、ATAC 的 confirmation 未执行。
MNase 登记的 500-arm 短模板指标未评估。默认未冻结。

## 历史快照与执行计划（2026-11-06 09:50，Asia/Shanghai）

本节在写入时取代当时的历史“当前”和“下一步”。它现在本身是历史，原样保留。
此次重新查询了 #184～#206 的开放状态、准确 head 的普通工作流结果，
并读取 #201～#206 的 CI job 列表；本地科学状态只读核对，不接管原 session。

main 仍为 `fa0d8401cbc3fecb765fed6d7e01e18bb8fbc96d`。
23 个 PR 全部 open、未合并：#194/#195 非 Draft，其余 21 个为 Draft。
17 个 head 的 CI、Lint、Lock Check 成功；#194、#195、#201、#204、#205、
#206 的 CI 失败，Lint 与 Lock Check 成功。失败不能因“共享基线”而记为通过。

### 准确 head 清单

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
| [#201](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/201) | `4bf6791b38146bc5348dab22977dccd684d5a444` | main | [failure](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40996506173) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40996506078) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40996506069) |
| [#202](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/202) | `634160e0d7c671db466a22a24bbe2287782b9af7` | #200 | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/41111755339) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/41111755211) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/41111755539) |
| [#203](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/203) | `4ce8c6ee7384083b2b6a8e191dd8e383f76437fe` | #202 | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/41373058173) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/41373058168) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/41373058171) |
| [#204](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/204) | `ae08bd2544ae47a3691cd30ce299ebaaf96c5a68` | main | [failure](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/42009798710) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/42009798683) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/42009798739) |
| [#205](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/205) | `be4a2db56811adc1e9702382234cdddb5819b998` | main | [failure](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/42028980114) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/42028979917) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/42028980154) |
| [#206](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/206) | `ce86c52cc489825e04e52cef6b268e0a678ceba7` | main | [failure](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/42144947651) | [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/42144947720) / [success](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/42144947667) |

#201 行记录本次修改前的远端 head；本次文档更新形成的新 head 不继承该行 CI。
#201～#206 的四项条件 job（real-execution、platform-real-execution、
container-smoke、bulk-rnaseq-real-execution）均 skipped，不是实际执行通过。
更早 PR 的执行分层证据保留在下方历史记录。

### 工作线与验收重点

| 工作线 | 范围 | 依赖 | 下一步 | 验收标准 | 关联 PR |
| --- | --- | --- | --- | --- | --- |
| 科学批次收口 | 汇总已执行 discovery、TF 跨靶点确认与尚未评估项 | 原科学 session 的比较结果、登记规则、输入/运行身份 | 完成 JUND 三臂后的 duplicate 与 max_insert 比较、审阅并交付结果；随后汇总各 assay 的限制与未评估项 | 比较口径与登记一致，缺失项显式保留；允许 no conclusion；执行成功不等于默认获准 | #192、#196～#198、#200、#202、#203；TF confirmation 尚无已核实的交付 PR |
| Workflow skills | 已验收显式调用的维护与集成适配 | #199、最终源码及 #204 启动行为 | 集成时核对运行前检查、解释器选择和 profile 指令；只按实际差异补验证 | 区分独立 bootstrap 检查与 Snakemake 进程内检查；不承诺完整工具链身份 | #199、#204 |
| 文档治理 | 单一当前入口、准确状态与历史证据 | 各 PR、CI 与科学比较交接 | 在 #201 更新当前快照、依赖和验收清单；保留旧记录 | 23 条 PR 的分支交付与 main 状态可区分；不倒改历史运行 SHA | #201；#194/#195 独立审阅 |
| 测试精简 | 首批去重与第二批构建复用收口 | 保留契约、隔离证据、共享 CI 基线 | #205/#206 保持范围；在更新基线后验收新 head，不自动启动第三批 | 断言不弱化、可独立执行；耗时收益限于实际测量区域 | #205、#206；#195 为审计来源 |
| 开发方式优化 | 任务范围、独立 worktree、写入边界及交接 | 并行 session 经验、#195 建议 | 将已采用的约定纳入交接；后续独立评估是否还需文档变更 | 无重复计算或交叉写入；源码、进程、产物与负责人可追溯 | #195、#201、#204 |
| CI 与 Gate 运维 | 共享基线修复和最终组合验收 | #190、#193、候选真实闭包差异 | 优先审阅 #190 的纳入方案；得到集成授权后更新受影响分支并复验 | 新准确 head 的必要 CI 成功；按闭包决定 Gate，历史通过不追认新身份 | #190、#193、#194、#195、#201、#204～#206 |

### 科学计算的最新边界

- TF discovery（NR2F2/HepG2）已交付 #202；ATAC Treg discovery 已交付
  #203。两个准确 head 的普通 CI 成功，科学结论均为 no conclusion。
  [ATAC 记录](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/blob/4ce8c6ee7384083b2b6a8e191dd8e383f76437fe/docs/development/sci7-public-data/atac-treg-discovery.md)
  明确保留同供体、26 bp 读长与 confirmation 未评估的限制。
- TF confirmation 的本地 `run-state.json` 已为 `ALL_ARMS_DONE`：
  exclude-x2000、retain-x2000、exclude-x500 均记录 exit 0。
  最后一臂结束于 2026-11-05 23:52（Asia/Shanghai），日志为 45/45。
  运行源码为 `4ce8c6ee7384083b2b6a8e191dd8e383f76437fe`，
  不能用以后结果文档的提交替换该身份。
- 09:50 只读核查时，比较目录已有共同峰区间，max_insert 日志正在记录
  500-arm 模板扫描；尚未取得最终比较报告。此状态不证明比较成功完成，
  也不据日志更新时间断言某个进程仍在运行。由原科学 session 完成比较和交付。
- JUND/K562（CEBPB deletion 背景）相对 NR2F2/HepG2 是跨靶点稳健性检查，
  不是同一靶点的独立重复。三臂完成后不自动新增全流程运行。
- CUT&Tag、H3K36me3、MNase Rep1 的已有结果保留。MNase 原登记的
  500-arm 指标仍未评估；ATAC confirmation 等未执行项在科学总表中明确列出，
  不把候选登记写成已经完成，也不为强求默认冻结自动扩张矩阵。
- 用户此前对范围内科学下载/计算的授权继续有效，不重新逐文件索取授权。
  本次计划更新只读核对科学任务，不启动、停止或清理它们。默认冻结单独裁决。

### 两批测试与启动修复

- #204 `ae08bd2`：Snakemake 在自身进程检查包来源，仅顶层包无法定位时
  才 fallback；包初始化 ImportError 明确失败。已补实际 Snakemake 负例，
  本轮定点审阅通过；这不是完整工具链、子进程或科学验收保证。
- #205 `be4a2db`：只删除两项与保留案例输入/断言相同的 MNase 校验；
  本地 41 → 39 passed。不声称显著性能收益。
- #206 `ce86c52`：模块内各构建一次 wheel/sdist，各测试使用独立副本，
  安装、venv 和子进程断言保留；PEP 517 构建仍独立执行。9 项测试均在
  [本次 shard 1](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/42144947651/job/127654044537)
  中通过。作者报告三次同环境测量约减少 3.6 秒，仅属于该模块，
  不外推为整个 shard 加速。本轮代码审阅与 diff 检查通过。

#206 主 CI 已从排队/运行中变为 failure：shard 1 为 24 failed / 3766 passed，
24 个失败均位于 bamCoverage fixture，日志为 `conda: command not found`；
shard 2 为 4 failed / 4027 passed，缺少 `/usr/bin/python3.12`；
browser 前置 smoke 同样因 conda/profile 失败，fast-checks 与 coverage 连带失败。
#201 历史 head 的本次日志复核也显示同类 24/4/1 失败。
#204/#205 既有失败核对结果与此一致；不要复制维护补丁或通过删除测试规避。

### 下一步顺序与交接

1. **并行收口：** 原科学 session 完成 TF 比较及其 Draft PR；文档 session
   维护本清单。#204/#205/#206 不再自动扩大范围。
2. **共享基线优先：** 审阅 #190 的纳入方案。实际合并或分支集成仍需授权；
   当前未做 merge、rebase 或 cherry-pick。基线变化后重新验证红灯 PR 的新 head，
   不把现有红灯改记为历史噪声。
3. **按依赖形成候选：** 延续下方 #184→#191、#192/#193 的顺序。
   科学结果链扩展为 #192→#196→#197→#200→#202→#203；
   #198/#199 仍分别依赖 #192。#204/#205/#206 均直接基于 main，彼此没有
   堆叠依赖，#206 不包含 #205 的删除。最终组合时核对 #204 的 Snakefile
   启动检查与 SCI-1～6 的科学身份逻辑，不能整文件覆盖。
4. **组合验收：** 在明确候选 SHA 上执行影响范围要求的检查，复核 skills、
   比较脚本和源码身份。不因文档或测试精简机械重跑全部真实文库。
   #204 修改共享 provenance 文件，相关 runtime/资格闭包是否受影响须实际比较。
5. **Gate 与交付：** #193 的 Protected Bulk Gate
   [run 40832642755](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/actions/runs/40832642755)，attempt 1，
   本次查询仍为 success，head 为 `39a5484119d688a3658c7e5f7d359704ee274866`。
   它是已完成的历史正式验收；最终候选是否需再运行按闭包差异决定。
   合并、发布和默认冻结保持分别裁决。

本次计划更新没有运行测试套件、科学流程或 Gate，没有修改其他任务分支。
以下 2026-10-30 快照作为历史保留。

## 历史快照（2026-10-30）

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

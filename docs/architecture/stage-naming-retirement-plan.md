# 重命名 stage 时代遗留命名（已批准）

> 状态：**已批准执行**。执行日期：2026-09-21。
> 已定决策：硬改不留别名；新命名与前端 UI 已有标签对齐。

## 背景

项目按 Stage 1/2/3/4b/5/53+ 迭代开发，遗留了一批阶段编号命名。这些名字已进入公共契约（config schema、output-contract、manifest 产物词汇表），需要一次性清除。

## 重命名映射

| 旧名 | 新名 | 类型 |
|------|------|------|
| `stage4b` | `replicate_analysis` | config 键（约 156 处 / 50 文件） |
| `stage5` | `chipseq_idr` | config 键（约 170 处 / 56 文件） |
| `stage3_qc_summary` | `project_qc_summary` | 规则名 + 输出文件 `results/multiqc/project_qc_summary.tsv` + manifest output_type（28 处 / 13 文件） |
| `stage5b_summary` | `chipseq_idr_summary` | 规则名 + 脚本文件名 + 日志名 `<exp>.chipseq_idr.summary.log`（18 处 + MANIFEST.in） |

防撞名已确认：`chipseq_idr_summary` 未被占用（`idr_reproducibility_summary` 及其三个 wrapper 是另一套体系，不动）。

## 执行步骤（按依赖顺序）

### 1. workflow/（规则与 schema）
- `workflow/Snakefile:66-67,134`：`VALIDATED_CONFIG.get("stage4b")` → `"replicate_analysis"`，`"stage5"` → `"chipseq_idr"`，变量名 `STAGE4B`/`STAGE5` → `REPLICATE_ANALYSIS`/`CHIPSEQ_IDR`
- `workflow/schemas/config.schema.yaml` + `config.schema.json`：改键名、更新 description 里的旧引用（yaml/json 双格式手工同步，无生成器）
- `workflow/rules/qc.smk:912-925`：规则改名 + 输出路径 + 头注释（qc.smk:17）
- `workflow/rules/idr.smk:397,412,419`：规则改名、日志名、脚本路径、头注释（idr.smk:14）
- `workflow/rules/targets.smk:153`：目标路径
- `workflow/rules/consensus.smk:130,414`：注释里的 "legacy stage5 IDR" 措辞
- 顺带清理 config.yaml 和各规则文件头注释中的 "Stage 3/4b/5/53+" 阶段语言，改为语义描述

### 2. scripts/
- `git mv scripts/stage5b_summary.py scripts/chipseq_idr_summary.py`，更新其 docstring usage 自引用
- `scripts/aggregate_qc_summary.py:4,10`：docstring 中的规则名/路径引用
- `scripts/prepare_public_validation_inputs.py` + 对应测试：`stage5_idr` 列名 → `chipseq_idr`
- `MANIFEST.in:17`：打包清单同步新脚本名

### 3. src/（平台代码）
- `config/validator.py`：键名校验 + 跨键约束错误消息（"stage5 requires stage4b" → "chipseq_idr requires replicate_analysis"）
- `config/reproducibility.py:145-150`：`chipseq_narrow: null` 推断来源改读 `chipseq_idr`
- `adapters/encode.py:209-210`：映射元组简化为同名映射
- `samples/replicates.py`、`samples/loader.py`、`cli/_validator.py`、`cli/results_visibility_fixture.py`：键名与错误消息
- `manifest/make.py`：`stage4b`/`stage5` 键（387,396,749-758,840,880,908,931,940,1382-1405）+ `stage3_qc_summary` output_type（1324-1336）
- `artifacts/artifact-inventory.yaml` + `docs/architecture/artifact-inventory.yaml`（镜像双份，行号一致，同步改）：producing_rule、path_template、manifest_output_type、config_gate 文本

### 4. config/config.yaml
- `stage4b: true` → `replicate_analysis: true`，`stage5: false` → `chipseq_idr: false`，注释语义化
- 注意：这是本地已修改文件，改完要保证 validate 通过

### 5. docs/（用户向文档）
- `output-contract.md`（stage3 行 149、stage5b 行 94-96、stage5 行 86/191、stage4b 行 181-194）
- `configuration.md`（stage4b 150-253、stage5 151-255、stage3 115）
- `reproducibility-policy.md`（17 处，含 "public configuration keys" 声明）、`idr-contract.md`、`assay-policy.md`、`sample-sheet.md`、`qc-interpretation.md`
- **不动** `docs/release-checks/`（历史发布记录）和 `docs/operations/snakemake-lint-warnings.txt`（lint 快照，如需可重新生成）

### 6. test/（约 30 文件）
- 所有引用的键名/规则名/产物名更新（重点：test/config/、test/workflow/、test/manifest/test_make_manifest.py、test/adapters/、test/samples/test_replicates.py 含错误消息文本断言、test/api/、test/services/、test/browser/、test/profiles/*/config.yaml ×9）
- DAG 快照重生成：`python3 -m pytest test/test_dag_snapshots.py --update-snapshots`
- 前端负向断言 `/stage4b|stage5/`（frontend/e2e + src/routes/__tests__）语义仍然成立（UI 本就不显示旧名），但需确认断言文本不被新名误伤

### 7. ADR
- 按 AGENTS.md 契约规则，在 `docs/architecture/` 加一篇短 design note：记录配置键重命名（破坏性变更）、新旧名映射、迁移方式（用户改 config 键名即可）

## 验证（全部通过后才算完成）

1. `./.local/envs/ci-fast/bin/python -I -S scripts/checkout_bootstrap.py --repository-root . validate --config config/config.yaml`
2. `snakemake -s workflow/Snakefile --configfile config/config.yaml -n`（DAG 构建）
3. `python3 -m pytest test/config test/workflow test/manifest test/adapters test/samples test/api/test_routes_workflows.py test/services/test_execution_planner.py test/test_dag_snapshots.py test/browser/test_results_visibility_fixture.py -q`（用 ci-fast 环境）
4. `python3 -m ruff check <改动文件>` + `ruff format --check`
5. `npm --prefix frontend test -- --run`（如前端断言受影响）
6. `git diff --check` 干净
7. 全仓 grep 确认生效代码、配置、schema、产物路径、规则名和 output_type 中旧名零残留；报告按“生效引用零残留 + 例外逐项清单”给出证据。

用户于 2026-09-21 确认的合法 grep 例外：

- `docs/release-checks/`：冻结的历史发布记录。
- `docs/operations/snakemake-lint-warnings.txt`：冻结的 lint 快照。
- 本计划与本次新增 ADR：保留旧→新映射的迁移记录。
- 前端负向断言：继续断言旧名不出现在 UI 或提交内容中。
- `docs/development/bug-log-2026-09-local-bring-up.md`：历史 bug 记录。
- 既有后端 schema、渲染结果及错误脱敏断言：保留旧名不出现的检查。

2026-09-21 后续决定：平台尚未正式投产，配置随代码同 PR 迁移；保留 validator
忽略未知键（含旧键）的既有行为，不新增旧键拒绝检查或拒绝测试。撤销此前相应例外，理由见本次 ADR。

## 执行方式

触点约 350 处 / 60+ 文件，为机械化但需保持一致性的重构。由 coder 子代理按上述步骤执行，主会话负责复核 diff 与验证结果。

## 明确不做

- 不改 `stage_bulk_rnaseq_runtime_assets.py`（"stage" 是动词"暂存"）
- 不改 `docs/release-checks/` 历史文档
- 不留任何向后兼容别名/垫片
- 不动 frontend/src/api/generated/（本次零命中，无需 regenerate）
- 不动本地未提交的 `workflow/rules/common.smk`（fastqc_trimmed）

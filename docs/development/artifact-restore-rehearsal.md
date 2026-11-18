# 单份结果本地恢复演练（IRENE-26）

## 结论与适用边界

当前可以恢复的是**结果字节，并在原布局重新建立结果索引**。普通文件复制会改变描述符身份；不能把同字节副本放回路径后仍声称旧 `artifact_generation` / `artifact_revision` 有效。本演练实测：

- 缺原文件时，授权下载明确返回 `409 RUN_ARTIFACT_DOWNLOAD_CONFLICT`；
- 同字节副本放回原路径后，旧身份仍返回 409，原因是 `ARTIFACT_DOWNLOAD_REVISION_MISMATCH`；
- 复用 `ArtifactExtractionService.extract` 重新提取后，`run_id` 和确定性的 `artifact_id` 不变，generation 和 descriptor-bound revision 更新；
- 旧 generation 返回 409，新身份的全文下载、WashU data/index Range 读取成功；
- SQLite 副本中的测试会话仍可授权读取，匿名请求仍返回 401。

这是一条内部服务级恢复路径，不是现成的一键管理员恢复命令。没有修改生产 schema、下载接口或文件解析规则，也没有把移动到新长期盘写成已有能力。

## 恢复前必须列清的对象

| 对象 | 本演练使用的身份或位置 | 恢复语义 |
| --- | --- | --- |
| 逻辑结果 | `run_id`、`artifact_id`、generation、revision | run/artifact ID 保持不变；重建索引后 generation/revision 必须重新读取 |
| 原始字节 | `workspaces/<run_id>/results/browser/sample_peaks.narrowPeak` | 备份与恢复副本均核对 SHA-256；本 fixture 只有两条合成 narrowPeak 记录 |
| SQLite | 独立 `live/platform.db` 与 `backup/platform.db` | 停止 app 后用 SQLite backup API 创建副本，`PRAGMA integrity_check` 必须为 `ok` |
| 结果索引元数据 | SQLite 中的相对路径、大小、output type、revision、`_platform_path_identity` | 不用手工 SQL 改写；通过既有提取服务重新生成 |
| WashU 衍生物 | `browser-derivatives/<run_id>/<artifact_id>/<revision>/track.gz[.tbi]` 与 provenance | 备份并核对摘要；恢复后用新 revision 重新准备，不把旧衍生物当成新身份 |
| 读取权限 | SQLite 中的合成管理员会话 | 新 app 不重新播种账号；授权请求成功，匿名 data Range 请求为 401 |

Artifact publication 是结果元数据投影，不复制或恢复原始字节；历史 publication 也不会让旧 generation 重新成为当前结果。部署契约中的 root-owned 数据库备份、writer-stop witness 和 pre-PONR 恢复属于另一条生产操作边界，本演练不模拟或替代它们。

下载接口本身声明 `Accept-Ranges: none`。这里的 Range 验收走 WashU byte 路由：准备好的 `track.gz` 和 `track.gz.tbi` 均以 `Range: bytes=0-3` 请求，期望 206、正确的 `Content-Range` 和 BGZF 文件头。

## 可复现命令

在仓库根目录使用任务专用的锁定环境。环境不存在时，先按 [Local Runtime](local-platform-runtime.md#run-locally) 的 `ci-fast` 说明创建，不得在共享环境安装临时依赖。

```bash
./.local/envs/ci-fast/bin/python -I -S scripts/checkout_bootstrap.py \
  --repository-root . pytest test/services/test_artifact_restore_rehearsal.py -v
```

pytest 的 `tmp_path` 提供独立临时目录；测试只导入本 checkout，复用 `test_artifact_extraction.py` 的合成 adapter，不启动 Redis、worker、Snakemake、服务器或科学计算，也不读取真实数据。

## 演练步骤与实测检查点

`test_single_artifact_restore_requires_original_layout_and_reindexing` 在同一绝对 `live/workspaces` 布局执行以下步骤：

1. 创建一个已成功结束的合成 run，写入两条 narrowPeak 记录，通过真实 `ArtifactExtractionService` 建立结果索引。
2. 记录 run/artifact/generation/revision 和原文件 SHA-256；经 HTTP 下载原文，准备真实 BGZF/tabix 衍生物，并分别验证授权 data/index Range 为 206。
3. 关闭 app，使用 SQLite backup API 创建隔离数据库副本并检查完整性；复制原文件、data、index、provenance，逐份核对副本摘要。
4. 删除 live SQLite 和原文件，只把 SQLite 副本放回原数据库路径。重开 app 后，数据库中的旧 generation/revision 与备份前一致，但原文件下载返回 409，reason 为 `ARTIFACT_DOWNLOAD_SOURCE_MISSING_OR_UNSAFE`。
5. 把同字节原文件副本放回原路径，确认 SHA-256 相同；旧身份下载仍返回 409，reason 为 `ARTIFACT_DOWNLOAD_REVISION_MISMATCH`。这一步防止把“文件存在且字节相同”误报为旧结果身份已恢复。
6. 再次调用既有提取服务。断言 run/artifact ID 不变，新 generation/revision 与旧值不同，旧 generation 下载返回 `ARTIFACT_DOWNLOAD_GENERATION_CHANGED`。
7. 用新身份下载原文并核对 SHA-256，重新准备衍生物，再核对 data/index Range 与演练前一致。最后以无 cookie 请求确认 Range 读取仍为 401，且错误响应不泄露临时路径。

因此，恢复验收必须同时包含**字节摘要一致**和**新结果身份下的授权读取成功**；只检查文件存在、只检查 SQLite 行存在，或只检查旧衍生物仍可读取都不充分。

## 最小下一步

增加一个经审阅的管理员“从原布局重建结果索引”入口，复用 `ArtifactExtractionService` 的 attempt/generation 并发保护，明确返回新的 generation/revision 并促使客户端刷新。不要把直接 UPDATE SQLite 当作恢复工具。

若后续需要迁移长期盘，应另行设计内容寻址的结果字节备份、完整性清单和重新绑定流程；本次同路径副本演练不为该能力背书。

## 未覆盖项

- 真实长期盘、跨文件系统/跨机器迁移、整 workspace 或多份结果的一致性恢复；
- 生产 systemd writer-stop、root-owned backup slot、PONR、升级/回滚及并发写入；
- 大文件的全内容哈希、备份截断/掉盘/崩溃窗口、权限或属主变化；
- 真实 GBM 数据、科学结果有效性、真实浏览器和服务器部署；
- 公共管理员恢复 API（目前未提供）及旧页面自动刷新。

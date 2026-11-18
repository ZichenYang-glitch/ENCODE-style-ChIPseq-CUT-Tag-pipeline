# 临时运行盘迁到长期盘的旧链接演练（IRENE-33）

## 结论与适用边界

同一条下载链接可以在复制、核验、切换之后继续读到迁移前封存的字节。链接里的
`generation` 和 `revision` 不改，SQLite 里的 run 与 artifact 身份也不改。切换完成前，读者仍读临时盘上的原文件；切换完成后，读者改读长期盘上已核验的副本，即使临时文件随后被改写或删除。

这不是目录迁移，也不改公开目录 API、Alembic 或生成客户端。`create_app` 的
`relocation_root` 默认是 `None`；不传入时，下载行为与这次改动之前相同。真实批量归档仍要单独确认保留哪些字节、保留多久、以及目标盘是哪一块。本演练不删除临时原件，也不选择生产磁盘。

依据是 IRENE-26 在 `41a4234dce17ac9f7df68cb38d311edd0a484eec` 上的原布局恢复实测，以及 IRENE-29 在 `ab9bbd1a7957defbc98929730c9e0af93767ddac` 上的精确修订读取契约。原布局恢复会更换 generation，所以旧链接不能靠重建索引跨盘继续读。登记契约按精确修订读字节，且不把改路径写成原地更新。这里的做法是：迁移前封存大小、SHA-256 和权限位，复制到另一个目录，核验后再切换；下载坐标保持为迁移前保存的那一组。

## 两个临时目录

| 对象 | 位置 | 迁移语义 |
| --- | --- | --- |
| 临时运行盘 | `scratch/workspaces/<run_id>/results/...` | 只放合成 narrowPeak，不读真实队列或 GBM 原件 |
| 长期盘 | `durable/payload/<run_id>/results/...` | 复制后的常规文件；权限位与封存时一致 |
| 阶段日志 | `durable/journals/<run_id>/<artifact_id>.json` | 只记阶段、逻辑身份、大小、SHA-256 和权限位，不记绝对路径 |
| 保存的链接 | `/api/v1/runs/<run_id>/artifacts/<artifact_id>/download` | 查询参数是迁移前的 `generation` 与 `revision` |

合成正文是两条 narrowPeak 记录，共 50 字节。捕获前把源文件权限设为 `0640`。超过 16 MiB 的文件不进入这条路径。

阶段按顺序是 `captured`、`copying`、`copied`、`verified`、`switched`。`resume` 每次只推进一个未完成阶段。复制中断后日志停在 `copying`，再次 `resume` 会整份重写长期副本。封存时先用下载路径的索引修订核验同一份读取：描述符修订或有界内容修订必须等于 SQLite 里保存的修订，路径身份也必须一致。摘要是这份已核对的字节，不是稍后另读的一份。源文件在封存前被改成等长的另一份正文时，修订对不上，不写阶段日志，也不发布副本。源字节或权限在复制时与封存值不同，则停在 `copying` 且不发布副本。核验或切换发现长期副本不符，则不把阶段推进到 `switched`。

## 验收动作

在仓库根目录使用任务自己的锁定环境。环境不存在时，先按 [Local Runtime](local-platform-runtime.md#run-locally) 的 `ci-fast` 说明创建。

```bash
./.local/envs/ci-fast/bin/python -I -S scripts/checkout_bootstrap.py \
  --repository-root . pytest \
  test/services/test_artifact_disk_relocation_rehearsal.py -v
```

pytest 的 `tmp_path` 提供 `scratch` 和 `durable` 两个临时目录。测试只导入本 checkout，复用 `test_artifact_extraction.py` 的合成 adapter，不启动 Redis、worker、Snakemake、服务器或科学计算。

必须同时看到：

1. 迁移前用授权会话下载，正文 SHA-256 与稍后封存的摘要相同；同一链接的匿名请求为 401。封存前若把源文件改成等长的另一份正文，封存被拒绝且不写阶段日志；这条已保存的链接返回 `ARTIFACT_DOWNLOAD_REVISION_MISMATCH`，响应正文不是改写后的字节。
2. 复制被中断时，长期盘上可以留下半份文件，但保存的链接仍返回临时盘上的 50 字节原件。未切换的阶段日志不会把这份半文件当成下载正文。
3. 复制完成后、切换前，即使长期副本被改坏，保存的链接仍返回临时盘原件。核验失败时阶段留在 `copied`。
4. 从完好的临时原件重新复制、核验并切换后，同一 URL、同一 `generation`、同一 `revision` 返回封存字节。SQLite 中的 generation 与 revision 不变。
5. 切换后改写临时文件、或删除临时文件，成功响应仍是封存字节。改长期副本的权限位返回 `ARTIFACT_DOWNLOAD_RELOCATION_PERMISSION_MISMATCH`；改字节返回 `ARTIFACT_DOWNLOAD_RELOCATION_BYTES_MISMATCH`；把长期路径换成符号链接返回 `ARTIFACT_DOWNLOAD_RELOCATION_UNAVAILABLE`，且响应不含链接目标。
6. 换一个 revision 查询同一路径返回 `ARTIFACT_DOWNLOAD_REVISION_CHANGED`。匿名请求在切换后仍是 401。错误和 401 正文不包含两个临时目录的绝对路径。
7. 阶段日志和实现都不删除临时原件。测试只在切换之后主动删除临时文件，用来确认读者已经改读长期副本。

## 本次实测的读者行为

`test_saved_link_reads_durable_bytes_after_resumable_switch` 在一次进程里保存链接，然后按上面的顺序中断并继续。观察到的结果是：

| 时机 | 保存的链接 |
| --- | --- |
| 封存前，源文件被改成等长的另一份正文 | 409，`ARTIFACT_DOWNLOAD_REVISION_MISMATCH`；不写阶段日志 |
| 复制中断，临时文件仍在 | 200，正文为迁移前的 50 字节 |
| 已复制但未切换，长期副本被改坏 | 200，仍是临时盘原件 |
| 核验因此失败 | 阶段保持 `copied`，读者仍是 200 原件 |
| 重新复制、核验并切换 | 200，摘要与封存值相同，generation 与 revision 未变 |
| 切换后改写或删除临时文件 | 200，仍是封存字节 |
| 切换后只改长期副本权限 | 409，`ARTIFACT_DOWNLOAD_RELOCATION_PERMISSION_MISMATCH` |
| 切换后改长期副本字节 | 409，`ARTIFACT_DOWNLOAD_RELOCATION_BYTES_MISMATCH` |
| 切换后把长期文件换成符号链接 | 409，`ARTIFACT_DOWNLOAD_RELOCATION_UNAVAILABLE` |
| 匿名请求，切换前和切换后 | 401 |
| 使用另一个 revision | 409，`ARTIFACT_DOWNLOAD_REVISION_CHANGED` |

源文件在封存前被改成等长正文时，`test_equal_length_edit_before_capture_keeps_revision_conflict` 看到封存拒绝，原因是索引修订已经对不上这次读到的字节。阶段日志不出现。保存的链接仍是下载路径原来的 `ARTIFACT_DOWNLOAD_REVISION_MISMATCH`，不是长期副本。

源文件在封存后被改写时，`test_scratch_change_before_copy_does_not_publish` 看到阶段停在 `copying`，长期副本没有发布。即便随后放进一份半文件，未切换的日志对下载查询仍是 absent，不会把那份字节交给读者。

## 未覆盖项

- 真实长期盘、跨机器、整 workspace、多份结果或大于 16 MiB 的全内容复制；
- 保留期限、批量归档范围、目标盘选择，以及切换后何时删除临时原件；
- 生产 writer-stop、root-owned 备份、PONR、并发写入和崩溃窗口的掉电证明；
- WashU 衍生物、目录登记表、公开目录路由、Alembic 和生成客户端；
- 真实 GBM 数据、科学有效性和服务器部署。

实现位于 `src/encode_pipeline/services/artifact_disk_relocation.py`。下载服务只在调用方显式传入 `relocation_root` 且阶段已经是 `switched` 时读取长期副本。

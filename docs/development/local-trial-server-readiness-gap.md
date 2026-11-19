# 本地试用到服务器部署的差距

状态：IRENE-37 设计，2026-11-19 按主干核对。核对基线是
`852202096db519a5e6ffdb39418fdb14e371887e`。本文只登记已有证据、尚未执行的
检查和验收条件。它不部署、不改生产主机、不打标签、不发布安装包，也不改
schema、技能或运行配置。

服务器或正式共享试用仍须另一次明确授权。那次授权沿用现有发布清单和部署
决定里的主机窗口，不从本文自动开始。

## 已有清单各自证明什么

| 已有文档 | 已覆盖的范围 |
| --- | --- |
| [发布清单](../../RELEASE_CHECKLIST.md) | 从干净提交做出 v0.4.0 / v0.4.1 安装包、离线 bundle、schema 08 到 15 的升级说明，以及受权主机窗口的命令边界 |
| [本地试用清单](../local-trial-checklist.md) | 精确 v0.4.0 标签在一台 Linux x86_64 systemd 或启用 systemd 的 WSL2 上的试用步骤 |
| [v0.4.0 候选证据](../release-checks/v0.4.0-release-candidate.md) | 日期提交之前的候选包、两次独立安装和 schema 08 到 15 的迁移记录；文中写明当时尚未取得真实主机部署验收 |
| [部署与运行时决定](../architecture/supported-deployment-runtime-decision.md) | 唯一支持拓扑：一台主机、systemd、SQLite、文件存储、独立 ENCODE 与 Bulk 运行时；升级、回滚和备份顺序 |
| [本地运行说明](local-platform-runtime.md) | v0.3.0 `20260717_08` 升到 v0.4.0 `20260827_15` 的备份与迁移步骤 |
| [Peri 运行资格](peri-runtime-qualification/README.md) | 固定 binary 在受限 A1/C4 运行时内被接受 |
| [WashU 分发说明](wuepgg-62.0.2-distribution.md) | `wuepgg@62.0.2` 的许可文本、打包路径，以及商业使用仍是未决产品决定 |
| [局域网认证决定](../architecture/lan-authentication-decision.md) | 一个受信实验室局域网内的 `administrator` 与 `member`；没有项目 ACL，也没有租户隔离 |
| [运行恢复决定](../architecture/run-recovery-administrator-operations-decision.md) | 管理员按精确身份失败或重新排队；SQLite 仍是生命周期权威 |

PostgreSQL、对象存储、多机执行和第二套 Viewer 继续留在
[部署决定](../architecture/supported-deployment-runtime-decision.md) 的不支持清单里。
只有出现下面「容量」一节写明的实测瓶颈，或另有书面需求时，才单独立项。

## 安装包、版本与可恢复构建

### 已有证据

已发布且可下载的正式资产只有两套，均含 wheel、sdist 和只覆盖这两个文件的
`SHA256SUMS`：

| 发布 | 提交 | wheel SHA-256 | sdist SHA-256 |
| --- | --- | --- | --- |
| [v0.4.0](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/releases/tag/v0.4.0) | `ba29c256d782d60ec4ad3e2d8c2fadddc11aae8c` | `eb131621f9cd7acf542c2ef0228f3fcfbbc57f06bcec5e89f6b655a868aceb00` | `5a129aef0097969ffa1230e4a57e3f4cac2c960ccf60bdb0a2b7e87e1b5a12b1` |
| [v0.4.1](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/releases/tag/v0.4.1) | `813bccb75ba2b107eb598a61ecbad98016a498f8` | `0e005276f1dd97d87e27778412e078f73680f99310b977521d153ebef60296b6` | `715cca9d952eca5e8b816c8fa40a6494865ba4f3098cf802400e1fb7c1f73847` |

上表是 2026-11-19 从 GitHub Release 资产读到的 `SHA256SUMS`。本轮没有下载
wheel 或 sdist，也没有安装。v0.4.0 候选记录里的日期前校验和与上表不同，
那份记录不能代替已发布资产。

v0.4.1 发布说明写明该补丁不改变数据库 schema。v0.4.1 提交是当前基线的祖先。
基线仍把 `encode_pipeline.__version__` 和 API 版本写成 `0.4.1`，但迁移清单
`src/encode_pipeline/persistence/migration-execution-inventory-1.0.0.json`
的 head 已是 `20261119_18`，清单 SHA-256 为
`c7de12f1710199bac0c76ca82ec7b0fc2b61da16f304818006f4d9bed503a790`。
v0.4.1 之后才进入主干的修订是 `20261110_16`、`20261118_17` 和
`20261119_18`。因此已发布的 v0.4.1 安装包不是这份基线的可恢复构建。

可恢复构建沿用现有定义，不另造一种包：

- sdist 在检出之外能独立打出 wheel，并与直接安装的发布 wheel 到达同一版本和同一迁移 head；
- 同一 schema head 的两个平台 bundle 之间可以按部署决定回滚；
- 新 writer 启动之前，可以用迁移前的已校验 SQLite 备份恢复。新 writer 启动是
  既有决定中的 point of no return，之后不把备份恢复写成无损降级。

v0.4.0 候选证据记录了该候选自己的两次检出外安装，以及用已发布 v0.3.0 wheel
造出的 schema 08 记录迁到 `20260827_15` 后行仍保留、备份字节不变。仓库里没有
把同一次主机窗口写成已执行的回执。发布清单要求的 `status`、`doctor`、
`verify`、同 schema 升级回滚和 ENCODE 运行时切换，仍是待授权主机窗口的动作。

### 待测

这些动作留到下一次服务器候选冻结之后，按发布清单执行。本轮不跑。

1. 选定一个干净的精确提交。它的包版本、API 版本、前端版本和迁移清单 head
   必须是同一身份。当前基线版本串仍是 `0.4.1`、head 已是 `20261119_18`，
   所以它还不是该候选。
2. 按发布清单从该提交重新构建 wheel 与 sdist，再从新下载的 `SHA256SUMS`
   做 `sha256sum --check`。
3. 在两个全新的 Python 3.12 环境里分别安装发布 wheel，以及由 sdist 重新
   打出的 wheel。两次都做 `pip check`、版本与 CLI、打包前端和一份新 SQLite
   库，并确认新库停在该候选自己的 head。
4. 从该候选实际承认的上一发布 schema 做一次备份、迁移、point of no return
   之前的恢复，以及新 writer 启动后对旧 schema 的失败关闭。回滚只在声明了
   同一 head 的 bundle 之间证明。
5. 一次受权主机窗口：安装、`status` / `doctor` / `verify`、一次同 schema
   平台升级与回滚。发布清单已经排除全量故障注入、第二台主机、浸泡和压测。

### 验收条件

服务器或共享试用只安装满足上一节五项的候选。已发布的 v0.4.0 或 v0.4.1
资产可以继续服务它们自己的 schema 15 部署；把它们装上当前主干所在的机器，
不能当作当前目录、作者草稿或实验绑定已经随包交付。

## Peri 资格

### 已有证据

[资格记录](peri-runtime-qualification/README.md) 把受限运行时钉在 binary
SHA-256 `5cc82caf6e04578f19ff5790fbc4eb8c9e115793ff6950a699c9e2ad28c50f47`。
对应 Peri 提交 `c7d154c9985cbe08ba25e964cb7468d15014e6d0`，八文件补丁
SHA-256 `06272e1a66423a9e562509a38d044cd8506a1c10e61ef1fd3fe239107680587b`。
`review-decision.json` 的决定是
`accepted_for_fixed_restricted_A1_C4_runtime`。同一文件写明这不包括整个
A1、真实供应商、生产容量或长期 fork。

重新资格验证的触发条件已经写在该记录里：只有 Peri 提交、八文件补丁或隔离
配置改变时，才重跑 23 项 native 矩阵或 ACP 编译任务。把同一摘要的文件放到
另一台机器上，不构成新的矩阵。

当前基线的产品针仍是 `src/encode_pipeline/services/peri_runtime/identity.py`
里的同一摘要。操作者经 `HELIXWEAVE_PERI_BINARY` 提供文件；摘要不符则保持
不可用且不启动。安装包不携带 Peri binary。

[PR #232](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/232)
的 `93a7631d7387438898051949c15db2d90836191a` 在 2026-11-19 的主持复核中
通过了离线传输：合成字段、撤权后不再采用，以及建议来源。该提交不是当前
基线的祖先。它不使用供应商凭据，也不能作为云端模型已经联调的证据。真实
供应商仍是单独的后续事项。

### 待测

- 服务器候选若继续使用上述摘要：只核对操作者提供的文件摘要，并复用资格
  记录中的两条 pytest（摘要不符保持不可用；摘要相符时在真实 cgroup 和本地
  stub 上读两个工具并拒绝旧代次）。不重跑 23 项矩阵。
- 服务器候选若更换 Peri 提交、补丁或隔离配置：先按资格记录整项重做，再改
  产品针。未重做之前，旧摘要的接受结论停在旧身份上。
- 共享试用若要调用真实模型：先完成供应商、预算、外发字段和凭据边界的书面
  决定，再做联机验收。只提供目录和已有结果的服务器不把这项列为门槛。

### 验收条件

不调用真实模型的服务器，在 binary 摘要与产品针一致、不一致时失败关闭，
并且没有把离线传输写成联机供应商时，沿用现有资格。更换身份而未重做资格，
或凭 #232 宣称云端模型可用，都不满足 Peri 重新资格验证。

## WashU 分发

### 已有证据

[分发说明](wuepgg-62.0.2-distribution.md) 固定 npm `wuepgg@62.0.2` 及其
registry integrity，并比较过上游 `LICENSE.txt`。许可排除商业使用和销售；
副本须带协议和 Washington University in St. Louis 的版权声明。依赖
`wuepgg3-track@9.0.3` 没有自己的 `license` 字段。HelixWeave 仓库许可没有
因此改变。该说明明确把「某一次使用是否属于非商业」留作产品决定，本身不
构成华盛顿大学的授权。

当前基线的打包前端清单
`src/encode_pipeline/frontend_assets/asset-manifest.json` 含
`third-party-notices/wuepgg-62.0.2-LICENSE.html`，文件 SHA-256
`ba86b108e86063a3307752c5f6bfd4ae3fdfba59b69cff91b0870da8f7ec1dd2`。清单
身份是 `sha256-b3b911aeda4746c0a81a8c8e1606dc14404c34ff5a153cd88c80f35179d45197`，
`frontend_version` 仍是 `0.4.1`。这是检出内的打包结果，不是新的 Release
资产。源码树里的许可副本本身不到达已构建前端；到达路径是该说明里的
`frontend/public/` 与前端打包器。

### 待测

在把构建后的前端发给开发者以外的人之前，用该候选实际打出的包确认许可 HTML
仍在打包清单内，并且浏览器能以 `text/html` 打开它。轨道读取、外部代理和
独立 Viewer 构建仍由 Viewer 自己的交付验收，不在这次分发决定里重做。

### 验收条件

共享服务器上的前端只在书面决定确认该次使用落在现有非商业许可之内时分发。
决定未写下，或使用超出该许可时，停止分发并先取得单独授权。许可 HTML 必须
来自这次分发的同一个打包清单。源码树中有一份副本，或当前基线的
`0.4.1` 清单里有这个文件，都不能代替下一次候选包的检查。

## 备份、恢复与权限

### 已有证据

部署决定要求：停 writer、确认进程消失、做一次 root 拥有的不可变 SQLite
备份、再迁移、核对完整性与目标 head，然后才切换状态并启动。point of no
return 之前可以恢复该备份；之后自动恢复备份被禁止。`status`、`doctor`、
`verify` 的在线读检查不能代替这次离线备份。

运行恢复决定提供的是管理员失败和重新排队。它不恢复文件字节，也不制造
成功的科学结果。

[PR #233](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/233)
的 `41a4234dce17ac9f7df68cb38d311edd0a484eec` 在隔离的临时 SQLite 和文件
副本上演练了单份结果恢复。2026-11-19 的主持复核通过了该演练，并写明它不
包含整库备灾、长期盘迁移或服务器部署。该提交不是当前基线的祖先。长期盘
迁移仍是后续事项，现有恢复方法是放回原布局。

认证决定的权限边界是：同一受信局域网、两个角色、成员可读共享项目的运行
与产物，恢复和参考变更只属于管理员。明文 HTTP 不能保护 bearer cookie，
不是生产安全模式；离开本机回环的部署使用该决定里的 HTTPS `Secure` cookie。
决定同时写明这不产生租户隔离或按项目的 ACL。

### 待测

共享真实数据进入主机之前，用书面回答两件事，并只演练与回答相符的现有方法：

1. 访问范围。列出将使用这台主机的组。若任一组不得读取另一组的项目、运行
   或文件，现有两角色模型不够，先完成隔离再放入数据。若书面结论是一个
   实验室、不需要组间隔离，则沿用现有角色，不在部署窗口里加 ACL。
2. 恢复目标。写明要恢复的是单份已发布结果、SQLite 生命周期库，还是两者，
   以及 point of no return 之后接受的丢失。然后在一次性副本上演练这个目标：
   单份结果沿用 #233 的原布局核对；库备份沿用部署决定的停写、备份、校验和
   point of no return 前恢复。演练使用候选自己的 schema head。整库灾备和新
   长期盘不从这次演练升级为已有能力。

离开开发者本机的试用还要核对 HTTPS cookie 策略已经按认证决定启用。

### 验收条件

没有书面访问结论和恢复目标时，不把共享真实数据放上服务器。需要组间隔离
而隔离尚未落地时，同样停止。演练必须发生在一次性副本上，并且核对的是
演练前后同一结果身份或同一备份修订；生产库上的一次复制不能代替它。

## 容量

仓库里没有当前基线的服务器容量、并发或吞吐测量。部署决定和 v0.4.1 发布
说明都把生产规模吞吐留在证明范围之外。本地试用清单也排除性能基准。

容量数字只有在受权主机窗口里、针对那次候选测到之后才能写入。没有测量结果
不是瓶颈。出现下列任一情况时再单独立项，而不是提前改拓扑：

- 已写下的恢复目标在单机 SQLite 或现有文件存储上无法按部署决定完成；
- 一次受权试用里，现有 WashU 路径不能打开一份被点名的结果，并且需要另一套
  Viewer；
- 书面需求点名 PostgreSQL、对象存储、多机或第二 Viewer。

在此之前保持一台主机、SQLite、现有文件存储和当前 WashU 构建。

## 技能

现场试用若记录到技能与所装候选不一致，再改对应技能。当前三个技能是
[helixweave-run](../../.agents/skills/helixweave-run/SKILL.md)、
[helixweave-qc](../../.agents/skills/helixweave-qc/SKILL.md) 和
[helixweave-troubleshoot](../../.agents/skills/helixweave-troubleshoot/SKILL.md)。
本轮没有现场反馈，所以不改它们。

## 复核这份清单

在干净的上述基线上：

```bash
git rev-parse HEAD
git rev-parse 'v0.4.0^{commit}' 'v0.4.1^{commit}'
git merge-base --is-ancestor 93a7631d7387438898051949c15db2d90836191a HEAD
git merge-base --is-ancestor 41a4234dce17ac9f7df68cb38d311edd0a484eec HEAD
python3 -c 'import pathlib,re; text=pathlib.Path("src/encode_pipeline/__init__.py").read_text(); print(re.search(r"__version__ = \"([^\"]+)\"", text).group(1))'
```

两条 `merge-base --is-ancestor` 在该基线上的退出码都是 1。版本命令读取
`src/encode_pipeline/__init__.py`，打印 `0.4.1`，不要求已安装的包。迁移
head 以清单文件的 `heads` 为准，当前为 `20261119_18`。发布校验和用对应 Release 的 `SHA256SUMS` 复核；本轮没有
把 wheel 装进环境。

文档链接用 `test/docs/test_internal_links.py` 检查。该测试不安装产品、
不迁移数据库，也不接触服务器。

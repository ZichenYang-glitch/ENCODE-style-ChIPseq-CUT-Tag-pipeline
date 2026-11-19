# Peri 工作台的运行确认

IRENE-45 实施 [#245 的仅创建与取消决定](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/blob/e4dc005e4bae62e5097f92f400971e54a33eba9b/docs/architecture/agent-run-operation-boundaries.md#实施决定仅向-peri-开放创建与取消)。

## 入口与服务

成员在配置草稿中连接 Peri，发出结构化的运行请求。创建使用该草稿当前修订校验所得的 snapshot 和成员填写的 tags。取消使用成员指定的同工作流 run ID。平台预览返回动作和具体目标，成员另点确认。更换动作、snapshot、tags 或 run ID 需要新的预览和确认。

现有作者会话的建议请求只接收已声明的建议和解释意图，现有人用运行入口接受成员直接操作。工作台通过作者会话下的两条浏览器内部 POST 路由完成预览和确认。这两条路由复用成员与 CSRF 校验，沿用已就绪的 Peri 连接记录和当前 owner epoch，并绑定登录、会话、连接代数、动作和目标。确认由服务端用宿主持有的连接材料签发，浏览器和模型均无法取得签名密钥。连接关闭、重新接入或登录撤销后，旧确认失效。同一确认可在签发后十分钟内用原目标重试。

确认复用现有创建和取消的生命周期服务及安全审计。创建返回 run ID 和首次创建或重放的 HTTP 状态。取消 202 表示已请求停止。工作台用现有 GET 回读该 ID，并跟踪停止请求直到权威终态。浏览器内部路由沿用作者会话的认证和错误边界，公开 OpenAPI 保持原契约。建议工具、模型传输、现有人用按钮、预检和启动沿用现有路径。

## 验收

用真实 FastAPI 应用、隔离 SQLite、合成 Peri 连接和记录调用的队列验证。不派发科学工作流。创建首次返回 201、原目标重试返回 200 和同一 ID，执行分配与入队计数为零。取消 A 时 B 不变，running 的 202 经 GET 仍为 running，直到写入权威终态。核对登录撤销、重新连接、错误会话、换目标、改 tags 和普通 Agent 请求。浏览器验收从工作台发出请求，检查显示目标和独立确认按钮，核对确认后读取的确切 run ID。

在仓库根目录用项目的 `ci-fast` Python 环境和 `.nvmrc` 指定的 Node 版本执行。

```bash
python3 -I -S scripts/checkout_bootstrap.py --repository-root . pytest test/api/test_peri_run_actions.py
npm --prefix frontend ci
npm --prefix frontend exec -- playwright install chromium
npm --prefix frontend run build
python3 -I -S scripts/checkout_bootstrap.py --repository-root . pytest test/api/peri_run_actions_browser.py
```

Chromium 用例发出真实 HTTP 请求，丢弃第一次创建响应后重试，并等待确切 run 的 GET 返回停止终态。测试使用的 Peri 连接和队列确认均为本地合成。

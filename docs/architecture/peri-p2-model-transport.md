# P2 模型传输契约（离线）

契约版本：`p2-transport/1`。状态：IRENE-25 的第一份本地交付。

本文件落实 [目录、Agent 与 Viewer 边界](catalog-agent-viewer-boundaries.md) 第 5 节的 P2 要求：把资格 stub、测试 transport、未来真实供应商分成三个角色。本轮没有选定供应商、外发字段或额度，不发起真实模型请求，也不读取或传递凭据和实验数据。

实现在 `src/encode_pipeline/services/peri_runtime/model_transport.py`。资格探针仍是 `stub.py` 里的 `ModelStub`。本契约不改 `services/authoring.py`、公开 API 或目录迁移。

## 1. 三种角色

| 角色 | 代码 | 现在做什么 |
| --- | --- | --- |
| 资格 stub | `open_qualification_model_transport` → `ModelStub` | 在本机私有 socket 上配合固定 Peri binary 做资格探针。`PeriAuthoringRuntime` 继续把这个开启函数传给 `NativeSession`。 |
| 测试 transport | `SyntheticSuggestionTransport` | 进程内交换合成字段。响应决定建议内容；只接受允许字段；超时、取消或撤权后不投递、不采用。 |
| 真实供应商 | `open_vendor_model_transport` | 直接拒绝，reason 为 `vendor_unselected`。参数被忽略。不读环境变量，不打开 socket，不构造客户端。 |

三者的角色名分别是 `qualification_stub`、`test_transport`、`vendor_transport`。测试 transport 没有 `ModelStub` 的 `script` / socket 行为，资格 stub 也没有 `deliver` / `adopt`。供应商入口不返回其中任何一个。

资格 stub 的通过不能写成生产模型可用，也不能写成这份建议契约已联机验收。

## 2. 测试 transport 的一次交换

调用方提供：

- `call_id`，形如 `call_` 加 8–32 位小写字母或数字。
- `allowed`：本次允许的合成 token，1–8 个，互不重复。
- `projected`：实际交给 responder 的合成字段，必须是 `allowed` 的子集。
- `authorize()`、`clock()`、`timeout_seconds`，以及可选的 `cancelled()`。

Token 形如 `syn.alpha`。值是 1–64 位的 `A-Za-z0-9_.:-`。这是离线夹具文法，不是目录字段名，也不是供应商 schema。8 个字段是夹具形状，不是花费额度。

`deliver` 按这个顺序处理：

1. 请求不合法则拒绝，reason `malformed_request`，不调用 responder。投影里出现允许集之外的 token 时 reason 为 `field_not_allowed`，同样不调用 responder。
2. 调用 responder 之前检查授权和取消。已撤权则 `authorization_revoked`；授权回调抛错则 `authorization_failed`；已取消则 `cancelled`。这三种情况都不调用 responder。
3. responder 只收到 `TransportProjection`：`call_id` 和 `projected`。它看不到允许集里未投影的字段。
4. responder 返回后再次读时钟。耗时大于 `timeout_seconds` 则丢弃返回值，reason `timeout`。等于时限仍算及时。时钟倒退或非有限值则 `malformed_clock`，同样丢弃返回值。
5. 返回之后再次检查授权和取消。失败则丢弃返回值。
6. 响应至少包含一个字段，token 必须全部属于 `allowed`，值必须符合夹具文法。任一 token 越界则整份拒绝，reason `field_not_allowed`，不保留其余字段。
7. 通过后字段暂存在这个 transport 实例里。返回的 `Delivery` 只有 `delivery_id`、`call_id`、`response_id` 和 `source`，不含字段值。

`adopt` 再次检查授权和取消。通过则返回 `Suggestion` 并删掉实例内的暂存；失败则先删掉暂存再拒绝。同一次投递只能采用一次，再次采用的 reason 是 `not_delivered`。撤权后即使授权恢复，也不能用旧 `Delivery` 再采用，需要一次新的 `deliver`。

建议字段等于被接受的响应字段，顺序也保留。因此同一请求配上两份不同响应，会得到两份不同建议。建议不复制 `projected` 里的值，除非响应自己写了同样的值。`source` 固定为 `test_transport`。契约版本写在建议上。

`suggestion_changes` 把已采用建议变成 `{"field_token", "value"}` 元组，供以后的 authoring 所有者接入。它不导入 authoring 服务。手写的、带路径或其他非法值的建议会被拒绝，reason `malformed_suggestion`，异常文本只有这个 reason。

用量记录只有 `attempted`、`responder_calls`、`delivered`、`adopted` 四个计数。它不是供应商配额，也不授权任何花费。

拒绝对象和 `TransportRejected` 的文本只有 reason code。回调里的异常消息不会进入 reason。

## 3. 离线能证明的事

用合成字段和注入的时钟、授权回调，可以重复证明：

- 三种角色是不同的对象；供应商入口在设了假的 endpoint 或 key 时仍立即拒绝。
- 更换响应值，已采用建议的 `field_token` / `value` 跟着变，并且不同于请求投影里的值。
- responder 看见的投影只含调用方给出的允许字段；越界 token 不会被投递或采用。
- 时钟越过本次时限后，即使 responder 已经返回，也不投递、不能采用。
- 撤权发生在调用前、响应返回后、或 `adopt` 之前时，相应步骤停止；`adopt` 前撤权会丢掉已暂存的字段。
- 取消与撤权使用同一道门，reason 为 `cancelled`。
- 拒绝文本不包含被拒的合成值。

复现：

```bash
./.local/envs/ci-fast/bin/python -I -S scripts/checkout_bootstrap.py \
  --repository-root . \
  pytest test/services/test_peri_model_transport.py
```

解释器须是指向本检出的可编辑安装，创建方式见
[本地平台运行时](../development/local-platform-runtime.md)。这条命令不启动 Peri，
不访问网络，也不读取模型凭据。

## 4. 仍要联机验收的事

下列事项没有做，也不能用本文件或测试 transport 的结果代替：

- 选定一家供应商，并按它的官方请求与响应契约实现客户端。
- 由 operator 配置 endpoint 和凭据，并证明它们不来自用户输入或模型指令。
- 冻结真实外发字段和额度；记录真实 token 用量。
- 对真实连接验证超时、取消、断连和撤权，包括晚到响应不会被采用。
- 把建议写入目录提案，并走平台的采用事务、回执和撤权顺序。
- 用资格 stub 或本测试 transport 的计数宣称生产模型可用。

## 5. 共享文件与下一步

本轮不改这些共享面：

- `src/encode_pipeline/services/authoring.py`，以及 persistence / API 里的 authoring 路由。P1b 正在使用公共 authoring 服务。
- 公开 HTTP 契约、OpenAPI 和目录迁移。
- `NativeSession` 的 `model_transport` 参数。它的类型仍是 `Callable[[Path], ModelStub]`，只打开资格 stub。

下一步要等 authoring 所有者接入，而不是在这里重做一套采用流程。接线时新增一个窄调用即可：

```python
outcome = transport.deliver(request)
if isinstance(outcome, Refusal):
    raise ...
suggestion = transport.adopt(outcome)
changes = suggestion_changes(suggestion)
```

`changes` 的每一项只有 `field_token` 和 `value`。authoring 所有者决定是否把它们送进现有提案校验。不要把 `SyntheticSuggestionTransport` 或未来的供应商客户端塞进 `model_transport=`：那个参数是资格探针的 socket 开启函数。真实供应商传输需要自己的 operator 配置对象，并单独做第 4 节的联机验收。

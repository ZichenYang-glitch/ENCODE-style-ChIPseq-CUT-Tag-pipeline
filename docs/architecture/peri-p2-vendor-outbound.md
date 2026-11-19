# P2 真实供应商外发（设计）

状态：IRENE-31 的设计契约。契约版本沿用离线传输的 `p2-transport/1`。
本文件不新增传输实现，不打开供应商连接，不读取凭据或实验数据。

离线传输已经由 IRENE-25 交付并经 Tibo 复核，提交
`93a7631d7387438898051949c15db2d90836191a`，Draft
[PR #232](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/pull/232)。
该提交不在本分支上。专题文档是：

- [目录、Agent 与 Viewer 边界](catalog-agent-viewer-boundaries.md) 第 5 节。
- [离线传输契约](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/blob/93a7631d7387438898051949c15db2d90836191a/docs/architecture/peri-p2-model-transport.md)
  与 [model_transport.py](https://github.com/ZichenYang-glitch/ENCODE-style-ChIPseq-CUT-Tag-pipeline/blob/93a7631d7387438898051949c15db2d90836191a/src/encode_pipeline/services/peri_runtime/model_transport.py)。
- [集成 ADR 的数据边界](peri-agent-integration-decision.md) 与
  [契约附录的最小外发](peri-agent-integration-contracts.md)。

验收输入在 [peri-p2-vendor-outbound-cases.json](peri-p2-vendor-outbound-cases.json)。
`success`、`rejection`、`timeout` 使用离线 `SyntheticSuggestionTransport` 的文法。
成功例在该 stub 上可以交付并采用；拒绝例与超时例在调用结果里只留下 reason code。
这三例都停在本地 stub，不经过供应商入口。`owner_grant_exceeded`、
`budget_exhausted`、`budget_unproven` 和 `timeout_consumes_reservation`
是发送前判定的合成输入，不交给该 stub 的 `deliver`。

## 1. 留给负责人的决定

| 决定 | 本文件的值 | 未填写时的拒绝码 |
| --- | --- | --- |
| 首家供应商、官方 API、模型名、endpoint | 空 | `vendor_unselected` |
| 允许外发的字段 | 空。夹具 token `syn.*` 只说明形状 | `outbound_fields_unset` |
| 预算：币种、金额的最小货币单位、请求数、输入/输出 token 上限 | 空 | `budget_unset` |
| 凭据边界：由谁配置、存放位置、能否进入 Peri 进程 | 空 | `credential_boundary_unset` |

四个值有任一为空，未来的供应商客户端在打开连接之前拒绝，reason 只有上表的代码，
按表中的顺序命中即停。已实现的 `open_vendor_model_transport` 仍对任何参数抛出
`vendor_unselected`：它忽略参数，不读环境变量，不打开 socket，不构造客户端。
后三个拒绝码尚未实现。在它们落地之前，空决定的可观察行为仍是不发送请求。

本地 stub 的成功例不填这些决定。它证明的是投影、授权期和时限，不构成外发许可。

## 2. 供应商无关的请求字段

平台交给传输层的对象沿用离线 `TransportRequest`。它是供应商无关的请求，各家 HTTP
正文留到负责人选定官方契约后再写。

| 字段 | 离线契约 | 外发时 |
| --- | --- | --- |
| `call_id` | `call_` 加 8–32 位小写字母或数字 | 作为不透明关联号发出。它不是凭据 |
| `allowed` | 1–8 个互异 token | 调用方声称的候选集，留在平台。它不是许可 |
| `projected` | `allowed` 的子集。token 形如 `syn.alpha`，值是 1–64 位 `A-Za-z0-9_.:-` | 唯一允许出现的业务载荷 |
| 授权期上的外发清单 | 签发该 `grant_...` 时冻住的 token 集合 | 外发许可。请求的 `allowed` 不能替换它 |
| 授权期 | `authorize()` 返回 `grant_` 加 8–32 位小写字母或数字 | 只记在本地来源记录。不发给供应商 |
| `timeout_seconds` | 有限且大于 0 的秒数 | 平台自己的截止时间。在官方契约选定前，不把它写成供应商参数 |
| `cancelled()` | 可选 | 平台本地门。不发给供应商 |
| `contract_version` | `p2-transport/1` | 记在建议和来源记录上。不写进提示词 |

外发体只有 `call_id` 与 `projected`。样本名、路径、序列、供体标识、自由文本、
未投影字段、工具目录、完整 schema、endpoint 和凭据都留在平台。

资格 stub 是本分支 `stub.py` 的 `ModelStub`，与建议传输是两个角色：

| | 资格 `ModelStub` | 建议 `SyntheticSuggestionTransport` |
| --- | --- | --- |
| 位置 | 私有 Unix socket，`POST /v1/chat/completions` | 进程内 `deliver` / `adopt` |
| 凭据 | 固定 `Bearer synthetic-a1-key` | 不持有钥匙 |
| 体积 | 正文 1–262144 字节，最多 256 次请求 | 1–8 个字段 |
| 目录 | `tools` 必须正好是 `SearchExtraTools` 与 `ExecuteExtraTool` | responder 只看见 `TransportProjection` |
| 撤权 | `authorize()` 为假时记录 `authorization_revoked`，业务标记不进入已记录请求 | 授权期在调用前、返回后、采用前三次核对 |
| 超时 | 连接 8 秒；没有脚本响应时等待 8 秒 | 调用耗时大于 `timeout_seconds` 则丢弃响应。采用截止是交付时钟加上同一时限 |

资格 stub 的 200 响应不表示建议已被采用，也不表示供应商可用。

## 3. 授权与撤权时点

一次交换按这个顺序做决定。撤权看授权期字符串，布尔“当前允许”不能代替某一期。

1. 请求形状不合法则 `malformed_request`。投影里有本次请求 `allowed` 之外的 token 则 `field_not_allowed`。这两步只核对请求自己的形状，不读授权，不调用 responder，也不把 `allowed` 当成负责人的许可。
2. 四项负责人决定有空项则按第 1 节拒绝，不打开连接。本文件仍把这四项留空，所以供应商客户端仍然不发送。
3. 外发字段以冻在当前授权期上的负责人清单为准。清单在签发该 `grant_...` 时固定；换一期就要换一份清单。请求里的 `allowed` 只是调用方声称的候选集，不能代替这份清单，也不能把别的期的清单借来用。`allowed` 和 `projected` 中的每个 token 都必须属于当前期的清单。任一超出则 `outbound_field_not_granted`，供应商调用数为 0，不打开连接。验收例 `owner_grant_exceeded`：期 `grant_aaaaaaaa` 的清单只有 `syn.alpha`，请求的 `allowed` 含有 `syn.alpha` 与 `syn.beta`，并投影 `syn.beta`。这例的清单是合成验收输入，不填写第 1 节的外发字段决定。
4. 发送前做第 5 节的余额证明。`budget_exhausted` 与 `budget_unproven` 都不打开连接，供应商调用数为 0。
5. 调用前读取授权期。`None` 为 `authorization_revoked`；回调抛错或期字符串不合文法为 `authorization_failed`；已取消为 `cancelled`。三种情况都不调用 responder。
6. 调用返回后读时钟。耗时大于 `timeout_seconds` 则丢弃返回值，reason `timeout`。等于时限仍算及时。时钟倒退或非有限值为 `malformed_clock`，同样丢弃返回值。已经发出后才超时的调用仍算一次发送。
7. 返回后的授权期必须仍是进入调用时的那一期，取消也必须仍未发生。撤权、换成别的期或取消都丢弃返回值，不保留字段。
8. 响应字段必须全部属于 `allowed`，并且符合夹具文法。任一 token 越界则整份 `field_not_allowed`。
9. `adopt` 核对这次投递自己的截止时间。晚于截止时间为 `timeout`，并删除暂存。等于截止时间仍可采用。当前授权期必须等于投递记下的期：撤权是 `authorization_revoked`，换成新的期是 `authorization_period_mismatch`。新的期只采用它自己的新 `deliver`。同一次投递只采用一次，再次采用为 `not_delivered`。

`attempted` 在 `deliver` 进入时加一。只有真正调用了 responder，`responder_calls` 才加一。超时发生在调用之后，所以超时例的 `responder_calls` 为 1，`delivered` 仍为 0。

## 4. 脱敏样例

外发投影是允许字段的副本。下面的候选对象只让 `syn.alpha=from_request` 进入外发体。
`DonorA-liver`、`/data/private/sample.fastq`、`sk-live-example`、`/etc/secret` 和
“ignore previous rules”留在平台。`sk-live-example` 是夹具标记，不放入环境变量，也不发给供应商。

现有 `RedactionPolicy` 是补防线，正则结果如下，均由本分支代码得出：

| 输入 | 输出 |
| --- | --- |
| `see /data/private/sample.fastq` | `see <FILE_PATH>` |
| `api_key=sk-live-example`，另有 `path=/data/private/sample.fastq` 与 `code=SAMPLE_SELECTION_REQUIRED` | 钥匙变为 `<REDACTED>`，路径变为 `<FILE_PATH>`，`code` 保留 |
| 键名为 `token`、值为 `syn.alpha` | 值变为 `<REDACTED>`。敏感键名单包含 `token` |
| 自由文本 `token sk-live-example stays` | 原文保留。自由文本里的钥匙形状不会被这条正则去掉 |

因此不能把建议投影先交给 `RedactionPolicy` 再外发：键名 `token` 会被当成秘密擦掉，
而自由文本里的标记又会留下来。外发门是第 3 节冻在授权期上的负责人清单。
拒绝对象和 `TransportRejected` 的文本只有 reason code。`hidden_value` 与
`late_value` 不进入拒绝文本。

## 5. 费用上限

费用记录的形状是：`status`、`currency`、`max_amount_minor`、`max_requests`、
`max_input_tokens`、`max_output_tokens`。当前 `status` 为 `owner_unset`，其余五项为空。
本文件不填写金额、币种或 token 单价。

离线 `UsageRecord` 只有 `attempted`、`responder_calls`、`delivered`、`adopted`。
这四项是进程内计数，不是配额，也不授权花费。预算为空则在第 3 节第 2 步得到
`budget_unset`，responder 不被调用。

上限填写之后仍然只是封顶，不是余额。打开连接前要有一份整数证明：本次预留的
请求数、输入 token、输出 token 和金额最小单位，每一项都小于或等于账上剩余。
任一剩余小于本次预留，reason 为 `budget_exhausted`。预留或剩余缺字段、不是非负整数，
或账本读不到，reason 为 `budget_unproven`。两种情况都不发请求，供应商调用数为 0。
金额比较用整数最小货币单位。

已经发出、随后超时的调用消耗本次预留。这笔预留不退回，这次超时也不能当成免费重试。
再发送一次，要按扣减后的余额重新证明。`timeout_consumes_reservation` 里，前一次发送
已经把四项剩余降到 0，下一份同样的预留得到 `budget_exhausted`，供应商调用数为 0。
stub 上的 `timeout` 例没有打开供应商连接，也不产生余额；它只说明晚到响应被丢弃。

## 6. 响应来源记录

成功例采用之后，来源记录包含：

| 字段 | 成功例的值 |
| --- | --- |
| `contract_version` | `p2-transport/1` |
| `source` | `test_transport` |
| `call_id` | `call_stub0001` |
| `response_id` | `resp_stub0001` |
| `grant_period` | `grant_aaaaaaaa` |
| `delivered_at` | `1.0` |
| `adopt_deadline` | `2.0`（交付时钟 `1.0` 加上时限 `1`） |
| `field_tokens` | `syn.alpha` |
| `vendor` | 空 |
| `snapshot_id` | 空 |
| `cost_cap_status` | `owner_unset` |

建议字段等于被接受的响应字段，本例是 `syn.alpha=from_response`。请求投影里的
`from_request` 不复制进建议，除非响应自己写了同样的值。

离线 stub 今天返回的 `Suggestion` 含有契约版本、`call_id`、`response_id`、字段和
`source`。授权期、交付时钟、截止时间、供应商、snapshot 和费用状态是这份设计要求将来
随建议保存的来源记录，当前 stub 对象尚未带上这些属性。来源记录不包含提示词、原始响应正文、
路径或凭据。

用户采用并生成新 snapshot 仍走现有 authoring 与 snapshot 路径。本设计不写入 snapshot，
所以 `snapshot_id` 保持为空。真实响应、用户采用和新 snapshot 要等第 1 节四项决定填写之后
再验收。

## 7. 本地 stub 的三项验收

夹具时钟和授权期按第 3 节的读取次数列出。成功例读取三次授权期和三次时钟：调用前、
返回后、采用时。拒绝例在投影检查处停止，授权与时钟读取次数为 0。超时例读取一次授权期
和两次时钟，在返回后的时钟比较处停止。

| 例 | 输入要点 | 预期 |
| --- | --- | --- |
| `success` | 时限 1。时钟 `0.0`、`1.0`、`1.0`。授权期始终 `grant_aaaaaaaa`。投影 `syn.alpha=from_request`。响应 `resp_stub0001` / `syn.alpha=from_response` | `deliver` 得到 `source=test_transport` 的投递。`adopt` 的字段是 `from_response`。用量 `1/1/1/1`。来源记录见第 6 节 |
| `rejection` | 允许集只有 `syn.alpha`，投影是 `syn.hidden=hidden_value` | reason `field_not_allowed`。responder 调用 0 次。用量 `1/0/0/0`。拒绝文本不含 `hidden_value` 或 `syn.hidden` |
| `timeout` | 时限 1。时钟 `0.0` 然后 `5.0`。响应值 `late_value` | reason `timeout`。responder 调用 1 次，交付数与采用数均为 0。拒绝文本不含 `late_value` |

发送前四例不进入上面的 stub。第 1 节的四项决定在这些例子里仍然为空；例子自带的清单和余额只用于判定。

| 例 | 输入要点 | 预期 |
| --- | --- | --- |
| `owner_grant_exceeded` | 期 `grant_aaaaaaaa` 的清单只有 `syn.alpha`。`allowed` 含 `syn.alpha` 与 `syn.beta`，并投影 `syn.beta` | `outbound_field_not_granted`。供应商调用数 0，不打开连接 |
| `budget_exhausted` | 清单允许 `syn.alpha`。四项剩余都是 0，本次预留各为 1 | `budget_exhausted`。供应商调用数 0，不打开连接 |
| `budget_unproven` | 同一清单。费用证明不可建立 | `budget_unproven`。供应商调用数 0，不打开连接 |
| `timeout_consumes_reservation` | 前一次已经发出并超时，四项剩余都是 0。下一次仍预留各 1 | `budget_exhausted`。供应商调用数 0。前一次预留不退，不能免费重试 |

回放离线 stub 时只使用 `success`、`rejection`、`timeout`，每个例子新建一个
`SyntheticSuggestionTransport`。本分支没有这个模块。资格 `ModelStub` 不能用来宣称
这三例已经联机通过。发送前的四例由本地判定核对，不调用该 stub 的 `deliver`。

本单将来的真实外发只验一次另行获准的调用：第 1 节四项决定已经填写，并且第 5 节的
余额证明通过。那一次使用成功例的非敏感投影，得到真实响应，用户采用后生成新 snapshot。
拒绝和超时留在 stub。主动在真实连接上制造超时不在本单范围内，除非另有授权和预算。

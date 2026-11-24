# Use the draft page with a bound Peri model

This is an optional local owner deployment. The API and browser do not get a provider key. The qualified ACP process still performs admission offline. For real suggestions the owner runs a separate one-call worker built from Peri’s pinned **peri-model** crate. It loads the specified Peri settings file in memory and executes the active profile’s model adapter; it does not launch a personal agent with tools, plugins or conversation history. Neither the qualification stub nor the offline test transport generates the result in this mode.

The first supported case is **Suggest values → one integer field → one integer preference**. The model receives a fixed instruction, the approved field token and that integer only. An unsupported request, binding protocol, missing grant, malformed response, expiry or inadequate budget fails closed. There is no vendor fallback. Manual draft editing remains available. Run creation and cancellation are unchanged.

## Install and enable on the owner

Use the pinned Peri checkout, its Rust toolchain and populated offline cargo cache. Build in an empty directory outside the source tree. Cargo uses its offline cache; on a fresh build the pinned upstream wrapper also downloads and hash-checks three crate archives before patching them:

```sh
sh integrations/peri-bound/build.sh "$PERI_SOURCE" "$BUILD_DIR"
export HELIXWEAVE_PERI_MODEL_WORKER="$BUILD_DIR/target/debug/examples/helixweave-bound"
export HELIXWEAVE_PERI_MODEL_SETTINGS="$PERI_SETTINGS"
export HELIXWEAVE_PERI_MODEL_LEDGER="$OWNER_STATE/bound-model.sqlite"
python -m encode_pipeline.services.peri_runtime.model_grant init
python -m encode_pipeline.services.peri_runtime.model_grant binding
```

Use absolute paths. Protect the owner directory, executable, ledger and settings from other accounts. The ledger is independent of the platform database; init refuses to replace it. Keep it on local persistent disk and do not restore older balances to active grants. Set these three variables **on the owner process** and restart it using the site's existing owner service. No binding is configured in the web frontend. No variables means the original qualification setup; a partial configuration fails. With live mode selected, an unbound or unauthorized turn never silently returns a scripted suggestion.

Supported bindings currently have explicit `anthropic` protocol, HTTPS base URL, active alias, provider, model and API key in Peri settings (including a DeepSeek Anthropic-compatible binding). Other protocols, credential resolvers and implicit model defaults are not supported. The worker checks a binding fingerprint including the settings file revision. Editing or replacing settings invalidates old grants. The protocol adapter has no network retries, redirects or proxy discovery. HTTP loopback requires an explicit worker test flag used only by the fixture.

## Authorize one connected draft

Connect the assistant on the draft page. Copy its authoring session ID and the non-secret fingerprint from `binding`. Supply the authorized caps and tariff **for this binding**. The example has no default amounts:

```sh
python -m encode_pipeline.services.peri_runtime.model_grant grant \
  --platform-db "$PLATFORM_SQLITE" --sid "$SID" --binding "$BINDING" \
  --field config.threads --expires "$UNIX_EXPIRY" --requests "$MAX_REQUESTS" \
  --input-tokens "$TOTAL_INPUT" --output-tokens "$TOTAL_OUTPUT" \
  --output-per-call "$MAX_OUTPUT" --currency "$CURRENCY" --minor-units "$MAX_COST" \
  --input-rate "$MAX_INPUT_RATE" --output-rate "$MAX_OUTPUT_RATE"
```

Use one currency throughout a grant. Rates are **upper bounds in billing-currency minor units per million tokens**, including any applicable surcharges. Round rates up. Do not authorize an account/model whose worst-case prices or framing are unknown. This slice supports token-priced bindings only: the fixed single-message body is bounded to 1024 bytes; the operator must attest that 4096 input tokens covers that exact binding’s maximum billed input, including framing. Output is capped at the provider. No tools, history, sample IDs, storage paths or rate schedules are sent. The maximum cost reserved per call is the ceiling of `(4096 × input_rate + output_per_call × output_rate) / 1,000,000`.

The CLI reads the current instance, generation and owner epoch from the connected draft. The field and binding are immutable for that instance. An immediate SQLite transaction reserves one request, 4096 input tokens, the full output allowance and maximum cost before the call. Concurrent requests cannot share the last reservation. Failed, timed-out or abandoned attempts do not refund it. `reservations` records turn, requested model, reserved units and delivery; it does not persist prompts, results or credentials. The proposal is validated by the usual schema and draft checks. Adoption rechecks the grant and binding; the delivery window is at most 35 seconds, and never extends the grant or existing proposal expiry.

To revoke, use the instance ID printed by `grant`:

```sh
python -m encode_pipeline.services.peri_runtime.model_grant revoke --instance "$INSTANCE"
```

A new connection/owner epoch needs a new explicit grant. For rollback revoke the grant, remove all three owner variables and restart the owner. Keep the ledger as evidence.

## Verify without a paid request

The fast regressions cover grant creation and revocation, field and budget refusals, request counts and proposal adoption:

```sh
python -I -S scripts/checkout_bootstrap.py --repository-root . pytest -q test/services/test_bound_peri_model.py test/services/test_model_grant.py
```

For the full local test use a built wheel frontend, the repo’s Node/Playwright/Chromium setup, and set `HW_TEST_PERI_MODEL_WORKER` and `HW_TEST_PERI_NATIVE` to the absolute worker and admitted ACP binary. Run:

```sh
python -I -S scripts/checkout_bootstrap.py --repository-root . pytest -qs test/api/peri_bound_model_browser.py
```

The loopback provider receives one real protocol request via Peri’s adapter, returns 7 for a preference of 3, and Chromium checks suggestion and adoption. A second request with the depleted grant fails without reaching the provider. A separate wire test verifies that redirects, rate limits and a changed binding cannot trigger a second transmission. This verifies the page path, not a commercial provider’s service. A paid smoke test requires a separate explicit field grant and spending ceiling; use a synthetic draft.

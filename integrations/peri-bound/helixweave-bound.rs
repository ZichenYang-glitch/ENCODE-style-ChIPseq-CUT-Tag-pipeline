// Runs on the Peri side. Only this process opens the explicitly bound settings.
// stdin/stdout is a one-request protocol; never report provider error bodies.
use peri_model::{
    AnthropicConfig, AnthropicModel, Model, ModelMessage, ModelRequest, ModelRuntimeConfig,
    RetryConfig,
};
use serde_json::{json, Value};
use std::{
    env, fs,
    io::{self, Read},
    os::unix::fs::MetadataExt,
};
use tokio_util::sync::CancellationToken;
use url::Url;

fn token(v: &Value) -> Option<&str> {
    v.as_str().filter(|s| {
        !s.is_empty()
            && s.len() <= 128
            && s.bytes()
                .all(|b| b.is_ascii_alphanumeric() || b"_-.:/".contains(&b))
    })
}

async fn run() -> Option<Value> {
    let path = env::args().nth(1)?;
    // Freeze the non-secret revision from the same open descriptor as the binding.
    // Even a same-model key rotation invalidates the current grant.
    let mut file = fs::File::open(path).ok()?;
    let revision =
        |m: fs::Metadata| format!("{}:{}:{}:{}", m.dev(), m.ino(), m.ctime(), m.ctime_nsec());
    let version = revision(file.metadata().ok()?);
    let mut bytes = Vec::new();
    (&mut file).take(1_048_577).read_to_end(&mut bytes).ok()?;
    if bytes.len() > 1_048_576 || revision(file.metadata().ok()?) != version {
        return None;
    }
    let data: Value = serde_json::from_slice(&bytes).ok()?;
    let c = &data["config"];
    let alias = token(&c["active_alias"])?;
    let profile = &c["profiles"][alias];
    let id = token(&profile["provider"])?;
    let mut providers = c["providers"].as_array()?.iter().filter(|p| p["id"] == id);
    let p = providers.next()?;
    if providers.next().is_some() || p["type"] != "anthropic" || !p["auth"].is_null() {
        return None;
    }
    let model = if profile["model"].is_null() {
        token(&p["models"][alias])?
    } else {
        token(&profile["model"])?
    };
    let endpoint = Url::parse(p["baseUrl"].as_str()?).ok()?;
    let local = env::args().nth(2).as_deref() == Some("--loopback-test")
        && endpoint.host_str() == Some("127.0.0.1")
        && endpoint.scheme() == "http";
    if !(local || endpoint.scheme() == "https")
        || endpoint.host_str().is_none()
        || endpoint.password().is_some()
        || !endpoint.username().is_empty()
        || endpoint.query().is_some()
        || endpoint.fragment().is_some()
    {
        return None;
    }
    let key = p["apiKey"].as_str()?.trim();
    if key.is_empty() {
        return None;
    }
    let binding = json!({"provider":id,"alias":alias,"model":model,"endpoint":endpoint.as_str(),"protocol":"anthropic","revision":version});
    let mut raw = String::new();
    io::stdin().take(4097).read_to_string(&mut raw).ok()?;
    if raw.len() > 4096 {
        return None;
    }
    let r: Value = serde_json::from_str(&raw).ok()?;
    if r == json!({"op":"binding"}) {
        return Some(binding);
    }
    if r["op"] != "suggest" || r["binding"] != binding {
        return None;
    }
    let field = token(&r["field"])?;
    let value = i32::try_from(r["value"].as_i64()?).ok()?;
    let limit = u32::try_from(r["output_tokens"].as_u64()?).ok()?;
    if !(1..=1024).contains(&limit) {
        return None;
    }
    let prompt=format!("Suggest one integer for this field. Respect the preference if reasonable. Return only the integer. {}",json!({"field":field,"preference":value}));
    let request = ModelRequest::new(vec![ModelMessage::user_text(prompt)]).with_max_tokens(limit);
    let config = AnthropicConfig::new(endpoint, key, model)
        .without_cache()
        .with_max_tokens(limit)
        .with_runtime(
            ModelRuntimeConfig::default().with_retry(RetryConfig::default().with_max_attempts(1)),
        );
    let adapter = AnthropicModel::single_request(config).ok()?;
    let prepared = adapter.prepare_stream(request).ok()?;
    let checkpoint = prepared.checkpoint();
    // Bound the entire wire body, not just the user's preference.
    if serde_json::to_vec(&checkpoint["body"]).ok()?.len() > 1024 {
        return None;
    }
    let url = checkpoint["endpoint"].as_str()?.to_owned();
    use futures::StreamExt;
    let mut stream = prepared.start(CancellationToken::new()).ok()?;
    let mut text = String::new();
    while let Some(event) = stream.next().await {
        match event.ok()? {
            peri_model::ModelStreamEvent::TextDelta { text: delta } => {
                text.push_str(&delta);
                if text.len() > 64 {
                    return None;
                }
            }
            peri_model::ModelStreamEvent::Completed(response) => {
                if !matches!(response.stop_reason(), peri_model::StopReason::EndTurn) {
                    return None;
                }
                let answer = if text.is_empty() {
                    response.assistant_text()?
                } else {
                    text
                };
                return Some(
                    json!({"binding":binding,"value":answer.trim().parse::<i32>().ok()?,"request_model":model,"endpoint":url}),
                );
            }
            peri_model::ModelStreamEvent::Interrupted { .. }
            | peri_model::ModelStreamEvent::ToolCallDelta { .. } => return None,
            _ => {}
        }
    }
    None
}

#[tokio::main(flavor = "current_thread")]
async fn main() {
    if let Some(v) = run().await {
        println!("{v}");
    } else {
        std::process::exit(1);
    }
}

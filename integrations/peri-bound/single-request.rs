// Site-specific transport constraints for the one-call HelixWeave grant.
impl AnthropicModel {
    pub fn single_request(config: AnthropicConfig) -> ModelResult<Self> {
        let client = reqwest::Client::builder()
            .no_proxy()
            .redirect(reqwest::redirect::Policy::none())
            .retry(reqwest::retry::never())
            .timeout(std::time::Duration::from_secs(30))
            .build()
            .map_err(|_| ModelError::protocol(crate::ProtocolErrorKind::Provider))?;
        Ok(Self {
            config,
            transport: Arc::new(ReqwestTransport::new(client.clone())),
            client,
        })
    }
}

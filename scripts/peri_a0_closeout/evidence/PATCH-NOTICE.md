# 实验补丁声明

`candidate-testfix.patch` 是针对 Peri 提交 `c7d154c9985cbe08ba25e964cb7468d15014e6d0` 的实验差异，不是 Peri 官方版本，也不是 HelixWeave 决定长期维护的 fork。

Peri 在该提交使用 Apache License 2.0。分发或继续修改补丁所覆盖的 Peri 源码时，保留上游版权声明和 Apache-2.0 许可文本。本仓库不内置 Peri 源码树、构建缓存或运行 binary。

已接受补丁的 SHA256 是 `06272e1a66423a9e562509a38d044cd8506a1c10e61ef1fd3fe239107680587b`。`candidate-final.patch` 是同一实验里测试桩修正之前的补丁，生产差异与最终补丁一致，binary SHA256 同为 `60144bbe0919992b8913b9031696b2e30b1e9f82d95048ae3b07c887d5d85517`。细节见 `testfix-identity.json`。

from encode_pipeline.services.agent_redaction import RedactionPolicy
_HOME_PREFIX = "/" + "home"
_DATA_PREFIX = "/" + "data"
_WINDOWS_PREFIX = "C:" + "\\\\"
_WINDOWS_BACKSLASH = "\\"
def test_redacts_sensitive_keys():
    policy = RedactionPolicy()
    result = policy.redact_value({"api_key": "secret123", "name": "chipseq"})
    assert result.value == {"api_key": "<REDACTED>", "name": "chipseq"}
    assert result.sensitive_keys_redacted == 1
    assert result.absolute_paths_redacted == 0

def test_does_not_redact_padded_or_substring_sensitive_keys():
    policy = RedactionPolicy()
    result = policy.redact_value({
        "  password  ": "secret",
        "my_api_key": "secret",
        "credentials": "secret",
    })
    assert result.value == {
        "  password  ": "secret",
        "my_api_key": "secret",
        "credentials": "secret",
    }
    assert result.sensitive_keys_redacted == 0

def test_redacts_nested_sensitive_keys():
    policy = RedactionPolicy()
    result = policy.redact_value({"nested": {"api_key": "x"}})
    assert result.value == {"nested": {"api_key": "<REDACTED>"}}
    assert result.sensitive_keys_redacted == 1

def test_redacts_mixed_list_and_nested_paths():
    policy = RedactionPolicy()
    result = policy.redact_value(["/tmp/a", {"password": "x"}])
    assert result.value == ["<FILE_PATH>", {"password": "<REDACTED>"}]
    assert result.absolute_paths_redacted == 1
    assert result.sensitive_keys_redacted == 1

def test_non_string_scalars_pass_through():
    policy = RedactionPolicy()
    result = policy.redact_value({"count": 42, "empty": None, "flag": True})
    assert result.value == {"count": 42, "empty": None, "flag": True}
    assert result.sensitive_keys_redacted == 0
    assert result.absolute_paths_redacted == 0

def test_redacts_absolute_paths():
    policy = RedactionPolicy()
    path = _HOME_PREFIX + "/user/project/file.tsv"
    result = policy.redact_text(f"Load data from {path}.")
    assert result.value == "Load data from <FILE_PATH>."
    assert result.absolute_paths_redacted == 1
    assert result.sensitive_keys_redacted == 0

def test_redacts_unix_paths():
    policy = RedactionPolicy()
    home_path = _HOME_PREFIX + "/user/project/file.tsv"
    tmp_path = "/tmp/run/config.yaml"
    result = policy.redact_text(f"See {home_path} and {tmp_path}.")
    assert result.value == "See <FILE_PATH> and <FILE_PATH>."
    assert result.absolute_paths_redacted == 2

def test_redacts_unix_paths_with_spaces_unicode_and_hidden_files():
    policy = RedactionPolicy()
    tmp_path = "/tmp/my file.tsv"
    home_path = _HOME_PREFIX + "/用户/文件.tsv"
    bashrc_path = "/tmp/.bashrc"
    result = policy.redact_text(
        f"Path {tmp_path} is bad; visit {home_path} now; "
        f"config at {bashrc_path} or ..."
    )
    assert result.value == (
        "Path <FILE_PATH> is bad; visit <FILE_PATH> now; "
        "config at <FILE_PATH> or ..."
    )
    assert result.absolute_paths_redacted == 3

def test_does_not_redact_api_paths_or_issue_paths():
    policy = RedactionPolicy()
    result = policy.redact_text("POST /api/v1/workflows/x/agent/chat and field config.samples")
    assert result.value == "POST /api/v1/workflows/x/agent/chat and field config.samples"
    assert result.absolute_paths_redacted == 0

def test_does_not_redact_bare_api():
    policy = RedactionPolicy()
    result = policy.redact_text("GET /api")
    assert result.value == "GET /api"
    assert result.absolute_paths_redacted == 0

def test_redacts_windows_paths():
    policy = RedactionPolicy()
    local_path = _WINDOWS_PREFIX + "Users\\name\\file.tsv"
    unc_path = "\\\\server\\share\\file.tsv"
    result = policy.redact_text(
        f"Files at {local_path} and {unc_path}."
    )
    assert result.value == "Files at <FILE_PATH> and <FILE_PATH>."
    assert result.absolute_paths_redacted == 2

def test_redacts_windows_path_with_spaces():
    policy = RedactionPolicy()
    path = _WINDOWS_PREFIX + "Users\\My Name\\file.tsv"
    result = policy.redact_text(f"Run {path}")
    assert result.value == "Run <FILE_PATH>"
    assert result.absolute_paths_redacted == 1

def test_redacts_issue_paths_conditionally():
    policy = RedactionPolicy()
    issue = {
        "code": "FRIP_LOW",
        "message": "Path /tmp/run is wrong.",
        "severity": "error",
        "path": "config.samples",
        "source": "samples",
        "technical_message": None,
        "hint": None,
        "context": {"file": _HOME_PREFIX + "/data.tsv"},
    }
    result = policy.redact_issue(issue)
    assert result.value["path"] == "config.samples"
    assert result.value["message"] == "Path <FILE_PATH> is wrong."
    assert result.value["context"] == {"file": "<FILE_PATH>"}
    assert result.absolute_paths_redacted == 2

def test_redacts_filesystem_issue_path():
    policy = RedactionPolicy()
    issue = {
        "code": "FILE_NOT_FOUND",
        "message": "Missing.",
        "severity": "error",
        "path": _HOME_PREFIX + "/user/missing.tsv",
        "source": "runtime",
        "technical_message": None,
        "hint": None,
        "context": {},
    }
    result = policy.redact_issue(issue)
    assert result.value["path"] == "<FILE_PATH>"

def test_redaction_does_not_mutate_input():
    policy = RedactionPolicy()
    original = {"api_key": "secret", "nested": {"path": "/tmp/x"}}
    snapshot = {"api_key": "secret", "nested": {"path": "/tmp/x"}}
    policy.redact_value(original)
    assert original == snapshot

def test_redacts_deeply_nested_structure_without_recursion_error():
    policy = RedactionPolicy()
    nested: object = {"api_key": "secret"}
    for _ in range(200):
        nested = {"child": nested}
    result = policy.redact_value(nested)
    assert isinstance(result.value, dict)
    assert result.sensitive_keys_redacted == 1

def test_redacts_cyclic_structure_without_recursion_error():
    policy = RedactionPolicy()
    cyclic: dict[str, object] = {"api_key": "secret"}
    cyclic["self"] = cyclic
    result = policy.redact_value(cyclic)
    assert result.value["api_key"] == "<REDACTED>"
    assert result.value["self"] == "<CYCLIC_REFERENCE>"
    assert result.sensitive_keys_redacted == 1

def test_redacts_uppercase_word_after_path():
    policy = RedactionPolicy()
    result = policy.redact_text("See /tmp/a And then leave.")
    assert result.value == "See <FILE_PATH> And then leave."
    assert result.absolute_paths_redacted == 1

def test_redacts_beyond_max_depth_without_crash():
    policy = RedactionPolicy()
    nested: object = {"api_key": "secret"}
    for _ in range(250):
        nested = {"child": nested}
    result = policy.redact_value(nested)
    # The structure is redacted level by level until the max depth is reached,
    # at which point the remaining subtree is replaced with a placeholder.
    assert isinstance(result.value, dict)
    assert result.value is not nested
    deepest = result.value
    for _ in range(200):
        deepest = deepest["child"]
    assert deepest == {"child": "<MAX_DEPTH_REACHED>"}
    assert result.sensitive_keys_redacted == 0

def test_redacts_issues_list_preserves_codes_and_sources():
    policy = RedactionPolicy()
    issues = [
        {
            "code": "FRIP_LOW",
            "message": "Path /tmp/run is wrong.",
            "severity": "error",
            "path": "config.samples",
            "source": "samples",
            "technical_message": None,
            "hint": None,
            "context": {"file": _HOME_PREFIX + "/data.tsv"},
        },
        {
            "code": "FILE_NOT_FOUND",
            "message": f"Missing {_HOME_PREFIX + '/user/missing.tsv'}.",
            "severity": "error",
            "path": _HOME_PREFIX + "/user/missing.tsv",
            "source": "runtime",
            "technical_message": None,
            "hint": None,
            "context": {},
        },
    ]
    result = policy.redact_issues(issues)
    redacted = result.value
    assert redacted[0]["code"] == "FRIP_LOW"
    assert redacted[0]["severity"] == "error"
    assert redacted[0]["source"] == "samples"
    assert redacted[0]["path"] == "config.samples"
    assert redacted[0]["message"] == "Path <FILE_PATH> is wrong."
    assert redacted[0]["context"] == {"file": "<FILE_PATH>"}
    assert redacted[1]["code"] == "FILE_NOT_FOUND"
    assert redacted[1]["severity"] == "error"
    assert redacted[1]["source"] == "runtime"
    assert redacted[1]["path"] == "<FILE_PATH>"
    assert redacted[1]["message"] == "Missing <FILE_PATH>."
    assert result.absolute_paths_redacted == 4
    assert result.sensitive_keys_redacted == 0

def test_does_not_redact_api_followed_by_text():
    policy = RedactionPolicy()
    result = policy.redact_text("/api something and /api /tmp/x")
    assert result.value == "/api something and /api <FILE_PATH>"
    assert result.absolute_paths_redacted == 1

def test_handles_non_string_mapping_keys():
    policy = RedactionPolicy()
    result = policy.redact_value({1: "x", "api_key": "secret"})
    assert result.value == {1: "x", "api_key": "<REDACTED>"}
    assert result.sensitive_keys_redacted == 1

def test_does_not_redact_urls():
    policy = RedactionPolicy()
    result = policy.redact_text(
        "See https://example.com/path and http://example.com/path?x=1."
    )
    assert result.value == (
        "See https://example.com/path and http://example.com/path?x=1."
    )
    assert result.absolute_paths_redacted == 0

def test_does_not_redact_api_case_or_punctuation_variants():
    policy = RedactionPolicy()
    result = policy.redact_text("/API/v1 /Api?x=1 /api, /apiary /apiv1")
    assert result.value == "/API/v1 /Api?x=1 /api, <FILE_PATH>"
    assert result.absolute_paths_redacted == 1

def test_does_not_redact_windows_drive_lookalikes():
    policy = RedactionPolicy()
    result = policy.redact_text("id: x, note: y, path: /tmp/x, " + _WINDOWS_PREFIX + "tmp\\x")
    assert result.value == "id: x, note: y, path: <FILE_PATH>, <FILE_PATH>"
    assert result.absolute_paths_redacted == 2

def test_empty_string_and_container_passthrough():
    policy = RedactionPolicy()
    assert policy.redact_text("").value == ""
    assert policy.redact_value({}).value == {}
    assert policy.redact_value([]).value == []

def test_does_not_redact_relative_unix_paths():
    policy = RedactionPolicy()
    result = policy.redact_text(
        f"Load .{_DATA_PREFIX}/file.tsv, ..{_DATA_PREFIX}/file.tsv, ~/file, "
        f". /tmp/x, .. /tmp/x"
    )
    assert result.value == (
        f"Load .{_DATA_PREFIX}/file.tsv, ..{_DATA_PREFIX}/file.tsv, ~/file, "
        f". <FILE_PATH>, .. <FILE_PATH>"
    )
    assert result.absolute_paths_redacted == 2

def test_does_not_redact_embedded_windows_drive_letters():
    policy = RedactionPolicy()
    result = policy.redact_text(
        "abc:/tmp is not a drive; xD:" + _WINDOWS_BACKSLASH + "tmp is not a drive; but "
        + _WINDOWS_PREFIX + "tmp is."
    )
    assert result.value == (
        "abc:/tmp is not a drive; xD:" + _WINDOWS_BACKSLASH + "tmp is not a drive; but <FILE_PATH> is."
    )
    assert result.absolute_paths_redacted == 1

def test_path_with_multiple_space_separated_words_stops_at_plain_word():
    # The regex continues across spaces only while tokens look path-like
    # (contain "/", ".", or "-"). Plain words terminate the match.
    policy = RedactionPolicy()
    result = policy.redact_text("Path /tmp/my file name")
    assert result.value == "Path <FILE_PATH> file name"
    assert result.absolute_paths_redacted == 1

def test_preserves_delimiters_around_paths():
    policy = RedactionPolicy()
    result = policy.redact_text("(/tmp/x) \" /tmp/x \" [ /tmp/x ] { /tmp/x }")
    assert result.value == "(<FILE_PATH>) \" <FILE_PATH> \" [ <FILE_PATH> ] { <FILE_PATH> }"
    assert result.absolute_paths_redacted == 4

def test_preserves_parenthetical_text_after_path():
    policy = RedactionPolicy()
    result = policy.redact_text("Path /tmp/x (yes) is bad")
    assert result.value == "Path <FILE_PATH> (yes) is bad"
    assert result.absolute_paths_redacted == 1

def test_does_not_redact_protocol_relative_urls():
    policy = RedactionPolicy()
    result = policy.redact_text("see //example.com/path and //cdn.example.com/x.js")
    assert result.value == "see //example.com/path and //cdn.example.com/x.js"
    assert result.absolute_paths_redacted == 0

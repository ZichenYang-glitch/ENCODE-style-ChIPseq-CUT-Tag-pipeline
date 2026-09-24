"""H5: real platform chain for the Hi-TrAC adapter — API, Redis/RQ, worker.

Explicit ``real_execution`` tier; missing coordinates fail, never skip. The
chain is authenticated HTTP API -> validated snapshot -> create/submit run ->
dedicated real Redis/RQ -> original worker with the default runtime assembly
-> pinned tracPre2 with real scientific tools -> SQLite atomic artifact/QC
publication -> original API -> actual download. No in-process ASGI, no
fakeredis, no direct runner calls, no pre-filled result database.

The only injected control is a test-side SIGSTOP on the already-started
scientific process group (the same sanctioned seam as the H3 control
qualification): it holds real science in place so cancel/timeout land
deterministically. It never changes scientific bytes, assertions, or the
product. Unrelated sentinel processes must survive every scenario.

Additional coordinates beyond the H2/H3/H4 set:
  HELIXWEAVE_H5_REDIS_SERVER       redis-server executable (7.x)
  HELIXWEAVE_H5_REDIS_LD_LIBRARY   optional LD_LIBRARY_PATH for that binary
  HELIXWEAVE_H5_RESULTS            directory for the browser-phase handoff
"""

from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import signal
import sqlite3
import subprocess
import sys
import time

import pytest

from encode_pipeline.adapters.hitrac_preprocess.results import METRIC_KEYS
from test_qualification_real import (
    coordinates as shared_coordinates,
    expected_rows,
    rows,
)

from h5_stack import (
    ADMIN_PASSWORD,
    ADMIN_USERNAME,
    MEMBER_PASSWORD,
    MEMBER_USERNAME,
    WORKFLOW_ID,
    CaptureSmtp,
    HarnessError,
    RedisServer,
    Stack,
    identity_alive,
    process_children,
    process_tree,
    read_proc_identity,
    set_task_python,
)

pytestmark = pytest.mark.real_execution

# Re-export the shared real-coordinate fixture under its conventional name.
coordinates = shared_coordinates

RESULTS = (
    Path(os.environ.get("HELIXWEAVE_H5_RESULTS", ""))
    if os.environ.get("HELIXWEAVE_H5_RESULTS")
    else None
)

CANCEL_REASON = "User requested cancellation."


@pytest.fixture(scope="module")
def h5_coordinates(coordinates):
    missing = [
        key for key in ("HELIXWEAVE_H5_REDIS_SERVER",) if not os.environ.get(key)
    ]
    if missing:
        pytest.fail(f"H5 real platform coordinates missing: {missing}")
    values = dict(coordinates)
    values["REDIS_SERVER"] = os.environ["HELIXWEAVE_H5_REDIS_SERVER"]
    values["REDIS_LD_LIBRARY_PATH"] = os.environ.get("HELIXWEAVE_H5_REDIS_LD_LIBRARY")
    values["RUNTIME_BINDING_SHA256"] = hashlib.sha256(
        Path(values["RUNTIME_BINDING"]).read_bytes()
    ).hexdigest()
    set_task_python(os.environ.get("HELIXWEAVE_H5_PYTHON") or sys.executable)
    return values


@pytest.fixture(scope="module")
def redis_server(tmp_path_factory, h5_coordinates):
    root = tmp_path_factory.mktemp("h5-redis")
    server = RedisServer(
        root,
        h5_coordinates["REDIS_SERVER"],
        h5_coordinates["REDIS_LD_LIBRARY_PATH"],
    )
    yield server
    evidence = server.stop()
    _record(root, "redis-stop.json", evidence)


@pytest.fixture(scope="module")
def capture_smtp(tmp_path_factory):
    root = tmp_path_factory.mktemp("h5-smtp")
    sink = CaptureSmtp(root)
    yield sink
    evidence = sink.stop()
    _record(root, "smtp-stop.json", evidence)


def _make_stack(tmp_path, redis_server, capture_smtp, **kwargs):
    root = tmp_path / "stack"
    root.mkdir()
    instance = Stack(
        root,
        redis_server.url,
        queue_name=f"h5-{os.urandom(4).hex()}",
        smtp=capture_smtp,
        **kwargs,
    )
    return instance


@pytest.fixture()
def stack(tmp_path, redis_server, capture_smtp, h5_coordinates):
    instance = _make_stack(tmp_path, redis_server, capture_smtp)
    registration = instance.prepare_database(h5_coordinates)
    instance.start(h5_coordinates)
    instance.reference_registration = registration
    yield instance
    stopped = instance.stop()
    _record(tmp_path, "stack-stop.json", {"services": stopped})


def _record(root: Path, name: str, payload) -> None:
    root.mkdir(parents=True, exist_ok=True)
    (root / name).write_text(json.dumps(payload, indent=2, default=str) + "\n")


def _login(client, username=ADMIN_USERNAME, password=ADMIN_PASSWORD) -> None:
    response = client.post(
        "/api/v1/auth/login",
        json={"username": username, "password": password},
    )
    assert response.status_code == 200, response.text


def _csrf_headers(client) -> dict:
    for cookie in client.cookies.jar:
        if "csrf" in cookie.name.lower():
            return {"X-CSRF-Token": cookie.value}
    raise AssertionError("the session has no CSRF cookie")


def _post(client, url, **kwargs):
    return client.post(url, headers=_csrf_headers(client), **kwargs)


def _poll_status(client, run_id, wanted, timeout=300):
    deadline = time.monotonic() + timeout
    seen = []
    while time.monotonic() < deadline:
        response = client.get(f"/api/v1/runs/{run_id}")
        assert response.status_code == 200, response.text
        status = response.json()["run"]["status"]
        if not seen or seen[-1] != status:
            seen.append(status)
        if status == wanted:
            return seen
        if status in ("failed", "cancelled") and wanted not in (
            "failed",
            "cancelled",
        ):
            raise AssertionError(f"run {run_id} reached {status}: wanted {wanted}")
        time.sleep(0.2)
    raise AssertionError(f"run {run_id} never reached {wanted}; saw {seen}")


def _validate_and_create(client, revision_id, samples, options=None):
    """Validate -> create over the authenticated API; return run/snapshot IDs."""
    body = {
        "config": {},
        "samples": samples,
        "options": options or {"threads": 2, "mapq": 10},
        "reference_profile_revision_id": revision_id,
    }
    validated = _post(client, f"/api/v1/workflows/{WORKFLOW_ID}/validate", json=body)
    assert validated.status_code == 200, validated.text
    snapshot = validated.json()["snapshot"]
    assert snapshot["snapshot_id"], snapshot
    created = _post(
        client,
        f"/api/v1/workflows/{WORKFLOW_ID}/runs",
        json={"snapshot_id": snapshot["snapshot_id"]},
    )
    assert created.status_code == 201, created.text
    return created.json()["run"]["run_id"], snapshot["snapshot_id"]


def _submit_run(client, revision_id, samples, options=None):
    run_id, snapshot_id = _validate_and_create(client, revision_id, samples, options)
    preflight = _post(client, f"/api/v1/runs/{run_id}/preflight")
    assert preflight.status_code == 202, preflight.text
    _poll_status(client, run_id, "planned")
    started = _post(client, f"/api/v1/runs/{run_id}/start")
    assert started.status_code == 202, started.text
    return run_id, snapshot_id


def _tiny_samples(coordinates, sample_ids, scenario="positive"):
    directory = Path(coordinates["TINY_INPUTS"]) / scenario / "fastq"
    pairs = sorted(directory.glob("*_R1.fastq.gz"))
    assert len(pairs) == len(sample_ids)
    return [
        {
            "sample_id": sample_id,
            "fastq_1": str(path),
            "fastq_2": str(path.with_name(path.name.replace("_R1", "_R2"))),
        }
        for sample_id, path in zip(sample_ids, pairs)
    ]


def _download(client, run_id, artifact, generation):
    response = client.get(
        f"/api/v1/runs/{run_id}/artifacts/{artifact['artifact_id']}/download",
        params={"generation": generation, "revision": artifact["revision"]},
    )
    assert response.status_code == 200, response.text
    assert "attachment;" in response.headers["content-disposition"]
    return response.content


def _poll_artifacts(client, run_id, count, timeout=120):
    """Science SUCCEEDED precedes atomic publication; wait for the bundle."""
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        response = client.get(f"/api/v1/runs/{run_id}/artifacts", params={"limit": 100})
        assert response.status_code == 200, response.text
        artifacts = response.json()
        if len(artifacts["artifacts"]) == count:
            return artifacts
        time.sleep(0.5)
    raise AssertionError(f"run {run_id} never published {count} artifacts")


def _expected_metric_table(coordinates):
    return coordinates["design"]["scenarios"]["positive"][
        "conditional_expected_summary"
    ]


def _list_all_events(client, run_id):
    events = []
    cursor = None
    while True:
        params = {"limit": 100}
        if cursor:
            params["after"] = cursor
        response = client.get(f"/api/v1/runs/{run_id}/events", params=params)
        assert response.status_code == 200, response.text
        payload = response.json()
        events.extend(payload["events"])
        cursor = payload["next_cursor"]
        if not cursor:
            return events


def _sqlite_counts(database_path: Path, run_id: str) -> dict:
    """Read-only SQLite evidence; never a substitute for product writes."""
    connection = sqlite3.connect(f"file:{database_path}?mode=ro", uri=True, timeout=10)
    try:
        return {
            table: connection.execute(
                f"SELECT count(*) FROM {table} WHERE run_id=?", (run_id,)
            ).fetchone()[0]
            for table in (
                "run_artifacts",
                "run_qc_metrics",
                "artifact_publications",
            )
        }
    finally:
        connection.close()


def _sqlite_artifact_digests(database_path: Path, run_id: str) -> dict:
    """Read-only registered digests from the persisted artifact metadata.

    The public API projection is a deliberate whitelist without the digest;
    the registered sha256 lives in the persisted metadata only.
    """
    connection = sqlite3.connect(f"file:{database_path}?mode=ro", uri=True, timeout=10)
    try:
        rows = connection.execute(
            "SELECT artifact_id, artifact_metadata FROM run_artifacts WHERE run_id=?",
            (run_id,),
        ).fetchall()
        return {
            artifact_id: json.loads(metadata)["sha256"]
            for artifact_id, metadata in rows
        }
    finally:
        connection.close()


def _sqlite_result_state(database_path: Path, run_id: str) -> dict:
    """Read-only attempt/generation evidence for the delivered report."""
    connection = sqlite3.connect(f"file:{database_path}?mode=ro", uri=True, timeout=10)
    try:
        cursor = connection.execute(
            "SELECT * FROM run_result_states WHERE run_id=?", (run_id,)
        )
        names = [column[0] for column in cursor.description]
        rows = cursor.fetchall()
        if not rows:
            return {}
        state = dict(zip(names, rows[0]))
        return {
            key: state.get(key)
            for key in (
                "artifact_attempt_id",
                "artifact_attempt_status",
                "artifact_outcome",
                "artifact_generation",
                "qc_attempt_id",
                "qc_attempt_status",
                "qc_outcome",
                "qc_generation",
                "qc_artifact_generation",
            )
        }
    finally:
        connection.close()


def _sqlite_publication_generations(database_path: Path, run_id: str) -> list:
    connection = sqlite3.connect(f"file:{database_path}?mode=ro", uri=True, timeout=10)
    try:
        return sorted(
            row[0]
            for row in connection.execute(
                "SELECT DISTINCT artifact_generation FROM artifact_publications"
                " WHERE run_id=?",
                (run_id,),
            )
        )
    finally:
        connection.close()


def _sha256(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _succeeded_mail(sink: CaptureSmtp, run_id: str) -> list[dict]:
    return [
        message
        for message in sink.captured()
        if f"HelixWeave run SUCCEEDED: {run_id}" in message["data"]
    ]


def _mail_for(sink: CaptureSmtp, run_id: str) -> list[str]:
    subjects = []
    for message in sink.captured():
        if run_id in message["data"]:
            for line in message["data"].splitlines():
                if line.lower().startswith("subject:"):
                    subjects.append(line.split(":", 1)[1].strip())
    return subjects


def _wait_science_start(stack: Stack, run_id: str, timeout=120):
    """Wait for one real scientific delegate; then hold it with SIGSTOP.

    The control is the sanctioned H3 seam: the original tool has actually
    started, so the outer stop cannot race a finished run. It never alters
    bytes, arguments, or assertions.
    """
    workspace = stack.workspace_root / run_id
    directory = workspace / "hitrac-attempt/private/calls"
    worker = stack.worker_process()
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline:
        for path in sorted(directory.glob("*.start.json")):
            try:
                receipt = json.loads(path.read_text())
            except (json.JSONDecodeError, FileNotFoundError):
                continue
            if receipt["stage"] != "align" or not receipt["sample"]:
                continue
            wrapper = read_proc_identity(receipt["pid"])
            if wrapper is None:
                continue
            delegate_tree = process_tree(receipt["pid"])
            if len(delegate_tree) < 2:
                continue
            science_group = wrapper["group"]
            leader = read_proc_identity(science_group)
            if leader is None:
                continue
            worker_children = process_children(worker.process.pid)
            horse_identities = [
                identity
                for child in worker_children
                if (identity := read_proc_identity(child)) is not None
            ]
            if not horse_identities:
                continue
            os.killpg(science_group, signal.SIGSTOP)
            held_tree = process_tree(receipt["pid"])
            # Capture the full observation set while the tree is held: the
            # tool subtree, every horse subtree, and the separately recorded
            # science group leader, deduplicated by PID.
            horse_trees = []
            for horse in horse_identities:
                horse_trees.extend(process_tree(horse["pid"]))
            captured = {}
            for identity in [leader] + held_tree + horse_identities + horse_trees:
                captured[identity["pid"]] = identity
            return {
                "receipt_path": str(path),
                "receipt": receipt,
                "science_group_leader": leader,
                "held_tree": held_tree,
                "worker_horse_identities": horse_identities,
                "captured_tree_identities": list(captured.values()),
                "injected_signal": "SIGSTOP",
                "injected_monotonic": time.monotonic(),
            }
        if worker.process.poll() is not None:
            raise HarnessError("worker exited before science started")
        time.sleep(0.01)
    raise AssertionError("no scientific align receipt before the deadline")


def _assert_identities_gone(observation, timeout=30):
    """The held tree and every captured identity are gone before cleanup.

    Runs before the test's own finally restoration/cleanup, so the product
    alone is responsible for the reaping. Identity matching uses PID +
    starttime, which guards against PID reuse; a zombie still counts as
    present, never as reaped. The captured set is finite: it covers the held
    tool subtree, the horse subtrees, and the science group leader observed at
    capture time, and is not a claim that the product can recover any
    arbitrary rootless descendant.
    """
    identities = list(observation["held_tree"]) + list(
        observation["worker_horse_identities"]
    )
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline and any(map(identity_alive, identities)):
        time.sleep(0.02)
    survivors = [item["pid"] for item in identities if identity_alive(item)]
    assert not survivors, f"process identities still alive: {survivors}"
    captured = list(observation["captured_tree_identities"])
    deadline = time.monotonic() + timeout
    while time.monotonic() < deadline and any(map(identity_alive, captured)):
        time.sleep(0.02)
    captured_survivors = [item["pid"] for item in captured if identity_alive(item)]
    assert not captured_survivors, (
        f"captured process identities still alive: {captured_survivors}"
    )


def _resume_and_reap(observation):
    """Test-owned cleanup for the injected hold; never touches others."""
    if not observation:
        return
    leader = observation.get("science_group_leader") or {}
    leader_pid = leader.get("pid")
    if isinstance(leader_pid, int) and leader_pid > 1 and identity_alive(leader):
        try:
            os.killpg(leader_pid, signal.SIGCONT)
        except (ProcessLookupError, PermissionError):
            pass
    for identity in observation.get("held_tree", []):
        if identity_alive(identity):
            try:
                os.kill(identity["pid"], signal.SIGKILL)
            except (ProcessLookupError, PermissionError):
                pass


def _sentinel(tmp_path):
    process = subprocess.Popen(
        [sys.executable, "-I", "-B", "-c", "import time; time.sleep(600)"],
        cwd=tmp_path,
        start_new_session=True,
    )
    return process, read_proc_identity(process.pid)


def _assert_run_private_retained(stack: Stack, run_id: str) -> dict:
    workspace = stack.workspace_root / run_id
    attempt = workspace / "hitrac-attempt"
    evidence = {"workspace": str(workspace)}
    assert not (attempt / "complete.json").exists()
    outcome = attempt / "private/outcome.json"
    if outcome.exists():
        evidence["private_outcome"] = json.loads(outcome.read_text())
        assert evidence["private_outcome"]["status"] != "complete"
    return evidence


def test_h5_real_platform_chain(stack, h5_coordinates, capture_smtp, tmp_path):
    """Authenticated API -> RQ/worker -> real science -> SQLite -> download."""
    sample_ids = ("a  b", "a b ")
    observations = {"sample_ids": sample_ids}
    revision_id = stack.reference_registration["revision_id"]
    with stack.client() as client:
        _login(client)
        listed = client.get("/api/v1/workflows/")
        assert listed.status_code == 200, listed.text
        hitrac = next(
            w
            for w in listed.json()["workflows"]
            if w["metadata"]["workflow_id"] == WORKFLOW_ID
        )
        assert hitrac["availability"]["execution"] == "available"

        samples = _tiny_samples(h5_coordinates, sample_ids)
        input_sha = {
            key: _sha256(Path(key).read_bytes())
            for row in samples
            for key in (row["fastq_1"], row["fastq_2"])
        }
        run_id, snapshot_id = _submit_run(client, revision_id, samples)
        observations.update(run_id=run_id, snapshot_id=snapshot_id)
        transitions = _poll_status(client, run_id, "succeeded", timeout=600)
        observations["status_transitions"] = transitions

        artifacts = _poll_artifacts(client, run_id, 5)
        qc_response = client.get(
            f"/api/v1/runs/{run_id}/qc-metrics", params={"limit": 100}
        )
        assert qc_response.status_code == 200
        qc = qc_response.json()
        assert len(artifacts["artifacts"]) == 5
        assert len(qc["qc_metrics"]) == 30
        assert artifacts["artifact_generation"] and qc["qc_generation"]
        observations["artifact_generation"] = artifacts["artifact_generation"]
        observations["qc_generation"] = qc["qc_generation"]
        assert _sqlite_counts(stack.root / "platform.db", run_id) == {
            "run_artifacts": 5,
            "run_qc_metrics": 30,
            "artifact_publications": 5,
        }
        result_state = _sqlite_result_state(stack.root / "platform.db", run_id)
        assert result_state["artifact_attempt_status"] == "succeeded"
        assert result_state["qc_attempt_status"] == "succeeded"
        assert (
            result_state["artifact_outcome"]
            == result_state["qc_outcome"]
            == ("succeeded")
        )
        assert result_state["artifact_generation"] == artifacts["artifact_generation"]
        assert result_state["qc_generation"] == qc["qc_generation"]
        assert (
            result_state["qc_artifact_generation"] == (artifacts["artifact_generation"])
        )
        observations["result_state"] = result_state
        # Append-only publications are bound to the same generation.
        assert _sqlite_publication_generations(stack.root / "platform.db", run_id) == [
            artifacts["artifact_generation"]
        ]
        observations["publication_generation"] = artifacts["artifact_generation"]

        wanted = _expected_metric_table(h5_coordinates)
        assert {m["sample_id"] for m in qc["qc_metrics"]} == set(sample_ids)
        for sample_id in sample_ids:
            table = {
                m["metric_key"]: m
                for m in qc["qc_metrics"]
                if m["sample_id"] == sample_id
            }
            assert set(table) == set(METRIC_KEYS)
            assert [float(table[key]["value"]) for key in METRIC_KEYS] == (
                pytest.approx(wanted, rel=0, abs=1e-12)
            )

        workspace = stack.workspace_root / run_id
        output = workspace / "hitrac-attempt/output"
        expected = expected_rows(h5_coordinates["design"])
        complete = json.loads((workspace / "hitrac-attempt/complete.json").read_text())
        observations["complete_identity_sha256"] = complete["identity"]["sha256"]
        tokens = {}
        for token, sample in complete["results"]["samples"].items():
            assert sample["all"] == 8 and sample["noBg"] == 5
            assert sample["metrics"] == pytest.approx(wanted, rel=0, abs=1e-12)
            assert Counter(
                map(tuple, rows(output / token / f"{token}_all.bedpe.gz"))
            ) == Counter(expected)
            assert Counter(
                map(tuple, rows(output / token / f"{token}_unique.bedpe.gz"))
            ) == Counter(expected[i] for i in (0, 3, 4, 5, 6))
        # Internal tokens map back to the exact user sample IDs: positional
        # order binds s000001.. to the submitted row order.
        for artifact in artifacts["artifacts"]:
            name = artifact["name"]
            if not name.endswith(".bedpe.gz"):
                continue
            token = name.split("_")[0]
            position = int(token[1:]) - 1
            assert artifact["metadata"]["sample_id"] == sample_ids[position]
            tokens.setdefault(token, artifact["metadata"]["sample_id"])
        assert sorted(tokens.values()) == sorted(sample_ids)

        registered = _sqlite_artifact_digests(stack.root / "platform.db", run_id)
        assert set(registered) == {a["artifact_id"] for a in artifacts["artifacts"]}
        downloads = {}
        for artifact in artifacts["artifacts"]:
            content = _download(
                client, run_id, artifact, artifacts["artifact_generation"]
            )
            assert len(content) == artifact["size_bytes"]
            downloads[artifact["artifact_id"]] = _sha256(content)
            # The digest registered in the persisted artifact metadata matches
            # the downloaded bytes (the public API projection whitelists it
            # out, so the registered digest is read from the same database).
            assert (
                downloads[artifact["artifact_id"]]
                == registered[artifact["artifact_id"]]
            )
            # Downloaded bytes equal the projected original workspace bytes.
            projected = workspace / artifact["relative_path"]
            assert content == projected.read_bytes()
            if artifact["name"].endswith(".bedpe.gz"):
                continue
            assert artifact["name"] == "tracPre_summary.txt"
            # The summary is the original upstream file: 15 metrics + sample
            # index = 16 physical TSV columns, indexed by internal token. The
            # platform maps tokens back to exact user sample IDs only in
            # metadata/QC coordinates; the original bytes are never rewritten.
            text = content.decode("ascii")
            table = [line.split("\t") for line in text.splitlines()]
            assert all(len(row) == 16 for row in table)
            assert table[0][0] == ""
            assert sorted(row[0] for row in table[1:]) == sorted(tokens)
        observations["download_sha256"] = downloads

        # Private artifacts are never listed or downloadable.
        for private in ("bam", "fastq", "private-request", "call-plan"):
            denied = client.get(
                f"/api/v1/runs/{run_id}/artifacts/{private}/download",
                params={
                    "generation": artifacts["artifact_generation"],
                    "revision": artifacts["artifacts"][0]["revision"],
                },
            )
            assert denied.status_code == 404
        names = json.dumps(artifacts) + json.dumps(qc)
        assert str(stack.root) not in names
        assert str(tmp_path) not in names

        # Unauthenticated requests are refused.
        with stack.client() as anon:
            for suffix in ("artifacts", "qc-metrics"):
                refused = anon.get(f"/api/v1/runs/{run_id}/{suffix}")
                assert refused.status_code == 401

        # A member can read shared results but cannot use admin routes.
        created_member = _post(
            client,
            "/api/v1/auth/accounts",
            json={"username": MEMBER_USERNAME, "password": MEMBER_PASSWORD},
        )
        assert created_member.status_code == 200, created_member.text
        member_id = created_member.json()["account"]["user_id"]
        with stack.client() as member:
            response = member.post(
                "/api/v1/auth/login",
                json={"username": MEMBER_USERNAME, "password": MEMBER_PASSWORD},
            )
            assert response.status_code == 200
            assert member.get(f"/api/v1/runs/{run_id}/artifacts").status_code == 200
            assert member.get("/api/v1/auth/accounts").status_code == 403
            # A disabled member loses access entirely.
            disabled = _post(
                client,
                f"/api/v1/auth/accounts/{member_id}/status",
                json={"enabled": False},
            )
            assert disabled.status_code == 200, disabled.text
            assert member.get(f"/api/v1/runs/{run_id}/artifacts").status_code == 401

        # Exactly one success notification, carrying this run's QC summary.
        deadline = time.monotonic() + 30
        while time.monotonic() < deadline and not _succeeded_mail(capture_smtp, run_id):
            time.sleep(0.2)
        successes = _succeeded_mail(capture_smtp, run_id)
        assert len(successes) == 1
        body = successes[0]["data"]
        assert "Persisted QC summary:" in body
        assert body.count("\n- ") == 12
        assert "18 more metrics" in body
        assert successes[0]["recipients"] == ["rcpt TO:<h5-review@example.test>"]
        observations["success_email"] = {
            "subject_line": "HelixWeave run SUCCEEDED",
            "metric_lines": 12,
            "recipients": successes[0]["recipients"],
        }

        # Idempotent replay: the same snapshot returns the canonical run, and a
        # repeated start of a terminal run changes nothing.
        replayed = _post(
            client,
            f"/api/v1/workflows/{WORKFLOW_ID}/runs",
            json={"snapshot_id": snapshot_id},
        )
        assert replayed.status_code == 200, replayed.text
        assert replayed.json()["run"]["run_id"] == run_id
        restarted = _post(client, f"/api/v1/runs/{run_id}/start")
        assert restarted.status_code in (200, 202), restarted.text
        assert restarted.json()["run"]["status"] == "succeeded"
        assert _sqlite_counts(stack.root / "platform.db", run_id)["run_artifacts"] == 5

        # Identical legal inputs create a new run with a fresh attempt; the
        # original run and its results are never overwritten.
        rerun_id, rerun_snapshot = _validate_and_create(client, revision_id, samples)
        assert rerun_id != run_id and rerun_snapshot != snapshot_id
        preflight = _post(client, f"/api/v1/runs/{rerun_id}/preflight")
        assert preflight.status_code == 202, preflight.text
        _poll_status(client, rerun_id, "planned")
        started = _post(client, f"/api/v1/runs/{rerun_id}/start")
        assert started.status_code == 202, started.text
        _poll_status(client, rerun_id, "succeeded", timeout=600)
        rerun_artifacts = _poll_artifacts(client, rerun_id, 5)
        assert len(rerun_artifacts["artifacts"]) == 5
        assert (
            rerun_artifacts["artifact_generation"] != artifacts["artifact_generation"]
        )
        rerun_workspace = stack.workspace_root / rerun_id
        rerun_complete = json.loads(
            (rerun_workspace / "hitrac-attempt/complete.json").read_text()
        )
        rerun_output = rerun_workspace / "hitrac-attempt/output"
        for token, sample in rerun_complete["results"]["samples"].items():
            assert sample["all"] == 8 and sample["noBg"] == 5
            assert Counter(
                map(tuple, rows(rerun_output / token / f"{token}_all.bedpe.gz"))
            ) == Counter(expected)
        # The original run's results still resolve to the original generation.
        again = client.get(f"/api/v1/runs/{run_id}/artifacts", params={"limit": 100})
        assert again.json() == artifacts
        first_artifact = artifacts["artifacts"][0]
        assert (
            _download(client, run_id, first_artifact, artifacts["artifact_generation"])
            == (workspace / first_artifact["relative_path"]).read_bytes()
        )
        deadline = time.monotonic() + 30
        while (
            time.monotonic() < deadline
            and len(_succeeded_mail(capture_smtp, rerun_id)) < 1
        ):
            time.sleep(0.2)
        assert len(_succeeded_mail(capture_smtp, rerun_id)) == 1
        observations["rerun"] = {
            "run_id": rerun_id,
            "snapshot_id": rerun_snapshot,
            "artifact_generation": rerun_artifacts["artifact_generation"],
        }

        # Original inputs remain byte-identical and read-only.
        assert {key: _sha256(Path(key).read_bytes()) for key in input_sha} == input_sha
        observations["input_sha256_unchanged"] = True

    # The run/job/attempt evidence is durable: close every service (API and
    # worker both die, closing their SQLite handles), then read the same
    # database through a fresh API process.
    stopped = stack.stop()
    observations["first_stack_stop"] = stopped
    stack.start_api_only(h5_coordinates)
    try:
        with stack.client() as reopened:
            _login(reopened)
            detail = reopened.get(f"/api/v1/runs/{run_id}")
            assert detail.status_code == 200
            assert detail.json()["run"]["status"] == "succeeded"
            assert detail.json()["run"]["ended_at"] is not None
            listed = reopened.get(
                f"/api/v1/runs/{run_id}/artifacts", params={"limit": 100}
            ).json()
            assert listed == artifacts
            metrics = reopened.get(
                f"/api/v1/runs/{run_id}/qc-metrics", params={"limit": 100}
            ).json()
            assert metrics == qc
            content = _download(
                reopened, run_id, first_artifact, artifacts["artifact_generation"]
            )
            assert _sha256(content) == downloads[first_artifact["artifact_id"]]
            rerun_detail = reopened.get(f"/api/v1/runs/{rerun_id}")
            assert rerun_detail.json()["run"]["status"] == "succeeded"
        observations["reopen_verified"] = True
    finally:
        reopened_stop = stack.stop()
    observations["second_stack_stop"] = reopened_stop

    _record(tmp_path, "chain-observations.json", observations)
    if RESULTS is not None:
        RESULTS.mkdir(parents=True, exist_ok=True)
        (RESULTS / "browser-source.json").write_text(
            json.dumps(
                {
                    "database_url": stack.database_url,
                    "workspace_root": str(stack.workspace_root),
                    "reference_config": str(stack.root / "reference-profiles.json"),
                    "runtime_binding": h5_coordinates["RUNTIME_BINDING"],
                    "runtime_sha256": h5_coordinates["RUNTIME_BINDING_SHA256"],
                    "run_id": run_id,
                    "samples": list(sample_ids),
                    "artifact_sha256": downloads,
                    "admin": {"username": ADMIN_USERNAME, "password": ADMIN_PASSWORD},
                },
                indent=2,
            )
            + "\n"
        )
        (RESULTS / "browser-source.json").chmod(0o600)


def test_h5_real_platform_cancel(stack, h5_coordinates, capture_smtp, tmp_path):
    """Real API cancel -> RQ stop -> worker monitor -> SQLite acknowledgement."""
    observations = {}
    sentinel, sentinel_identity = _sentinel(tmp_path)
    revision_id = stack.reference_registration["revision_id"]
    try:
        with stack.client() as client:
            _login(client)
            samples = _tiny_samples(h5_coordinates, ("cancel one", "cancel two"))
            run_id, snapshot_id = _submit_run(client, revision_id, samples)
            observations.update(run_id=run_id, snapshot_id=snapshot_id)
            _poll_status(client, run_id, "running", timeout=120)
            science = _wait_science_start(stack, run_id)
            observations["science"] = science
            cancelled = _post(client, f"/api/v1/runs/{run_id}/cancel")
            assert cancelled.status_code == 202, cancelled.text
            observations["cancel_http"] = cancelled.status_code
            transitions = _poll_status(client, run_id, "cancelled", timeout=180)
            observations["status_transitions"] = transitions

            detail = client.get(f"/api/v1/runs/{run_id}").json()["run"]
            assert detail["status"] == "cancelled"
            assert detail["cancellation_reason"] == CANCEL_REASON
            assert detail["ended_at"] is not None
            events = _list_all_events(client, run_id)
            types = [event["event_type"] for event in events]
            observations["event_sequence"] = [
                {
                    "sequence": event["sequence"],
                    "event_type": event["event_type"],
                    "status": event["status"],
                }
                for event in events
            ]
            assert "cancellation_requested" in types
            assert "cancellation_acknowledged" in types
            assert types.index("cancellation_requested") < types.index(
                "cancellation_acknowledged"
            )
            acknowledged = events[types.index("cancellation_acknowledged")]
            job_id = acknowledged["context"]["job_id"]
            observations["job_id"] = job_id
            assert acknowledged["status"] == "cancelled"
            assert acknowledged["context"]["backend"] == "rq"
            assert acknowledged["context"]["queue_name"] == stack.queue_name

            # The matched horse and the held nested science tree are reaped;
            # the unrelated sentinel is untouched.
            _assert_identities_gone(science)
            observations["processes_gone"] = True
            assert identity_alive(sentinel_identity)
            observations["unrelated_sentinel_alive"] = True

            # RQ's own backend evidence agrees the job stopped.
            from rq.job import Job
            from redis import Redis

            connection = Redis.from_url(stack.redis_url, socket_timeout=2)
            try:
                job = Job.fetch(job_id, connection=connection)
                observations["rq_job"] = {
                    "id": job.id,
                    "status": job.get_status(refresh=True).value,
                    "origin": job.origin,
                }
                assert job.get_status().value == "stopped"
                assert job.origin == stack.queue_name
            finally:
                connection.close()

            # No result publication of any kind for the cancelled attempt.
            assert _sqlite_counts(stack.root / "platform.db", run_id) == {
                "run_artifacts": 0,
                "run_qc_metrics": 0,
                "artifact_publications": 0,
            }
            observations["result_rows"] = 0
            observations.update(_assert_run_private_retained(stack, run_id))

            # No success notification; the terminal cancel mail exists instead.
            deadline = time.monotonic() + 30
            while time.monotonic() < deadline and not _mail_for(capture_smtp, run_id):
                time.sleep(0.2)
            subjects = _mail_for(capture_smtp, run_id)
            observations["mail_subjects"] = subjects
            assert subjects == [f"HelixWeave run CANCELLED: {run_id}"]
            assert not _succeeded_mail(capture_smtp, run_id)

            input_sha = {
                key: _sha256(Path(key).read_bytes())
                for row in samples
                for key in (row["fastq_1"], row["fastq_2"])
            }
            observations["input_sha256"] = input_sha
    finally:
        _resume_and_reap(observations.get("science"))
        sentinel.terminate()
        sentinel.wait(timeout=5)

    # Persisted cancellation survives a full close: reopen via a fresh API.
    stopped = stack.stop()
    observations["first_stack_stop"] = stopped
    stack.start_api_only(h5_coordinates)
    try:
        with stack.client() as reopened:
            _login(reopened)
            detail = reopened.get(f"/api/v1/runs/{run_id}").json()["run"]
            assert detail["status"] == "cancelled"
            assert detail["cancellation_reason"] == CANCEL_REASON
            events = _list_all_events(reopened, run_id)
            assert [event["event_type"] for event in events] == types
            refused = reopened.get(
                f"/api/v1/runs/{run_id}/artifacts", params={"limit": 100}
            )
            assert refused.status_code == 200
            assert refused.json()["artifacts"] == []
        observations["reopen_verified"] = True
    finally:
        observations["second_stack_stop"] = stack.stop()
    _record(tmp_path, "cancel-observations.json", observations)


def test_h5_real_platform_timeout(tmp_path, redis_server, capture_smtp, h5_coordinates):
    """The existing job-timeout contract fires inside real nested science."""
    stack = _make_stack(tmp_path, redis_server, capture_smtp, job_timeout_seconds=30)
    registration = stack.prepare_database(h5_coordinates)
    stack.start(h5_coordinates)
    observations = {"job_timeout_seconds": 30}
    sentinel, sentinel_identity = _sentinel(tmp_path)
    try:
        with stack.client() as client:
            _login(client)
            samples = _tiny_samples(h5_coordinates, ("timeout one", "timeout two"))
            run_id, _ = _submit_run(client, registration["revision_id"], samples)
            observations["run_id"] = run_id
            _poll_status(client, run_id, "running", timeout=120)
            science = _wait_science_start(stack, run_id)
            observations["science"] = science
            transitions = _poll_status(client, run_id, "failed", timeout=180)
            observations["status_transitions"] = transitions

            detail = client.get(f"/api/v1/runs/{run_id}").json()["run"]
            assert detail["status"] == "failed"
            assert detail["error"] is not None
            assert detail["error"]["code"] == "RUN_EXECUTION_FAILED"
            observations["error"] = detail["error"]
            events = _list_all_events(client, run_id)
            failed_event = next(
                event
                for event in events
                if event["event_type"] == "status_changed"
                and event["status"] == "failed"
            )
            assert failed_event["issue"]["context"]["reason_code"] == (
                "PROCESS_RUNNER_TIMEOUT"
            )
            assert failed_event["context"]["reason_code"] == ("PROCESS_RUNNER_TIMEOUT")
            observations["failed_event"] = failed_event
            observations["event_sequence"] = [
                {
                    "sequence": event["sequence"],
                    "event_type": event["event_type"],
                    "status": event["status"],
                }
                for event in events
            ]
            assert "cancellation_requested" not in [
                event["event_type"] for event in events
            ]

            _assert_identities_gone(science)
            observations["processes_gone"] = True
            assert identity_alive(sentinel_identity)
            observations["unrelated_sentinel_alive"] = True
            assert _sqlite_counts(stack.root / "platform.db", run_id) == {
                "run_artifacts": 0,
                "run_qc_metrics": 0,
                "artifact_publications": 0,
            }
            observations.update(_assert_run_private_retained(stack, run_id))
            deadline = time.monotonic() + 30
            while time.monotonic() < deadline and not _mail_for(capture_smtp, run_id):
                time.sleep(0.2)
            subjects = _mail_for(capture_smtp, run_id)
            observations["mail_subjects"] = subjects
            assert subjects == [f"HelixWeave run FAILED: {run_id}"]
            assert not _succeeded_mail(capture_smtp, run_id)
    finally:
        _resume_and_reap(observations.get("science"))
        sentinel.terminate()
        sentinel.wait(timeout=5)
        observations["stack_stop"] = stack.stop()
    _record(tmp_path, "timeout-observations.json", observations)


@pytest.mark.parametrize(
    "scenario,sample_ids,reason_code,script_exit",
    [
        (
            "multisample_one_empty",
            ("multi one", "multi two"),
            "empty_pet_set",
            0,
        ),
        ("all_trimmed", ("trimmed only",), "call_missing", 1),
    ],
    ids=["policy-a-one-sample-empty", "original-script-exit-1"],
)
def test_h5_real_platform_science_rejected_never_publishes_partial(
    stack,
    h5_coordinates,
    capture_smtp,
    tmp_path,
    scenario,
    sample_ids,
    reason_code,
    script_exit,
):
    """Policy A and original-script failure refuse the whole attempt."""
    observations = {"scenario": scenario}
    revision_id = stack.reference_registration["revision_id"]
    with stack.client() as client:
        _login(client)
        samples = _tiny_samples(h5_coordinates, sample_ids, scenario=scenario)
        input_sha = {
            key: _sha256(Path(key).read_bytes())
            for row in samples
            for key in (row["fastq_1"], row["fastq_2"])
        }
        run_id, _ = _submit_run(client, revision_id, samples)
        observations["run_id"] = run_id
        _poll_status(client, run_id, "failed", timeout=600)
        detail = client.get(f"/api/v1/runs/{run_id}").json()["run"]
        assert detail["status"] == "failed"
        assert detail["error"]["code"] == "RUN_EXECUTION_FAILED"
        observations["error"] = detail["error"]

        # No partial success: zero artifacts, QC, and publications; no
        # completion marker; the original private outputs are retained.
        assert _sqlite_counts(stack.root / "platform.db", run_id) == {
            "run_artifacts": 0,
            "run_qc_metrics": 0,
            "artifact_publications": 0,
        }
        workspace = stack.workspace_root / run_id
        attempt = workspace / "hitrac-attempt"
        assert not (attempt / "complete.json").exists()
        outcome = json.loads((attempt / "private/outcome.json").read_text())
        observations["private_outcome"] = outcome
        assert outcome["status"] == "rejected"
        assert outcome["reason_code"] == reason_code
        execution = json.loads((attempt / "private/execution.json").read_text())
        assert execution["returncode"] == script_exit
        observations["script_returncode"] = script_exit
        # Everything the original script produced is retained, privately.
        assert (attempt / "output").is_dir()
        assert any((attempt / "output").iterdir())
        assert (attempt / "private/upstream.stdout").is_file()
        observations["private_outputs_retained"] = True

        listed = client.get(f"/api/v1/runs/{run_id}/artifacts", params={"limit": 100})
        qc = client.get(f"/api/v1/runs/{run_id}/qc-metrics", params={"limit": 100})
        assert listed.json()["artifacts"] == []
        assert qc.json()["qc_metrics"] == []
        subjects = _mail_for(capture_smtp, run_id)
        deadline = time.monotonic() + 30
        while time.monotonic() < deadline and not subjects:
            time.sleep(0.2)
            subjects = _mail_for(capture_smtp, run_id)
        observations["mail_subjects"] = subjects
        assert subjects == [f"HelixWeave run FAILED: {run_id}"]
        assert {key: _sha256(Path(key).read_bytes()) for key in input_sha} == input_sha
        observations["input_sha256_unchanged"] = True
    _record(tmp_path, "rejected-observations.json", observations)

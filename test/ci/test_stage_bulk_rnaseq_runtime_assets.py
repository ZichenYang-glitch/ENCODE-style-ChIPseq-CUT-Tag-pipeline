"""Behavioral contract for the explicit bulk RNA-seq asset staging tool."""

from __future__ import annotations

from dataclasses import replace
import gzip
import hashlib
import io
import json
from pathlib import Path
from pathlib import PurePosixPath
import stat
import tarfile
import zipfile

from jsonschema import Draft202012Validator
import pytest

import scripts.stage_bulk_rnaseq_runtime_assets as staging


def _sha256(content: bytes) -> str:
    return hashlib.sha256(content).hexdigest()


def _json_bytes(value: object) -> bytes:
    return (json.dumps(value, sort_keys=True, separators=(",", ":")) + "\n").encode()


def _tree_digest(entries: list[dict[str, object]]) -> str:
    ordered = sorted(entries, key=lambda entry: PurePosixPath(str(entry["path"])))
    return _sha256(json.dumps(ordered, sort_keys=True, separators=(",", ":")).encode())


def test_tree_identity_uses_the_production_path_component_order() -> None:
    entries = [
        {"path": "a.b", "size_bytes": 1, "sha256": "1" * 64},
        {"path": "a/file", "size_bytes": 1, "sha256": "2" * 64},
    ]

    assert staging._canonical_tree_digest(entries) == _tree_digest(entries)


def _tar_bytes(files: dict[str, tuple[bytes, int]], *, compressed: bool) -> bytes:
    output = io.BytesIO()
    if compressed:
        with gzip.GzipFile(fileobj=output, mode="wb", mtime=0) as compressed_output:
            with tarfile.open(fileobj=compressed_output, mode="w") as package:
                _add_tar_files(package, files)
    else:
        with tarfile.open(fileobj=output, mode="w") as package:
            _add_tar_files(package, files)
    return output.getvalue()


def _add_tar_files(
    package: tarfile.TarFile,
    files: dict[str, tuple[bytes, int]],
) -> None:
    for name, (content, mode) in sorted(files.items()):
        member = tarfile.TarInfo(name)
        member.size = len(content)
        member.mode = mode
        member.mtime = 0
        member.uid = 0
        member.gid = 0
        member.uname = ""
        member.gname = ""
        package.addfile(member, io.BytesIO(content))


def _zip_bytes(files: dict[str, bytes]) -> bytes:
    output = io.BytesIO()
    with zipfile.ZipFile(output, mode="w") as package:
        for name, content in sorted(files.items()):
            entry = zipfile.ZipInfo(name)
            entry.date_time = (1980, 1, 1, 0, 0, 0)
            entry.external_attr = 0o100644 << 16
            package.writestr(entry, content)
    return output.getvalue()


class _FakeRunner:
    def __init__(
        self,
        *,
        downloads: dict[str, bytes],
        image_coordinate: str,
        index_manifest: bytes,
        selected_manifest: bytes,
        image_digest: str,
        image_config_digest: str,
        docker_archive: bytes,
    ) -> None:
        self.downloads = downloads
        self.image_coordinate = image_coordinate
        self.index_manifest = index_manifest
        self.selected_manifest = selected_manifest
        self.image_digest = image_digest
        self.image_config_digest = image_config_digest
        self.docker_archive = docker_archive
        self.calls: list[tuple[tuple[str, ...], bool]] = []

    def __call__(
        self,
        argv: tuple[str, ...],
        *,
        capture_stdout: bool,
    ) -> bytes:
        self.calls.append((argv, capture_stdout))
        if argv[0] == staging.CURL_EXECUTABLE:
            destination = Path(argv[argv.index("--output") + 1])
            destination.write_bytes(self.downloads[argv[-1]])
            return b""
        if "imagetools" in argv:
            reference = argv[-1]
            if reference == self.image_coordinate:
                return self.index_manifest
            if reference == f"{self.image_coordinate}@{self.image_digest}":
                return self.selected_manifest
            raise AssertionError(f"unexpected raw manifest reference: {reference}")
        if argv[-3:-1] == ("image", "pull") or "pull" in argv:
            return b""
        if argv[-3:-1] == ("image", "inspect") or "inspect" in argv:
            return _json_bytes(
                [
                    {
                        "Id": self.image_config_digest,
                        "Os": "linux",
                        "Architecture": "amd64",
                        "RootFS": {"Type": "layers", "Layers": []},
                    }
                ]
            )
        if "save" in argv:
            destination = Path(argv[argv.index("--output") + 1])
            destination.write_bytes(self.docker_archive)
            return b""
        raise AssertionError(f"unexpected command: {argv}")


class _AdmissionRecorder:
    def __init__(self, *, fail_call: int | None = None) -> None:
        self.roots: list[Path] = []
        self.fail_call = fail_call

    def __call__(self, binding) -> None:
        self.roots.append(binding.root)
        if self.fail_call == len(self.roots):
            raise staging.StagingError("runtime_admission_failed")
        assert binding.root.is_dir()
        assert (binding.root / binding.container_lock).is_file()


@pytest.fixture
def staging_case(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    commit = staging.NFCORE_RNASEQ_COMMIT
    source_files = {
        "main.nf": (b"nextflow.enable.dsl=2\n", 0o644),
        "bin/tool": (b"#!/bin/sh\nexit 0\n", 0o755),
    }
    source_entries = [
        {
            "path": path,
            "size_bytes": len(content),
            "sha256": _sha256(content),
            "executable": bool(mode & 0o111),
        }
        for path, (content, mode) in sorted(source_files.items())
    ]
    source_archive = _tar_bytes(
        {f"rnaseq-{commit}/{path}": value for path, value in source_files.items()},
        compressed=True,
    )
    source_manifest = {
        "schema_version": "1.0.0",
        "project": "nf-core/rnaseq",
        "release": "3.26.0",
        "commit": commit,
        "source_archive_sha256": _sha256(source_archive),
        "file_count": len(source_entries),
        "total_size_bytes": sum(item["size_bytes"] for item in source_entries),
        "tree_sha256": _tree_digest(source_entries),
        "files": source_entries,
    }

    nextflow = b"#!/bin/sh\necho nextflow\n"
    jdk_files = {
        "LICENSE": (b"license\n", 0o644),
        "ASSEMBLY_EXCEPTION": (b"exception\n", 0o644),
        "bin/java": (b"#!/bin/sh\necho java >&2\n", 0o755),
    }
    jdk_entries = [
        {
            "path": path,
            "size_bytes": len(content),
            "sha256": _sha256(content),
            "executable": bool(mode & 0o111),
        }
        for path, (content, mode) in sorted(jdk_files.items())
    ]
    jdk_archive = _tar_bytes(
        {f"corretto/{path}": value for path, value in jdk_files.items()},
        compressed=True,
    )

    plugin_files = {"classes/Plugin.class": b"plugin-bytecode\n"}
    plugin_archive = _zip_bytes(plugin_files)
    plugin_meta = _json_bytes(
        {
            "version": "2.5.1",
            "requires": ">=25.04.0",
            "sha512sum": hashlib.sha512(plugin_archive).hexdigest(),
        }
    )
    plugin_entries = [
        {
            "path": path,
            "size_bytes": len(content),
            "sha256": _sha256(content),
        }
        for path, content in sorted(
            {
                **plugin_files,
                "nf-schema-2.5.1-meta.json": plugin_meta,
            }.items()
        )
    ]

    image_coordinate = "quay.io/helixweave/tiny:1.0"
    image_config = _json_bytes(
        {
            "architecture": "amd64",
            "os": "linux",
            "rootfs": {"type": "layers", "diff_ids": []},
        }
    )
    image_config_digest = f"sha256:{_sha256(image_config)}"
    selected_manifest = _json_bytes(
        {
            "schemaVersion": 2,
            "mediaType": "application/vnd.oci.image.manifest.v1+json",
            "config": {
                "mediaType": "application/vnd.oci.image.config.v1+json",
                "size": len(image_config),
                "digest": image_config_digest,
            },
            "layers": [],
        }
    )
    image_digest = f"sha256:{_sha256(selected_manifest)}"
    index_manifest = _json_bytes(
        {
            "schemaVersion": 2,
            "mediaType": "application/vnd.oci.image.index.v1+json",
            "manifests": [
                {
                    "mediaType": "application/vnd.oci.image.manifest.v1+json",
                    "size": len(selected_manifest),
                    "digest": image_digest,
                    "platform": {"os": "linux", "architecture": "amd64"},
                }
            ],
        }
    )
    docker_archive = _tar_bytes(
        {
            f"{image_config_digest.removeprefix('sha256:')}.json": (
                image_config,
                0o644,
            ),
            "manifest.json": (
                _json_bytes(
                    [
                        {
                            "Config": (
                                image_config_digest.removeprefix("sha256:") + ".json"
                            ),
                            "RepoTags": [image_coordinate],
                            "Layers": [],
                        }
                    ]
                ),
                0o644,
            ),
        },
        compressed=False,
    )

    inventory_entries = [
        {
            "process": process,
            "image_coordinate": image_coordinate,
            "source_file": "main.nf",
            "source_file_sha256": _sha256(source_files["main.nf"][0]),
        }
        for process in ("SALMON_QUANT", "STAR_ALIGN")
    ]
    inventory_sha = _sha256(
        json.dumps(
            inventory_entries,
            sort_keys=True,
            separators=(",", ":"),
        ).encode()
    )
    inventory = {
        "schema_version": "1.0.0",
        "project": "nf-core/rnaseq",
        "release": "3.26.0",
        "commit": commit,
        "process_count": len(inventory_entries),
        "entries_sha256": inventory_sha,
        "entries": inventory_entries,
    }

    contract_root = (
        Path(__file__).resolve().parents[2]
        / "src/encode_pipeline/contracts/nfcore_rnaseq"
    )
    container_schema = json.loads(
        (contract_root / "container-availability-lock-1.0.0.schema.json").read_bytes()
    )
    schema_pipeline = container_schema["properties"]["pipeline"]["properties"]
    schema_pipeline["source_tree_sha256"]["const"] = source_manifest["tree_sha256"]
    schema_pipeline["container_inventory_sha256"]["const"] = inventory_sha

    identity = {
        "schema_version": "1.0.0",
        "source": {
            "project": "nf-core/rnaseq",
            "release": "3.26.0",
            "commit": commit,
            "default_relative_path": f"source/nf-core-rnaseq-{commit}",
            "tree_sha256": source_manifest["tree_sha256"],
        },
        "nextflow": {
            "version": "25.04.3",
            "build": "5949",
            "default_relative_path": "nextflow/nextflow-25.04.3-dist",
            "size_bytes": len(nextflow),
            "sha256": _sha256(nextflow),
        },
        "jdk": {
            "distribution": "amazon-corretto",
            "version": "21.0.7.6.1",
            "runtime_version": "21.0.7+6-LTS",
            "archive_url": "https://fixtures.invalid/corretto.tar.gz",
            "default_archive_relative_path": "jdk/corretto.tar.gz",
            "archive_size_bytes": len(jdk_archive),
            "archive_sha256": _sha256(jdk_archive),
            "default_tree_relative_path": "jdk/corretto",
            "tree_file_count": len(jdk_entries),
            "tree_size_bytes": sum(item["size_bytes"] for item in jdk_entries),
            "tree_sha256": _tree_digest(jdk_entries),
            "java_relative_path": "bin/java",
            "java_size_bytes": len(jdk_files["bin/java"][0]),
            "java_sha256": _sha256(jdk_files["bin/java"][0]),
            "license_file": "LICENSE",
            "license_file_sha256": _sha256(jdk_files["LICENSE"][0]),
            "assembly_exception_file": "ASSEMBLY_EXCEPTION",
            "assembly_exception_file_sha256": _sha256(
                jdk_files["ASSEMBLY_EXCEPTION"][0]
            ),
        },
        "plugins": [
            {
                "id": "nf-schema",
                "version": "2.5.1",
                "default_archive_relative_path": "plugins/nf-schema-2.5.1.zip",
                "archive_size_bytes": len(plugin_archive),
                "archive_sha256": _sha256(plugin_archive),
                "default_meta_relative_path": ("plugins/nf-schema-2.5.1-meta.json"),
                "meta_size_bytes": len(plugin_meta),
                "meta_sha256": _sha256(plugin_meta),
                "default_tree_relative_path": "plugins/nf-schema-2.5.1",
                "tree_file_count": len(plugin_entries),
                "tree_sha256": _tree_digest(plugin_entries),
            }
        ],
        "containers": {"asset_format": "docker-archive+distribution-manifest"},
    }
    contract = staging.StagingContract(
        identity=identity,
        source_manifest=source_manifest,
        container_inventory=inventory,
        container_schema=container_schema,
    )
    downloads = {
        staging.SOURCE_ARCHIVE_URL: source_archive,
        staging.NEXTFLOW_DISTRIBUTION_URL: nextflow,
        staging.CORRETTO_ARCHIVE_URL: jdk_archive,
        staging.NF_SCHEMA_ARCHIVE_URL: plugin_archive,
        staging.NF_SCHEMA_META_URL: plugin_meta,
    }
    runner = _FakeRunner(
        downloads=downloads,
        image_coordinate=image_coordinate,
        index_manifest=index_manifest,
        selected_manifest=selected_manifest,
        image_digest=image_digest,
        image_config_digest=image_config_digest,
        docker_archive=docker_archive,
    )
    executable = tmp_path / "bin/docker"
    executable.parent.mkdir()
    executable.write_bytes(b"#!/bin/sh\nexit 1\n")
    executable.chmod(0o755)
    socket_path = tmp_path / "docker.sock"
    socket_path.write_bytes(b"not-a-real-socket\n")
    endpoint_validator = staging._verify_docker_endpoint
    monkeypatch.setattr(staging, "_verify_docker_endpoint", lambda _options: None)
    return {
        "contract": contract,
        "runner": runner,
        "asset_root": tmp_path / "runtime-assets",
        "docker_executable": executable,
        "docker_socket": socket_path,
        "endpoint_validator": endpoint_validator,
        "image_coordinate": image_coordinate,
        "image_digest": image_digest,
        "selected_manifest": selected_manifest,
        "container_schema": container_schema,
    }


def _options(case, **updates) -> staging.StageOptions:
    options = staging.StageOptions(
        asset_root=case["asset_root"],
        docker_executable=case["docker_executable"],
        docker_socket=case["docker_socket"],
    )
    return replace(options, **updates)


def test_committed_download_coordinates_are_fixed_to_verified_assets():
    contract = staging.load_committed_contract()
    jdk = contract.identity["jdk"]
    [plugin] = contract.identity["plugins"]

    assert staging.CORRETTO_ARCHIVE_URL == (
        "https://corretto.aws/downloads/resources/21.0.7.6.1/"
        "amazon-corretto-21.0.7.6.1-linux-x64.tar.gz"
    )
    assert jdk["version"] == "21.0.7.6.1"
    assert jdk["archive_size_bytes"] == 208603382
    assert jdk["archive_sha256"] == (
        "8bb627728d147e7507b2e38a5ef872549e895da50c2685d435c0d4c15ba95eb4"
    )
    assert staging.NF_SCHEMA_ARCHIVE_URL.endswith("/nf-schema-2.5.1.zip")
    assert staging.NF_SCHEMA_META_URL.endswith("/nf-schema-2.5.1-meta.json")
    assert plugin["archive_size_bytes"] == 478526
    assert plugin["archive_sha256"] == (
        "f3833d0c29a51dc5e3759e00a6af87fcb1e989c94d1c1e7995d06e3f7090c461"
    )
    assert plugin["meta_size_bytes"] == 354
    assert plugin["meta_sha256"] == (
        "acf2f2deb1a71940e5d27255508b3a31aca221638544b67edf81790235c491b5"
    )


def test_external_runner_uses_shell_false_and_a_secret_free_environment(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    command_home = tmp_path / "command-home"
    (command_home / "docker").mkdir(parents=True)
    observed = {}

    def fake_run(argv, **kwargs):
        observed["argv"] = argv
        observed.update(kwargs)
        return staging.subprocess.CompletedProcess(argv, 0, stdout=b"manifest")

    monkeypatch.setenv("HELIXWEAVE_TEST_SECRET", "must-not-leak")
    monkeypatch.setattr(staging.subprocess, "run", fake_run)

    output = staging._run_external(
        ("/usr/bin/docker", "version"),
        capture_stdout=True,
        command_home=command_home,
    )

    assert output == b"manifest"
    assert observed["shell"] is False
    assert observed["stdin"] is staging.subprocess.DEVNULL
    assert observed["stderr"] is staging.subprocess.DEVNULL
    assert "HELIXWEAVE_TEST_SECRET" not in observed["env"]
    assert set(observed["env"]) == {"DOCKER_CONFIG", "HOME", "LANG", "LC_ALL", "PATH"}


def test_endpoint_validation_rejects_a_non_socket(staging_case):
    with pytest.raises(staging.StagingError, match="docker_endpoint_invalid"):
        staging_case["endpoint_validator"](_options(staging_case))


def test_plan_is_the_default_and_never_invokes_network_or_mutation(staging_case):
    recorder = _AdmissionRecorder()

    report = staging.stage_runtime_assets(
        _options(staging_case),
        contract=staging_case["contract"],
        command_runner=staging_case["runner"],
        admission_verifier=recorder,
    )

    assert report.phase == "plan"
    assert report.container_process_count == 2
    assert report.unique_image_count == 1
    assert staging_case["runner"].calls == []
    assert recorder.roots == []
    assert not staging_case["asset_root"].exists()


@pytest.mark.parametrize(
    ("allow_network", "allow_mutation"),
    [(False, False), (True, False), (False, True)],
)
def test_stage_requires_separate_explicit_network_and_mutation_authority(
    staging_case,
    allow_network: bool,
    allow_mutation: bool,
):
    with pytest.raises(staging.StagingError, match="stage_authority_required"):
        staging.stage_runtime_assets(
            _options(
                staging_case,
                phase="stage",
                allow_network=allow_network,
                allow_mutation=allow_mutation,
            ),
            contract=staging_case["contract"],
            command_runner=staging_case["runner"],
            admission_verifier=_AdmissionRecorder(),
        )

    assert staging_case["runner"].calls == []
    assert not staging_case["asset_root"].exists()


def test_stage_builds_one_digest_closed_image_and_schema_valid_lock(staging_case):
    recorder = _AdmissionRecorder()

    report = staging.stage_runtime_assets(
        _options(
            staging_case,
            phase="stage",
            allow_network=True,
            allow_mutation=True,
        ),
        contract=staging_case["contract"],
        command_runner=staging_case["runner"],
        admission_verifier=recorder,
    )

    root = staging_case["asset_root"]
    assert root.is_dir()
    assert report.phase == "stage"
    assert report.container_process_count == 2
    assert report.unique_image_count == 1
    assert len(recorder.roots) == 2
    assert recorder.roots[0] != root
    assert recorder.roots[1] == root
    assert not any(root.parent.glob(f".{root.name}.stage-*"))

    lock = json.loads((root / "containers/availability-lock.json").read_bytes())
    assert (
        list(Draft202012Validator(staging_case["container_schema"]).iter_errors(lock))
        == []
    )
    assert [entry["process"] for entry in lock["entries"]] == [
        "SALMON_QUANT",
        "STAR_ALIGN",
    ]
    assert {entry["oci_digest"] for entry in lock["entries"]} == {
        staging_case["image_digest"]
    }
    assert len({entry["local_asset"] for entry in lock["entries"]}) == 1
    assert len({entry["distribution_manifest_asset"] for entry in lock["entries"]}) == 1
    [manifest_name] = {
        entry["distribution_manifest_asset"] for entry in lock["entries"]
    }
    assert (root / "containers/assets" / manifest_name).read_bytes() == staging_case[
        "selected_manifest"
    ]
    assert (
        stat.S_IMODE((root / "nextflow/nextflow-25.04.3-dist").stat().st_mode) == 0o755
    )

    calls = [call for call, _capture in staging_case["runner"].calls]
    assert (
        sum(
            "imagetools" in call and call[-1] == staging_case["image_coordinate"]
            for call in calls
        )
        == 1
    )
    assert sum("pull" in call for call in calls) == 1
    assert sum("save" in call for call in calls) == 1


def test_wrong_platform_manifest_fails_before_pull_and_cleans_candidate(staging_case):
    raw_index = json.loads(staging_case["runner"].index_manifest)
    raw_index["manifests"][0]["platform"] = {
        "os": "linux",
        "architecture": "arm64",
    }
    staging_case["runner"].index_manifest = _json_bytes(raw_index)

    with pytest.raises(staging.StagingError, match="container_platform_invalid"):
        staging.stage_runtime_assets(
            _options(
                staging_case,
                phase="stage",
                allow_network=True,
                allow_mutation=True,
            ),
            contract=staging_case["contract"],
            command_runner=staging_case["runner"],
            admission_verifier=_AdmissionRecorder(),
        )

    assert not staging_case["asset_root"].exists()
    assert not any(
        staging_case["asset_root"].parent.glob(
            f".{staging_case['asset_root'].name}.stage-*"
        )
    )
    calls = [call for call, _capture in staging_case["runner"].calls]
    assert not any("pull" in call or "save" in call for call in calls)


def test_download_identity_failure_is_atomic_and_path_redacted(staging_case):
    staging_case["runner"].downloads[staging.NEXTFLOW_DISTRIBUTION_URL] += b"tamper"

    with pytest.raises(staging.StagingError) as error:
        staging.stage_runtime_assets(
            _options(
                staging_case,
                phase="stage",
                allow_network=True,
                allow_mutation=True,
            ),
            contract=staging_case["contract"],
            command_runner=staging_case["runner"],
            admission_verifier=_AdmissionRecorder(),
        )

    assert error.value.code == "asset_identity_invalid"
    assert str(staging_case["asset_root"]) not in str(error.value)
    assert not staging_case["asset_root"].exists()
    assert not any(
        staging_case["asset_root"].parent.glob(
            f".{staging_case['asset_root'].name}.stage-*"
        )
    )


def test_existing_root_is_never_overwritten(staging_case):
    root = staging_case["asset_root"]
    root.mkdir()
    sentinel = root / "user-owned"
    sentinel.write_bytes(b"preserve\n")

    with pytest.raises(staging.StagingError, match="asset_root_exists"):
        staging.stage_runtime_assets(
            _options(
                staging_case,
                phase="stage",
                allow_network=True,
                allow_mutation=True,
            ),
            contract=staging_case["contract"],
            command_runner=staging_case["runner"],
            admission_verifier=_AdmissionRecorder(),
        )

    assert sentinel.read_bytes() == b"preserve\n"
    assert staging_case["runner"].calls == []


def test_atomic_publish_never_replaces_a_concurrently_created_root(tmp_path: Path):
    candidate = tmp_path / "candidate"
    candidate.mkdir()
    (candidate / "staged").write_bytes(b"candidate\n")
    existing = tmp_path / "existing"
    existing.mkdir()
    sentinel = existing / "user-owned"
    sentinel.write_bytes(b"preserve\n")

    with pytest.raises(staging.StagingError, match="asset_root_exists"):
        staging._publish_candidate(candidate, existing)

    assert (candidate / "staged").read_bytes() == b"candidate\n"
    assert sentinel.read_bytes() == b"preserve\n"


def test_post_publish_admission_failure_removes_only_new_root(staging_case):
    recorder = _AdmissionRecorder(fail_call=2)

    with pytest.raises(staging.StagingError, match="runtime_admission_failed"):
        staging.stage_runtime_assets(
            _options(
                staging_case,
                phase="stage",
                allow_network=True,
                allow_mutation=True,
            ),
            contract=staging_case["contract"],
            command_runner=staging_case["runner"],
            admission_verifier=recorder,
        )

    assert len(recorder.roots) == 2
    assert not staging_case["asset_root"].exists()


def test_verify_is_read_only_and_uses_final_root_admission(staging_case):
    staging.stage_runtime_assets(
        _options(
            staging_case,
            phase="stage",
            allow_network=True,
            allow_mutation=True,
        ),
        contract=staging_case["contract"],
        command_runner=staging_case["runner"],
        admission_verifier=_AdmissionRecorder(),
    )
    staging_case["runner"].calls.clear()
    before = {
        path.relative_to(staging_case["asset_root"]).as_posix(): (
            path.stat().st_size,
            path.stat().st_mtime_ns,
        )
        for path in staging_case["asset_root"].rglob("*")
    }
    recorder = _AdmissionRecorder()

    report = staging.stage_runtime_assets(
        _options(staging_case, phase="verify"),
        contract=staging_case["contract"],
        command_runner=staging_case["runner"],
        admission_verifier=recorder,
    )

    after = {
        path.relative_to(staging_case["asset_root"]).as_posix(): (
            path.stat().st_size,
            path.stat().st_mtime_ns,
        )
        for path in staging_case["asset_root"].rglob("*")
    }
    assert report.phase == "verify"
    assert recorder.roots == [staging_case["asset_root"]]
    assert staging_case["runner"].calls == []
    assert after == before


@pytest.mark.parametrize("unsafe", ["relative", "symlink_parent"])
def test_asset_root_path_is_fail_closed(tmp_path: Path, staging_case, unsafe: str):
    if unsafe == "relative":
        asset_root = Path("relative-assets")
    else:
        real_parent = tmp_path / "real-parent"
        real_parent.mkdir()
        linked_parent = tmp_path / "linked-parent"
        linked_parent.symlink_to(real_parent, target_is_directory=True)
        asset_root = linked_parent / "assets"

    with pytest.raises(staging.StagingError, match="asset_root_invalid"):
        staging.stage_runtime_assets(
            replace(_options(staging_case), asset_root=asset_root),
            contract=staging_case["contract"],
            command_runner=staging_case["runner"],
            admission_verifier=_AdmissionRecorder(),
        )

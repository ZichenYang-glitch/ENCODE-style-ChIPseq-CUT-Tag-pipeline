"""Administrator-only Project and Sample registry CLI contracts."""

from __future__ import annotations

from contextlib import contextmanager
import json
from pathlib import Path
from typing import Iterator

import pytest

from encode_pipeline.cli import admin


DATABASE_URL = "sqlite:////tmp/helixweave-admin-test.db"


class RecordingRegistry:
    def __init__(self) -> None:
        self.calls: list[tuple[str, object]] = []

    def create_project(self, *, display_name: str) -> dict[str, object]:
        self.calls.append(("create_project", display_name))
        return {
            "project_id": "prj_11111111111111111111111111111111",
            "display_name": display_name,
            "kind": "user",
            "created_at": "2026-07-26T08:00:00Z",
            "archived_at": None,
        }

    def list_projects(self, *, include_archived: bool) -> tuple[dict[str, object], ...]:
        self.calls.append(("list_projects", include_archived))
        return (
            {
                "project_id": "prj_11111111111111111111111111111111",
                "display_name": "Project One",
            },
        )

    def archive_project(self, *, project_id: str) -> dict[str, object]:
        self.calls.append(("archive_project", project_id))
        return {
            "project_id": project_id,
            "archived_at": "2026-07-26T09:00:00Z",
        }

    def import_samples(
        self,
        *,
        project_id: str,
        rows: tuple[admin.SampleImportRow, ...],
    ) -> tuple[dict[str, object], ...]:
        self.calls.append(("import_samples", (project_id, rows)))
        return tuple(
            {
                "sample_id": f"smp_{index:032x}",
                "project_id": project_id,
                "sample_key": row.sample_key,
                "display_name": row.display_name,
                "attributes": row.attributes,
            }
            for index, row in enumerate(rows, start=1)
        )

    def list_samples(self, *, project_id: str) -> tuple[dict[str, object], ...]:
        self.calls.append(("list_samples", project_id))
        return (
            {
                "sample_id": "smp_11111111111111111111111111111111",
                "project_id": project_id,
                "sample_key": "tumor-a",
            },
        )

    def revise_sample(
        self,
        *,
        sample_id: str,
        display_name: str,
        attributes: dict[str, object],
    ) -> dict[str, object]:
        self.calls.append(("revise_sample", (sample_id, display_name, attributes)))
        return {
            "sample_id": sample_id,
            "revision_id": "smpr_22222222222222222222222222222222",
            "revision": 2,
            "display_name": display_name,
            "attributes": attributes,
        }


def registry_factory(
    registry: RecordingRegistry,
    observed_database_urls: list[str],
):
    @contextmanager
    def factory(database_url: str) -> Iterator[RecordingRegistry]:
        observed_database_urls.append(database_url)
        yield registry

    return factory


def invoke(
    arguments: list[str],
    *,
    registry: RecordingRegistry,
    observed_database_urls: list[str],
) -> int:
    return admin.main(
        ["--database-url", DATABASE_URL, *arguments],
        registry_factory=registry_factory(registry, observed_database_urls),
    )


def test_project_create_list_and_archive_emit_machine_readable_json(capsys) -> None:
    registry = RecordingRegistry()
    database_urls: list[str] = []

    assert (
        invoke(
            ["project", "create", "--display-name", "Project One"],
            registry=registry,
            observed_database_urls=database_urls,
        )
        == 0
    )
    created = json.loads(capsys.readouterr().out)
    assert created["project_id"] == "prj_11111111111111111111111111111111"
    assert created["display_name"] == "Project One"

    assert (
        invoke(
            ["project", "list", "--include-archived"],
            registry=registry,
            observed_database_urls=database_urls,
        )
        == 0
    )
    assert json.loads(capsys.readouterr().out) == [
        {
            "display_name": "Project One",
            "project_id": "prj_11111111111111111111111111111111",
        }
    ]

    project_id = "prj_11111111111111111111111111111111"
    assert (
        invoke(
            ["project", "archive", project_id],
            registry=registry,
            observed_database_urls=database_urls,
        )
        == 0
    )
    assert json.loads(capsys.readouterr().out)["archived_at"].endswith("Z")
    assert registry.calls == [
        ("create_project", "Project One"),
        ("list_projects", True),
        ("archive_project", project_id),
    ]
    assert database_urls == [DATABASE_URL, DATABASE_URL, DATABASE_URL]


def test_sample_import_reads_tsv_once_and_preserves_structured_attributes(
    tmp_path: Path,
    capsys,
) -> None:
    intake = tmp_path / "samples.tsv"
    intake.write_text(
        "sample_key\tdisplay_name\tattributes_json\n"
        'tumor-a\tTumor A\t{"condition":"treated","replicate":1}\n'
        "control-a\tControl A\t\n",
        encoding="utf-8",
    )
    registry = RecordingRegistry()
    database_urls: list[str] = []
    project_id = "prj_11111111111111111111111111111111"

    assert (
        invoke(
            [
                "sample",
                "import",
                "--project-id",
                project_id,
                "--tsv",
                str(intake),
            ],
            registry=registry,
            observed_database_urls=database_urls,
        )
        == 0
    )

    imported = json.loads(capsys.readouterr().out)
    assert [row["sample_key"] for row in imported] == ["tumor-a", "control-a"]
    assert imported[0]["attributes"] == {"condition": "treated", "replicate": 1}
    assert imported[1]["attributes"] == {}
    assert registry.calls == [
        (
            "import_samples",
            (
                project_id,
                (
                    admin.SampleImportRow(
                        sample_key="tumor-a",
                        display_name="Tumor A",
                        attributes={"condition": "treated", "replicate": 1},
                    ),
                    admin.SampleImportRow(
                        sample_key="control-a",
                        display_name="Control A",
                        attributes={},
                    ),
                ),
            ),
        )
    ]


@pytest.mark.parametrize(
    ("contents", "error"),
    [
        ("display_name\nTumor A\n", "sample_key"),
        (
            "sample_key\tdisplay_name\tattribute_json\n"
            'tumor-a\tTumor A\t{"condition":"treated"}\n',
            "unknown column(s): attribute_json",
        ),
        (
            "sample_key\tsample_key\tdisplay_name\ntumor-a\tshadowed-key\tTumor A\n",
            "duplicate column(s): sample_key",
        ),
        (
            "sample_key\t\tdisplay_name\ntumor-a\tunexpected\tTumor A\n",
            "empty column name",
        ),
        (
            "sample_key\tdisplay_name\ntumor-a\tTumor A\tunexpected trailing value\n",
            "unexpected trailing field(s)",
        ),
        (
            "sample_key\tdisplay_name\tattributes_json\ntumor-a\tTumor A\t[]\n",
            "JSON object",
        ),
        (
            "sample_key\tdisplay_name\tattributes_json\ntumor-a\tTumor A\t{not-json}\n",
            "valid JSON object",
        ),
    ],
)
def test_sample_import_rejects_invalid_tsv_before_opening_database(
    tmp_path: Path,
    capsys,
    contents: str,
    error: str,
) -> None:
    intake = tmp_path / "samples.tsv"
    intake.write_text(contents, encoding="utf-8")
    registry = RecordingRegistry()
    database_urls: list[str] = []

    with pytest.raises(SystemExit, match="2"):
        invoke(
            [
                "sample",
                "import",
                "--project-id",
                "prj_11111111111111111111111111111111",
                "--tsv",
                str(intake),
            ],
            registry=registry,
            observed_database_urls=database_urls,
        )

    assert error in capsys.readouterr().err
    assert database_urls == []
    assert registry.calls == []


def test_sample_list_and_revise_use_stable_ids_and_full_revision_payload(
    capsys,
) -> None:
    registry = RecordingRegistry()
    database_urls: list[str] = []
    project_id = "prj_11111111111111111111111111111111"
    sample_id = "smp_11111111111111111111111111111111"

    assert (
        invoke(
            ["sample", "list", "--project-id", project_id],
            registry=registry,
            observed_database_urls=database_urls,
        )
        == 0
    )
    assert json.loads(capsys.readouterr().out)[0]["sample_id"] == sample_id

    assert (
        invoke(
            [
                "sample",
                "revise",
                sample_id,
                "--display-name",
                "Tumor A revised",
                "--attributes-json",
                '{"condition":"treated","replicate":2}',
            ],
            registry=registry,
            observed_database_urls=database_urls,
        )
        == 0
    )
    revised = json.loads(capsys.readouterr().out)
    assert revised["revision"] == 2
    assert registry.calls == [
        ("list_samples", project_id),
        (
            "revise_sample",
            (
                sample_id,
                "Tumor A revised",
                {"condition": "treated", "replicate": 2},
            ),
        ),
    ]


def test_admin_requires_explicit_database_and_has_no_server_path_options(
    capsys,
) -> None:
    with pytest.raises(SystemExit, match="2"):
        admin.main(["project", "list"])
    assert "--database-url" in capsys.readouterr().err

    with pytest.raises(SystemExit, match="0"):
        admin.main(["--help"])
    help_text = capsys.readouterr().out
    assert "--database-url" in help_text
    assert "--server-path" not in help_text
    assert "--workspace-root" not in help_text
    assert "--storage-root" not in help_text


def test_real_sqlite_registry_cli_round_trip_is_revision_preserving(
    tmp_path: Path,
    capsys,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"

    assert (
        admin.main(
            [
                "--database-url",
                database_url,
                "project",
                "create",
                "--display-name",
                "Pilot project",
            ]
        )
        == 0
    )
    project = json.loads(capsys.readouterr().out)
    assert project["kind"] == "user"
    project_id = project["project_id"]

    intake = tmp_path / "samples.tsv"
    intake.write_text(
        "sample_key\tdisplay_name\tattributes_json\n"
        'donor-01\tDonor 01\t{"condition":"initial"}\n',
        encoding="utf-8",
    )
    assert (
        admin.main(
            [
                "--database-url",
                database_url,
                "sample",
                "import",
                "--project-id",
                project_id,
                "--tsv",
                str(intake),
            ]
        )
        == 0
    )
    imported = json.loads(capsys.readouterr().out)
    sample_id = imported[0]["sample"]["sample_id"]
    original_revision_id = imported[0]["revision"]["sample_revision_id"]
    assert imported[0]["sample"]["stable_key"] == "donor-01"

    assert (
        admin.main(
            [
                "--database-url",
                database_url,
                "sample",
                "revise",
                sample_id,
                "--display-name",
                "Donor 01 reviewed",
                "--attributes-json",
                '{"condition":"reviewed"}',
            ]
        )
        == 0
    )
    revised = json.loads(capsys.readouterr().out)
    assert revised["sample_revision_id"] != original_revision_id
    assert revised["revision_number"] == 2

    assert (
        admin.main(
            [
                "--database-url",
                database_url,
                "sample",
                "list",
                "--project-id",
                project_id,
            ]
        )
        == 0
    )
    listed = json.loads(capsys.readouterr().out)
    assert listed[0]["sample"]["sample_id"] == sample_id
    assert listed[0]["revision"]["sample_revision_id"] == revised["sample_revision_id"]
    assert json.loads(listed[0]["revision"]["canonical_payload"]) == {
        "attributes": {"condition": "reviewed"},
        "display_name": "Donor 01 reviewed",
    }


def test_real_sqlite_batch_import_rolls_back_every_row_on_late_conflict(
    tmp_path: Path,
    capsys,
) -> None:
    database_url = f"sqlite:///{tmp_path / 'platform.db'}"
    assert (
        admin.main(
            [
                "--database-url",
                database_url,
                "project",
                "create",
                "--display-name",
                "Pilot project",
            ]
        )
        == 0
    )
    project_id = json.loads(capsys.readouterr().out)["project_id"]
    existing = tmp_path / "existing.tsv"
    existing.write_text(
        "sample_key\tdisplay_name\nexisting\tExisting sample\n",
        encoding="utf-8",
    )
    assert (
        admin.main(
            [
                "--database-url",
                database_url,
                "sample",
                "import",
                "--project-id",
                project_id,
                "--tsv",
                str(existing),
            ]
        )
        == 0
    )
    capsys.readouterr()

    conflicting = tmp_path / "conflicting.tsv"
    conflicting.write_text(
        "sample_key\tdisplay_name\n"
        "would-partially-write\tMust roll back\n"
        "existing\tConflicts after first row\n",
        encoding="utf-8",
    )
    assert (
        admin.main(
            [
                "--database-url",
                database_url,
                "sample",
                "import",
                "--project-id",
                project_id,
                "--tsv",
                str(conflicting),
            ]
        )
        == 1
    )
    error = json.loads(capsys.readouterr().err)
    assert error["error"]["code"] == "registry-command-failed"

    assert (
        admin.main(
            [
                "--database-url",
                database_url,
                "sample",
                "list",
                "--project-id",
                project_id,
            ]
        )
        == 0
    )
    listed = json.loads(capsys.readouterr().out)
    assert [entry["sample"]["stable_key"] for entry in listed] == ["existing"]

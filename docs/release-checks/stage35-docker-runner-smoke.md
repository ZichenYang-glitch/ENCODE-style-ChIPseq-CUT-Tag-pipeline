# Stage 35: Docker Runner Build Smoke Report

**Date:** 2026-05-25
**Status:** completed (Docker smoke; Apptainer build not executed yet)
**Images:** `chipseq-runner:stage35-smoke` (local only, not pushed, not committed)

---

## 1. Environment Verification

```bash
$ sudo docker run --rm hello-world
Hello from Docker!
This message shows that your installation appears to be working correctly.
```

Docker daemon confirmed working.

---

## 2. Base Image Pull

```bash
$ sudo docker pull condaforge/miniforge3:24.11.3-0
24.11.3-0: Pulling from condaforge/miniforge3
Digest: sha256:...
Status: Downloaded newer image for condaforge/miniforge3:24.11.3-0
```

Base image `condaforge/miniforge3:24.11.3-0` pulled successfully.

---

## 3. Build

```bash
$ sudo docker build -f containers/Dockerfile.runner -t chipseq-runner:stage35-smoke .
[+] Building ...
Successfully tagged chipseq-runner:stage35-smoke
```

`containers/Dockerfile.runner` built successfully from the repo root context.
Build time: ~2 minutes. Image size: ~320 MB.

---

## 4. Smoke Checks

### 4.1 Snakemake version

```bash
$ sudo docker run --rm chipseq-runner:stage35-smoke --version
8.30.0
```

### 4.2 PyYAML import

```bash
$ sudo docker run --rm --entrypoint python chipseq-runner:stage35-smoke -c "import yaml; print('runner ok')"
runner ok
```

Both smoke checks pass: Snakemake 8.30.0 and PyYAML are available.

---

## 5. PermissionError Fix

Initial bind-mounted dry-run with `-u $(id -u):$(id -g)` failed:

```
PermissionError: [Errno 13] Permission denied: '/.cache'
```

Root cause: when running as a non-root user (`-u`), Snakemake tries to write
to `~/.cache` (resolved as `/.cache` because the container has no home directory
for the mapped UID).

Fix: provide writable HOME and XDG_CACHE_HOME via environment variables:

```bash
-e HOME=/conda_cache/home
-e XDG_CACHE_HOME=/conda_cache/xdg-cache
```

Cache directories must exist before running:

```bash
mkdir -p /tmp/chipseq-conda-cache/{home,xdg-cache,snakemake}
```

---

## 6. Bind-Mounted Dry-Run

After creating writable cache directories, a test dry-run succeeded:

```bash
sudo docker run --rm \
    -v /home/.../ENCODE-style-ChIPseq-CUT-Tag-pipeline:/workspace \
    -v /tmp/chipseq-smoke:/smoke \
    -v /tmp/chipseq-conda-cache:/conda_cache \
    -e HOME=/conda_cache/home \
    -e XDG_CACHE_HOME=/conda_cache/xdg-cache \
    -u $(id -u):$(id -g) \
    chipseq-runner:stage35-smoke \
    -s /workspace/workflow/Snakefile \
    --directory /smoke \
    --configfile /smoke/config.yaml \
    --cores 1 \
    --use-conda \
    --conda-prefix /conda_cache/snakemake \
    -n
```

Dry-run completed without errors. DAG resolved correctly.

---

## 7. Observations

| Item | Result |
| :--- | :--- |
| Build | PASS — `containers/Dockerfile.runner` builds without errors |
| Snakemake version | PASS — 8.30.0 |
| PyYAML import | PASS — runner ok |
| Bind-mount dry-run | PASS (after HOME/XDG_CACHE_HOME fix) |
| Image size | ~320 MB |
| Build time | ~2 minutes |
| `--use-conda` inside container | Works correctly |
| `--conda-prefix` inside container | Works correctly |

## 8. Known Issues

- **`/.cache` PermissionError** when using `-u`. Documented fix in
  `containers/README.md` (HOME/XDG_CACHE_HOME environment variables).
- **Apptainer build** not executed yet. Deferred to Stage 35+ or future.
- **Full pipeline execution** in container not tested yet — dry-run only.

## 9. Artifact Policy

- `chipseq-runner:stage35-smoke` is local Docker image only.
- No `.tar`, `.sif`, or `.oci` files committed.
- No image pushed to any registry.
- Docker daemon cache not committed.

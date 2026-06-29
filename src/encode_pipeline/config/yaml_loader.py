"""YAML config loading with PyYAML and stdlib fallback."""

import os


def load_yaml(path: str) -> dict:
    """Load YAML config, preferring PyYAML if available, with stdlib fallback.

    PyYAML is available in the chipseq Conda environment (transitive
    dependency of Snakemake). The stdlib fallback handles bare-minimum
    YAML parsing for standalone use outside Conda.
    """
    try:
        import yaml
    except ImportError:
        return parse_config_minimal(path)

    with open(path) as fh:
        return yaml.safe_load(fh)


def parse_config_minimal(path: str) -> dict:
    """Minimal YAML config parser using stdlib only.

    Extracts flat keys plus the genome_resources and qc nested blocks.
    This covers the config structure used by validate_config().
    """
    config: dict = {}
    genome_resources: dict = {}
    qc: dict = {}
    tool_parameters: dict = {}
    cuttag: dict = {}
    seacr_sub: dict = {}
    section: str | None = None
    current_genome: str | None = None
    current_entry: dict | None = None
    current_tool: str | None = None
    current_tool_entry: dict | None = None

    def parse_scalar(value: str):
        value = value.strip().strip('"').strip("'")
        if value == "":
            return ""
        if value.lower() in ("true", "false"):
            return value.lower() == "true"
        if value.isdigit():
            return int(value)
        return value

    def save_current_genome() -> None:
        nonlocal current_genome, current_entry
        if current_genome is not None and current_entry is not None:
            genome_resources[current_genome] = current_entry
        current_genome = None
        current_entry = None

    def save_current_tool() -> None:
        nonlocal current_tool, current_tool_entry
        if current_tool is not None and current_tool_entry is not None:
            tool_parameters[current_tool] = current_tool_entry
        current_tool = None
        current_tool_entry = None

    with open(path) as fh:
        for line in fh:
            stripped = line.rstrip()

            # Skip comments and empty lines
            if not stripped or stripped.lstrip().startswith("#"):
                continue

            indent = len(line) - len(line.lstrip(" "))

            # Top-level keys start or end nested blocks.
            if indent == 0:
                save_current_genome()
                save_current_tool()

                if stripped.startswith("genome_resources:"):
                    section = "genome_resources"
                    continue
                if stripped.startswith("qc:"):
                    section = "qc"
                    continue
                if stripped.startswith("tool_parameters:"):
                    section = "tool_parameters"
                    continue
                if stripped.startswith("cuttag:"):
                    section = "cuttag"
                    continue

                section = None
                if ":" in stripped:
                    k, v = stripped.split(":", 1)
                    config[k.strip()] = parse_scalar(v)
                continue

            if section == "genome_resources" and indent == 2:
                # Genome key: "  hs:", "  mm10:", etc.
                key = stripped.strip().rstrip(":")
                if key and not key.startswith("#"):
                    save_current_genome()
                    current_genome = key
                    current_entry = {}
                continue

            if section == "genome_resources" and indent == 4:
                # Field within a genome entry
                field = stripped.strip()
                if ":" in field and current_genome is not None and current_entry is not None:
                    k, v = field.split(":", 1)
                    current_entry[k.strip()] = parse_scalar(v)
                continue

            if section == "qc" and indent == 2:
                field = stripped.strip()
                if ":" in field:
                    k, v = field.split(":", 1)
                    qc[k.strip()] = parse_scalar(v)
                continue

            if section == "cuttag" and indent == 2:
                key = stripped.strip().rstrip(":")
                if key == "seacr":
                    seacr_sub = {}
                elif ":" in stripped:
                    k, v = stripped.split(":", 1)
                    cuttag[k.strip()] = parse_scalar(v)
                continue

            if section == "cuttag" and indent == 4:
                field = stripped.strip()
                if ":" in field:
                    k, v = field.split(":", 1)
                    seacr_sub[k.strip()] = parse_scalar(v)
                continue

            if section == "tool_parameters" and indent == 2:
                # Tool block key: "  fastqc:", "  trim_galore:", etc.
                key = stripped.strip().rstrip(":")
                if key and not key.startswith("#"):
                    save_current_tool()
                    current_tool = key
                    current_tool_entry = {}
                continue

            if section == "tool_parameters" and indent == 4:
                # Field within a tool block
                field = stripped.strip()
                if ":" in field and current_tool is not None and current_tool_entry is not None:
                    k, v = field.split(":", 1)
                    current_tool_entry[k.strip()] = parse_scalar(v)
                continue

    # Save last genome entry
    save_current_genome()
    save_current_tool()

    if genome_resources:
        config["genome_resources"] = genome_resources
    if qc:
        config["qc"] = qc
    if tool_parameters:
        config["tool_parameters"] = tool_parameters
    if seacr_sub:
        cuttag["seacr"] = seacr_sub
    if cuttag:
        config["cuttag"] = cuttag

    return config

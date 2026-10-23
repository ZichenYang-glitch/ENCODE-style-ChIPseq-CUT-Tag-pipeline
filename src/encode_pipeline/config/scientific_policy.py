"""Internal scientific policy identity for the ENCODE-style workflow.

The policy revision identifies the shipped scientific behavior of the
workflow (filtering, duplicate, organelle, blacklist, and QC semantics).
It is an internal, non-user-overridable value: it is recorded in the run
identity marker and the effective configuration record so that corrected
workflow versions never mix outputs with directories produced by older
scientific behavior.
"""

# Bump when a shipped scientific policy change requires output isolation
# from runs produced by earlier behavior. 1 = SCI-1 run-identity mechanism;
# 2 = SCI-2 assay-based duplicate policy + MACS3 --keep-dup all.
SCIENTIFIC_POLICY_REVISION = 2

# Schema version of the .helixweave/analysis_identity.json marker.
# Version 1 markers were the unversioned {"identity", "analysis_profile"}
# payloads written before the marker format was versioned; they are
# rejected as outdated rather than reinterpreted.
ANALYSIS_MARKER_SCHEMA_VERSION = 2

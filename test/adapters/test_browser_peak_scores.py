"""WashU descriptions distinguish IDR rank from consensus support."""

from encode_pipeline.adapters.browser.catalog import describe_browser_artifact


WORKFLOW = "encode-style-chipseq-cuttag-atac-mnase"


def _describe(output_type, name):
    return describe_browser_artifact(
        workflow_id=WORKFLOW,
        output_type=output_type,
        name=name,
    )


def test_idr_peaks_expose_rank_including_chipseq_finals():
    for output_type, name in (
        ("idr_conservative", "conservative.narrowPeak"),
        ("idr_optimal", "optimal.narrowPeak"),
        (
            "atac_macs3_narrow_idr_final_peak",
            "EXP.atac.macs3.narrow.replicate_validated.idr.narrowPeak",
        ),
        (
            "cuttag_macs3_broad_idr_final_peak",
            "EXP.cuttag.macs3.broad.replicate_validated.idr.broadPeak",
        ),
    ):
        described = _describe(output_type, name)
        assert described is not None
        assert described.wash_u_type == "bed"
        assert described.role == "interval"
        assert described.score_meaning == "idr_rank"
        assert "IDR rank" in described.detail
        assert described.label == "IDR peaks"


def test_consensus_peaks_expose_biorep_support_including_new_finals():
    for output_type, name in (
        (
            "chipseq_macs3_narrow_consensus_peak",
            "EXP.chipseq.macs3.narrow.consensus.narrowPeak",
        ),
        (
            "chipseq_macs3_narrow_consensus_final_peak",
            "EXP.chipseq.macs3.narrow.replicate_validated.consensus.narrowPeak",
        ),
        (
            "atac_macs3_narrow_consensus_final_peak",
            "EXP.atac.macs3.narrow.replicate_validated.consensus.narrowPeak",
        ),
        (
            "cuttag_seacr_consensus_final_peak",
            "EXP.cuttag.seacr.stringent.replicate_validated.consensus.bed",
        ),
    ):
        described = _describe(output_type, name)
        assert described is not None
        assert described.wash_u_type == "bed"
        assert described.score_meaning == "consensus_biorep_support"
        assert "biological-replicate support" in described.detail
        assert described.label == "Consensus peaks"


def test_caller_peaks_keep_the_generic_score():
    described = _describe("macs3_peak", "S1_peaks.narrowPeak")
    assert described is not None
    assert described.score_meaning == "source_interval_score"
    assert described.label == "Peaks"


def test_summaries_are_not_tracks():
    assert (
        _describe(
            "chipseq_macs3_narrow_consensus_summary",
            "EXP.chipseq.macs3.narrow.consensus.summary.tsv",
        )
        is None
    )

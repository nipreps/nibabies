"""Tests for the subcortical MNI alignment workflow."""

from nibabies.workflows.bold.alignment import init_subcortical_mni_alignment_wf


def test_subcortical_alignment_wf_structure():
    wf = init_subcortical_mni_alignment_wf()
    names = set(wf.list_node_names())

    # per-ROI cropping, registration, and label merge are present
    assert 'roi_bbox' in names
    assert 'roi2atlas' in names
    assert 'merge_segs' in names
    # main QC reportlet always present; debug reportlet gated off by default
    assert 'subcortical_rpt' in names
    assert 'overlap_rpt' not in names

    # dead branches removed
    assert 'applyxfm_atlas' not in names
    assert 'vol_resample' not in names
    assert 'agg_rois' not in names
    assert 'fmt_agg_rois' not in names

    out_fields = wf.get_node('outputnode').inputs.copyable_trait_names()
    assert 'subcortical_seg' in out_fields
    assert 'out_report' in out_fields
    # the collision map is internal, not a workflow output
    assert 'subcortical_overlap' not in out_fields


def test_subcortical_alignment_wf_debug_reportlet():
    wf = init_subcortical_mni_alignment_wf(debug=True)

    assert 'overlap_rpt' in set(wf.list_node_names())
    assert 'out_report_overlap' in wf.get_node('outputnode').inputs.copyable_trait_names()

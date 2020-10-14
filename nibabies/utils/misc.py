# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Miscellaneous utilities."""

def fix_multi_source_name(in_files, modality='T1w'):
    """
    Make up a generic source name when there are multiple
    >>> fix_multi_source_name([
    ...     '/path/to/sub-045_ses-test_T1w.nii.gz',
    ...     '/path/to/sub-045_ses-retest_T1w.nii.gz'])
    '/path/to/sub-045_T1w.nii.gz'
    """
    import os
    from nipype.utils.filemanip import filename_to_list

    if len(in_files) == 1:
        return in_files[0]

    base, in_file = os.path.split(filename_to_list(in_files)[0])
    subject_label = in_file.split("_", 1)[0].split("-")[1]
    return os.path.join(base, f"sub-{subject_label}_{modality}.nii.gz")


def check_deps(workflow):
    """Make sure dependencies are present in this system."""
    from nipype.utils.filemanip import which
    return sorted(
        (node.interface.__class__.__name__, node.interface._cmd)
        for node in workflow._get_all_nodes()
        if (hasattr(node.interface, '_cmd')
            and which(node.interface._cmd.split()[0]) is None))
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Miscellaneous utilities."""


def fix_multi_source_name(in_files):
    """
    Make up a generic source name when there are multiple
    >>> fix_multi_source_name([
    ...     '/path/to/sub-045_ses-test_T1w.nii.gz',
    ...     '/path/to/sub-045_ses-retest_T1w.nii.gz'])
    '/path/to/sub-045_T1w.nii.gz'
    """
    import re
    from pathlib import Path
    from nipype.utils.filemanip import filename_to_list

    if not isinstance(in_files, (tuple, list)):
        return in_files
    elif len(in_files) == 1:
        return in_files[0]

    p = Path(filename_to_list(in_files)[0])
    # subject_label = p.name.split("_", 1)[0].split("-")[1]
    try:
        subj = re.search(r'(?<=^sub-)[a-zA-Z0-9]*', p.name).group()
        suffix = re.search(r'(?<=_)\w+(?=\.)', p.name).group()
    except AttributeError:
        raise AttributeError("Could not extract BIDS information")
    return str(p.parent / f"sub-{subj}_{suffix}.nii.gz")


def check_deps(workflow):
    """Make sure dependencies are present in this system."""
    from nipype.utils.filemanip import which
    return sorted(
        (node.interface.__class__.__name__, node.interface._cmd)
        for node in workflow._get_all_nodes()
        if (hasattr(node.interface, '_cmd')
            and which(node.interface._cmd.split()[0]) is None))


def cohort_by_months(template, months):
    """
    Produce a recommended cohort based on partipants age
    """
    cohort_key = {
        'MNIInfant': (
            # upper bound of template | cohort
            2,  # 1
            5,  # 2
            8,  # 3
            11,  # 4
            14,  # 5
            17,  # 6
            21,  # 7
            27,  # 8
            33,  # 9
            44,  # 10
            60,  # 11
        ),
        'UNCInfant': (
            8,  # 1
            12,  # 2
            24,  # 3
        ),
    }
    ages = cohort_key.get(template)
    if ages is None:
        raise KeyError("Template cohort information does not exist.")

    for cohort, age in enumerate(ages, 1):
        if months <= age:
            return cohort
    raise KeyError("Age exceeds all cohorts!")


def check_total_memory(recommended_gb):
    """
    Check total memory allocated to the process, and compare with a recommended value.
    If available memory is equal to or greater than recommended, return ``True``.
    Otherwise, return ``False``.
    """

    try:
        import psutil
    except ImportError:
        return

    tot = int(psutil.virtual_memory().total / 1024 ** 3)
    return tot >= recommended_gb


def combine_meepi_source(in_files):
    """
    Create a new source name when optimally
    combining multiple multi-echo EPIs
    >>> combine_meepi_source([
    ...     'sub-01_run-01_echo-1_bold.nii.gz',
    ...     'sub-01_run-01_echo-2_bold.nii.gz',
    ...     'sub-01_run-01_echo-3_bold.nii.gz',])
    'sub-01_run-01_bold.nii.gz'
    """
    import os
    from nipype.utils.filemanip import filename_to_list
    base, in_file = os.path.split(filename_to_list(in_files)[0])
    entities = [ent for ent in in_file.split('_') if not ent.startswith('echo-')]
    basename = '_'.join(entities)
    return os.path.join(base, basename)

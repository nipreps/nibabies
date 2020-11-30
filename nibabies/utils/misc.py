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
        subj = re.search('(?<=^sub-)[a-zA-Z0-9]*', p.name).group()
        suffix = re.search('(?<=_)\w+(?=\.)', p.name).group()
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
    cohort_key = {
        'MNIInfant': (2, 5, 8, 11, 14, 17, 21, 27, 33, 44, 60),
        'UNCInfant': (8, 12, 24),
    }
    ages = cohort_key.get(template)
    if ages is None:
        raise KeyError("Template cohort information does not exist.")

    for cohort, age in enumerate(ages, 1):
        if months <= age:
            return cohort

    raise KeyError("Age exceeds all cohorts!")
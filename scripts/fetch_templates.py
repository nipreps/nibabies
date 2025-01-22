#!/usr/bin/env python

'Preemptive caching of commonly used TemplateFlow templates'
import templateflow.api as tf


def fetch_MNI6():
    """
    Expected templates:

    tpl-MNI152NLin6Asym/tpl-MNI152NLin6Asym_res-01_T1w.nii.gz
    tpl-MNI152NLin6Asym/tpl-MNI152NLin6Asym_res-02_T1w.nii.gz
    tpl-MNI152NLin6Asym/tpl-MNI152NLin6Asym_res-01_desc-brain_mask.nii.gz
    tpl-MNI152NLin6Asym/tpl-MNI152NLin6Asym_res-02_desc-brain_mask.nii.gz
    tpl-MNI152NLin6Asym/tpl-MNI152NLin6Asym_res-02_atlas-HCP_dseg.nii.gz
    """
    template = 'MNI152NLin6Asym'

    tf.get(template, resolution=(1, 2), desc=None, suffix='T1w')
    tf.get(template, resolution=(1, 2), desc='brain', suffix='mask')
    # CIFTI
    tf.get(template, resolution=2, atlas='HCP', suffix='dseg')


def fetch_UNCInfant():
    """
    Expected templates:

    tpl-UNCInfant/cohort-1/tpl-UNCInfant_cohort-1_T1w.nii.gz
    tpl-UNCInfant/cohort-1/tpl-UNCInfant_cohort-1_label-brain_probseg.nii.gz
    tpl-UNCInfant/cohort-1/tpl-UNCInfant_cohort-1_label-brain_mask.nii.gz
    tpl-UNCInfant/cohort-1/tpl-UNCInfant_cohort-1_label-BrainCerebellumExtraction_mask.nii.gz
    """
    template = 'UNCInfant'

    tf.get(template, cohort=1, desc=None, suffix='T1w')
    tf.get(template, cohort=1, label='brain', suffix='probseg')
    tf.get(template, cohort=1, label='brain', suffix='mask')
    tf.get(template, cohort=1, label='BrainCerebellumExtraction', suffix='mask')


def fetch_fsaverage():
    """
    Expected templates:

    tpl-fsaverage/tpl-fsaverage_hemi-L_den-164k_desc-std_sphere.surf.gii
    tpl-fsaverage/tpl-fsaverage_hemi-R_den-164k_desc-std_sphere.surf.gii
    tpl-fsaverage/tpl-fsaverage_hemi-L_den-164k_desc-vaavg_midthickness.shape.gii
    tpl-fsaverage/tpl-fsaverage_hemi-R_den-164k_desc-vaavg_midthickness.shape.gii
    tpl-fsaverage/tpl-fsaverage_hemi-L_den-41k_sphere.surf.gii
    tpl-fsaverage/tpl-fsaverage_hemi-R_den-41k_sphere.surf.gii
    """
    template = 'fsaverage'

    tf.get(template, density='164k', desc='std', suffix='sphere')
    tf.get(template, density='164k', desc='vaavg', suffix='midthickness')
    tf.get(template, density='41k', desc=None, suffix='sphere', extension='.surf.gii')


def fetch_fsLR():
    """
    Expected templates:

    tpl-fsLR/tpl-fsLR_hemi-L_den-32k_desc-nomedialwall_dparc.label.gii
    tpl-fsLR/tpl-fsLR_hemi-L_den-32k_desc-vaavg_midthickness.shape.gii
    tpl-fsLR/tpl-fsLR_hemi-L_den-32k_sphere.surf.gii
    tpl-fsLR/tpl-fsLR_hemi-R_den-32k_desc-nomedialwall_dparc.label.gii
    tpl-fsLR/tpl-fsLR_hemi-R_den-32k_desc-vaavg_midthickness.shape.gii
    tpl-fsLR/tpl-fsLR_hemi-R_den-32k_sphere.surf.gii
    tpl-fsLR/tpl-fsLR_space-fsaverage_hemi-L_den-32k_sphere.surf.gii
    tpl-fsLR/tpl-fsLR_space-fsaverage_hemi-R_den-32k_sphere.surf.gii
    """
    tf.get('fsLR', density='32k')
    tf.get('fsLR', density='164k', desc='std', suffix='sphere')
    tf.get('fsLR', density='164k', suffix='midthickness')


def fetch_MNIInfant(cohort=1):
    """
    Expected templates:

    tpl-MNIInfant/cohort-1/tpl-MNIInfant_cohort-1_res-1_T1w.nii.gz
    tpl-MNIInfant/cohort-1/tpl-MNIInfant_cohort-1_res-1_T2w.nii.gz
    tpl-MNIInfant/cohort-1/tpl-MNIInfant_cohort-1_res-1_desc-brain_mask.nii.gz
    tpl-MNIInfant/cohort-1/tpl-MNIInfant_cohort-1_res-2_T1w.nii.gz
    tpl-MNIInfant/cohort-1/tpl-MNIInfant_cohort-1_res-2_T1w.nii.gz
    tpl-MNIInfant/cohort-1/tpl-MNIInfant_cohort-1_res-2_desc-brain_mask.nii.gz
    """
    template = 'MNIInfant'

    tf.get(template, cohort=cohort, suffix='T1w')
    tf.get(template, cohort=cohort, suffix='T2w')
    tf.get(template, cohort=cohort, desc='brain', suffix='mask')


def fetch_dhcpAsym(cohort=42):
    """
    Expected templates:

    tpl-dhcpAsym_cohort-42_hemi-L_den-32k_sphere.surf.gii
    tpl-dhcpAsym_cohort-42_hemi-R_den-32k_sphere.surf.gii
    tpl-dhcpAsym_cohort-42_space-fsaverage_hemi-L_den-41k_desc-reg_sphere.surf.gii
    tpl-dhcpAsym_cohort-42_space-fsaverage_hemi-R_den-41k_desc-reg_sphere.surf.gii
    """
    template = 'dhcpAsym'

    tf.get(template, cohort=cohort, density='32k', desc=None, suffix='sphere')
    tf.get(template, cohort=cohort, space='fsaverage', density='41k', desc='reg', suffix='sphere')


def fetch_MNI2009():
    template = 'MNI152NLin2009cAsym'

    tf.get(template, resolution=(1, 2), desc=None, suffix='T1w')
    tf.get(template, resolution=(1, 2), desc='brain', suffix='mask')
    tf.get(template, resolution=2, desc='fMRIPrep', suffix='boldref')
    tf.get(template, resolution=1, label='brain', suffix='probseg')


def main():
    fetch_MNI6()
    fetch_UNCInfant()
    fetch_fsaverage()
    fetch_fsLR()
    fetch_MNIInfant()
    fetch_dhcpAsym()
    fetch_MNI2009()


if __name__ == '__main__':
    main()

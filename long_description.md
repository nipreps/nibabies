Magnetic resonance imaging (MRI) requires a set of preprocessing steps before
any statistical analysis. In an effort to standardize preprocessing,
we developed [fMRIPrep](https://fmriprep.org/en/stable/) (a preprocessing tool
for functional MRI, fMRI), and generalized its standardization approach to
other neuroimaging modalities ([NiPreps](https://www.nipreps.org/)). NiPreps
brings standardization and ease of use to the researcher, and effectively
limits the methodological variability within preprocessing. fMRIPrep is designed
to be used across wide ranges of populations; however it is designed for (and
evaluated with) human adult datasets. Infant MRI (i.e., 0-2 years) presents
unique challenges due to head size (e.g., reduced SNR and increased partial
voluming and rapid shifting in tissue contrast due to myelination. These and
other challenges require a more specialized workflow. *NiBabies*, an open-source
pipeline extending from fMRIPrep for infant structural and functional MRI
preprocessing, aims to address this need.

The workflow is built atop [Nipype](https://nipype.readthedocs.io) and encompasses a large
set of tools from well-known neuroimaging packages, including
[FSL](https://fsl.fmrib.ox.ac.uk/fsl/fslwiki/),
[ANTs](https://stnava.github.io/ANTs/),
[FreeSurfer](https://surfer.nmr.mgh.harvard.edu/),
[AFNI](https://afni.nimh.nih.gov/),
[Connectome Workbench](https://humanconnectome.org/software/connectome-workbench),
and [Nilearn](https://nilearn.github.io/).
This pipeline was designed to provide the best software implementation for each state of
preprocessing, and will be updated as newer and better neuroimaging software becomes
available.

*NiBabies* performs basic preprocessing steps (coregistration, normalization, unwarping,
segmentation, skullstripping etc.) providing outputs that can be
easily submitted to a variety of group level analyses, including task-based or resting-state
fMRI, graph theory measures, surface or volume-based statistics, etc.
*NiBabies* allows you to easily do the following:

  * Take fMRI data from *unprocessed* (only reconstructed) to ready for analysis.
  * Implement tools from different software packages.
  * Achieve optimal data processing quality by using the best tools available.
  * Generate preprocessing-assessment reports, with which the user can easily identify problems.
  * Receive verbose output concerning the stage of preprocessing for each subject, including
    meaningful errors.
  * Automate and parallelize processing steps, which provides a significant speed-up from
    typical linear, manual processing.

[Repository](https://github.com/nipreps/nibabies)
[Documentation](https://nibabies.readthedocs.io/en/stable/)

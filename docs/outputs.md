# Outputs

NiBabies outputs conform to the BIDS Derivatives specification (see BIDS Derivatives, along with the upcoming BEP 011 and BEP 012). NiBabies generates three broad classes of outcomes:


# Processing level
As of version 24.0.0, NiBabies supports three levels of derivatives:

--level minimal: This processing mode aims to produce the smallest working directory and output dataset possible, while enabling all further processing results to be deterministically generated. Most components of the visual reports can be generated at this level, so the quality of preprocessing can be assessed. Because no resampling is done, confounds and carpetplots will be missing from the reports.

--level resampling: This processing mode aims to produce additional derivatives that enable third-party resampling, resampling BOLD series in the working directory as needed, but these are not saved to the output directory. The --me-output-echos flag will be enabled at this level, in which case the individual echos will be saved to the working directory after slice-timing correction, head-motion correction, and susceptibility distortion correction.

--level full: This processing mode aims to produce all derivatives that have previously been a part of the NiBabies output dataset. This is the default processing level.

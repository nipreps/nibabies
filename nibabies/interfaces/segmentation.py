from nipype.interfaces.base import (
    traits, TraitedSpec, File, Directory,
    CommandLine, CommandLineInputSpec
)

class MultiSegPipelineInputSpec(CommandLineInputSpec):
    t1_file = File(
        mandatory=True,
        argstr='--T1 %s',
        exists=True,
        help="Path of the T1 image to segment",
    )
    t2_file = File(
        mandatory=True,
        argstr='--T2 %s',
        exists=True,
        help="Path of the T2 image to segment",
    )
    mask_file = File(
        mandatory=True,
        argstr='--mask %s',
        exists=True,
        help="Path of the brain mask",
    )
    dwi_file = File(
        argstr='--DWI %s',
        exists=True,
        help="Path of the DWI image",
    )
    data_file = File(
        argstr='--data %s',
        exists=True,
        help="Path of the XML file that contains the data configuration",
    )
    parameters_file = File(
        argstr='--parameters %s',
        exists=True,
        help="Path of the XML file that contains the parameters configuration",
    )
    executables_file = File(
        argstr='--executables %s',
        exists=True,
        help="Path of the XML file that contains the executables configuration",
    )
    return_parameter_file = File(
        argstr='--returnparameterfile %s',
        help='Filename in which to write simple return parameters (int, float, '
             'int-vector, etc.) as opposed to bulk return parameters (image, '
             'geometry, transform, measurement, table).',
    )
    force = traits.Bool(
        argstr='--force',
        help='Runs the tool even if errors are detected while loading parameter files',
    )
    suffix = traits.Str(
        argstr='--suffix %s',
        help="Suffix added to output files (default: NP)",
    )
    prefix = traits.Str(
        argstr='--prefix %s',
        help="Prefix added to output files (default: neo)",
    )
    output_directory = Directory(
        argstr='--output %s',
        help="Path of the output directory where all files will be written",
    )


class MultiSegPipelineOutputSpec(TraitedSpec):
    output_dir = Directory(exists=True, help="Output directory")


class MultiSegPipeline(CommandLine):
    """
    A tool to create an sharp atlas specific atlas for Neoseg.
    For any suggestion or question, please contact cherel[at]unc.edu

    Author(s): Kevin Pham, Marie Cherel, Francois Budin, Juan Carlos Prieto,
    Martin Styner
    """
    _cmd = "MultiSegPipeline2.1.0 --noGui"
    input_spec = MultiSegPipelineInputSpec
    output_spec = MultiSegPipelineOutputSpec


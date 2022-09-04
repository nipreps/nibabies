# -*- coding: utf-8 -*-
# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""This module provides interfaces for workbench surface commands"""
import os

from nipype.interfaces.base import TraitedSpec, File, traits, CommandLineInputSpec
from nipype.interfaces.workbench.base import WBCommand
from nipype import logging

iflogger = logging.getLogger("nipype.interface")


class VolumeToSurfaceMappingInputSpec(CommandLineInputSpec):

    in_file = File(
        exists=True,
        mandatory=True,
        argstr="%s ",
        position=0,
        desc="the volume to map data from",
    )

    surface = File(
        exists=True,
        mandatory=True,
        argstr="%s ",
        position=1,
        desc="the surface to map the data onto",
    )

    out_file = File(
        name_source=["in_file"],
        name_template="%s.func.gii",
        keep_extension=False,
        argstr="%s ",
        position=2,
        desc="output - the output metric file",
    )
    
    mapping_method = traits.Enum(
        "trilinear",
        "enclosing",
        "cubic",
        "ribbon-constrained",
        "myelin-style",
        position=3,
        argstr="-%s ",
        desc="choose mapping method: trilinear, enclosing voxel, cubic spline, "
             "ribbon-constrained, or myelin-style",
    )

    inner_surf = File(
        argstr="%s ",
        position=4,
        desc="the inner surface of the ribbon",
    )

    outer_surf = File(
        argstr="%s ",
        position=5,
        desc="the outer surface of the ribbon",
    )

    roi_volume = File(
        argstr="-volume-roi %s ",
        position=6,
        desc="the roi volume file to use",
    )

    weighted = traits.Bool(
        position=7,
        argstr="-weighted ",
        desc="treat the roi values as weightings rather than binary",
    )

    subdiv_num = traits.Int(
        position=8,
        argstr="-voxel-subdiv %d ",
        desc="number of subdivisions while estimating voxel weights, "
        "default 3 if this option is not specified",
    )

    thin_columns = traits.Bool(
        position=9,
        argstr="-thin-columns ",
        desc="use non-overlapping polyhedra",
    )

    gaussian_scale = traits.Float(
        argstr="-gaussian %f ",
        position=10,
        desc="value to multiply the local thickness by, to get the gaussian sigma. "
             "reduce weight to voxels that aren't near <surface>",
    )

    interpolate_method = traits.Enum(
        "CUBIC",
        "ENCLOSING_VOXEL",
        "TRILINEAR",
        position=11,
        argstr="-interpolate %s ",
        desc="must be CUBIC, ENCLOSING_VOXEL, or TRILINEAR. "
             "instead of a weighted average of voxels, "
             "interpolate at subpoints inside the ribbon",
    )
    
    bad_vertices_out = File(
        argstr="-bad-vertices-out %s ",
        position=12,
        desc="output an ROI of which vertices didn't intersect any valid voxels",
    )

    output_weights = traits.Bool(
        position=13,
        argstr="-output-weights ",
        desc="the column number",
    )

    output_weights_vertex = traits.Int(
        position=14,
        argstr="%d ",
        desc="the column number",
    )

    output_weights_weights_out = File(
        argstr="%s ",
        position=15,
        desc="output - volume to write the weights to",
    )

    output_weights_text_out = File(
        argstr="%s ",
        position=16,
        desc="text file name. "
        "write the voxel weights for all vertices to a text file",
    )

    myelin_style_ribbon_roi = File(
        argstr="%s ",
        position=17,
        desc="roi volume of the cortical ribbon for this hemisphere",
    )

    myelin_style_thickness = File(
        argstr="%s ",
        position=18,
        desc="a metric file of cortical thickness",
    )

    myelin_style_sigma = traits.Float(
        argstr="%f ",
        position=19,
        desc="gaussian kernel in mm for weighting voxels within range",
    )

    legacy_bug = traits.Bool(
        position=20,
        argstr="-legacy-bug ",
        desc="emulate old Workbench v1.2.3 and earlier code that didn't follow a cylinder cutoff",
    )

    subvol_select = File(
        argstr="-subvol-select %s ",
        position=21,
        desc="the subvolume number or name",
    )


class VolumeToSurfaceMappingOutputSpec(TraitedSpec):
    out_file = File(exists=True, desc="the output metric file")
    bad_vertices_out = File(desc="ROI of which vertices didn't intersect any valid voxels")
    output_weights_weights_out = File(desc="volume to which the specified vertex's weights are written")
    output_weights_text_out = File(desc="text file of the voxel weights for all vertices")


class VolumeToSurfaceMapping(WBCommand):
    """
    You must specify exactly one mapping method.  Enclosing voxel uses the
    value from the voxel the vertex lies inside, while trilinear does a 3D
    linear interpolation based on the voxels immediately on each side of the
    vertex's position.

    The ribbon mapping method constructs a polyhedron from the vertex's
    neighbors on each surface, and estimates the amount of this polyhedron's
    volume that falls inside any nearby voxels, to use as the weights for
    sampling.  If -thin-columns is specified, the polyhedron uses the edge
    midpoints and triangle centroids, so that neighboring vertices do not
    have overlapping polyhedra.  This may require increasing -voxel-subdiv to
    get enough samples in each voxel to reliably land inside these smaller
    polyhedra.  The volume ROI is useful to exclude partial volume effects of
    voxels the surfaces pass through, and will cause the mapping to ignore
    voxels that don't have a positive value in the mask.  The subdivision
    number specifies how it approximates the amount of the volume the
    polyhedron intersects, by splitting each voxel into NxNxN pieces, and
    checking whether the center of each piece is inside the polyhedron.  If
    you have very large voxels, consider increasing this if you get zeros in
    your output.  The -gaussian option makes it act more like the myelin
    method, where the distance of a voxel from <surface> is used to
    downweight the voxel.  The -interpolate suboption, instead of doing a
    weighted average of voxels, interpolates from the volume at the
    subdivided points inside the ribbon.  If using both -interpolate and the
    -weighted suboption to -volume-roi, the roi volume weights are linearly
    interpolated, unless the -interpolate method is ENCLOSING_VOXEL, in which
    case ENCLOSING_VOXEL is also used for sampling the roi volume weights.

    The myelin style method uses part of the caret5 myelin mapping command to
    do the mapping: for each surface vertex, take all voxels that are in a
    cylinder with radius and height equal to cortical thickness, centered on
    the vertex and aligned with the surface normal, and that are also within
    the ribbon ROI, and apply a gaussian kernel with the specified sigma to
    them to get the weights to use.  The -legacy-bug flag reverts to the
    unintended behavior present from the initial implementation up to and
    including v1.2.3, which had only the tangential cutoff and a bounding box
    intended to be larger than where the cylinder cutoff should have been.

    >>> from nipype.interfaces.workbench import VolumeToSurfaceMapping
    >>> vtsm = VolumeToSurfaceMapping()
    >>> vtsm.inputs.in_file = 'sub-01_task-rest_run-1_bold.nii.gz'
    >>> vtsm.inputs.surface = 'sub-01.L.midthickness.surf.gii'
    >>> vtsm.inputs.ribbon_constrained = True
    >>> vtsm.inputs.inner_surf = 'sub-01.L.white.surf.gii'
    >>> vtsm.inputs.outer_surf = 'sub-01.L.pial.surf.gii'
    >>> vtsm.inputs.roi_volume = 'sub-01_task-rest_run-1_goodvoxels.nii.gz'
    >>> vtsm.cmdline
    'wb_command -volume-to-surface-mapping \
    sub-01_task-rest_run-1_bold.nii.gz \
    sub-01.L.midthickness.surf.gii \
    sub-01_task-rest_run-1_bold.func.gii \
    -ribbon-constrained sub-01.L.white.surf.gii sub-01.L.pial.surf.gii \
    -volume-roi sub-01_task-rest_run-1_goodvoxels.nii.gz
    """
    input_spec = VolumeToSurfaceMappingInputSpec
    output_spec = VolumeToSurfaceMappingOutputSpec
    _cmd = "wb_command -volume-to-surface-mapping "

    def _format_arg(self, opt, spec, val):
        if opt in [
            "inner_surf",
            "outer_surf",
            "roi_volume",
            "subdiv_num",
            "thin_columns",
            "gaussian_scale",
            "interpolate_method",
            "bad_vertices_out",
            "output_weights",
            "output_weights_text_out",
        ]:               
            if self.inputs.mapping_method.val is not "ribbon-constrained":
                raise ValueError(
                    "{} was set, but ribbon-constrained mapping method was not".format(opt)
                )

        if opt == "weighted":
            if (val == True) and (self.inputs.mapping_method.val is not "ribbon-constrained"):
                raise ValueError(
                    "weighted was set, but the required roi_volume was not"
                )

        if ((self.inputs.mapping_method.val == "ribbon-constrained")
            and not (self.inputs.inner_surf and self.inputs.outer_surf)
            ):
                raise ValueError(
                    "mapping method is ribbon-constrained but at least one of "
                    "inner_surf, outer_surf was not set"
                )

        if opt == "output_weights":
            if ((val == True) and not (self.inputs.output_weights_vertex
                                       and self.inputs.output_weights_weights_out)
            ):
                raise ValueError(
                    "output_weights was set but at least one of "
                    "output_weights_vertex, output_weights_weights_out was not set"
                )

        if ((self.inputs.mapping_method.val == "myelin-style")
             and not (self.inputs.myelin_style_ribbon_roi
                      and self.inputs.myelin_style_thickness
                      and self.inputs.myelin_style_sigma)
            ):
                raise ValueError(
                    "mapping method is myelin-style but at least one of myelin_style_ribbon_roi, "
                    "myelin_style_thickness, myelin_style_sigma was not set"
                )

        if opt == "legacy_bug":             
            if (val == True) and (not self.inputs.mapping_method.val is not "myelin-style"):
                raise ValueError(
                    "legacy_bug was set, but the mapping method is not myelin-style"
                )
        
        return super(VolumeToSurfaceMapping, self)._format_arg(opt, spec, val)

    def _list_outputs(self):
        outputs = super(VolumeToSurfaceMapping, self)._list_outputs()
        return outputs
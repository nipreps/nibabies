# emacs: -*- mode: python; py-indent-offset: 4; indent-tabs-mode: nil -*-
# vi: set ft=python sts=4 ts=4 sw=4 et:
"""Interfaces to generate reportlets."""

import logging
import os
import re
import time
from collections import Counter
from pathlib import Path

from nipype.interfaces.base import (
    BaseInterfaceInputSpec,
    Directory,
    File,
    InputMultiObject,
    SimpleInterface,
    Str,
    TraitedSpec,
    isdefined,
    traits,
)
from niworkflows.interfaces.reportlets import base as nrb

LOGGER = logging.getLogger('nipype.interface')

SUBJECT_TEMPLATE = """\
\t<ul class="elem-desc">
\t\t<li>Subject ID: {subject_id}</li>
\t\t<li>Session ID: {session_id}</li>
\t\t<li>Chronological age (months): {age}</li>
\t\t<li>Structural images: {num_t1w} T1-weighted, {num_t2w} T2-weighted</li>
\t\t<li>Anatomical reference space: {anat_ref}<li>
\t\t<li>Functional series: {n_bold:d}</li>
{tasks}
\t\t<li>Standard output spaces: {std_spaces}</li>
\t\t<li>Non-standard output spaces: {nstd_spaces}</li>
\t\t<li>Surface reconstruction method: {recon_method}</li>
\t\t<li>Surface reconstruction status: {recon_status}</li>
\t</ul>
"""

FUNCTIONAL_TEMPLATE = """\
\t\t<details open>
\t\t<summary>Summary</summary>
\t\t<ul class="elem-desc">
\t\t\t<li>Original orientation: {ornt}</li>
\t\t\t<li>Repetition time (TR): {tr:.03g}s</li>
\t\t\t<li>Phase-encoding (PE) direction: {pedir}</li>
\t\t\t<li>{multiecho}</li>
\t\t\t<li>Slice timing correction: {stc}</li>
\t\t\t<li>Susceptibility distortion correction: {sdc}</li>
\t\t\t<li>Registration: {registration}</li>
\t\t\t<li>Non-steady-state volumes: {dummy_scan_desc}</li>
\t\t</ul>
\t\t</details>
"""

ABOUT_TEMPLATE = """\t<ul>
\t\t<li>NiBabies version: {version}</li>
\t\t<li>NiBabies command: <code>{command}</code></li>
\t\t<li>Date preprocessed: {date}</li>
\t</ul>
</div>
"""


# TODO: Move to niworkflows
class SummaryOutputSpec(TraitedSpec):
    out_report = File(exists=True, desc='HTML segment containing summary')


class SummaryInterface(SimpleInterface):
    output_spec = SummaryOutputSpec

    def _run_interface(self, runtime):
        segment = self._generate_segment()
        fname = os.path.join(runtime.cwd, 'report.html')
        with open(fname, 'w') as fobj:
            fobj.write(segment)
        self._results['out_report'] = fname
        return runtime

    def _generate_segment(self):
        raise NotImplementedError


class SubjectSummaryInputSpec(BaseInterfaceInputSpec):
    t1w = InputMultiObject(File(exists=True), desc='T1w structural images')
    t2w = InputMultiObject(File(exists=True), desc='T2w structural images')
    subjects_dir = Directory(desc='Infant FreeSurfer subjects directory')
    subject_id = Str(mandatory=True, desc='Subject ID')
    session_id = Str(desc='Session ID')
    anatomical_reference = traits.Enum('T1w', 'T2w', mandatory=True)
    bold = InputMultiObject(
        traits.Either(File(exists=True), traits.List(File(exists=True))),
        desc='BOLD functional series',
    )
    std_spaces = traits.List(Str, desc='list of standard spaces')
    nstd_spaces = traits.List(Str, desc='list of non-standard spaces')
    recon_method = traits.Enum(
        None,
        'freesurfer',
        'infantfs',
        'mcribs',
        desc='surface reconstruction method',
    )
    age = traits.Int(desc='Chronological age in months at the time of session')


class SubjectSummaryOutputSpec(SummaryOutputSpec):
    # This exists to ensure that the summary is run prior to the first ReconAll
    # call, allowing a determination whether there is a pre-existing directory
    subject_id = Str(desc='Surface reconstruction subject ID')


class SubjectSummary(SummaryInterface):
    input_spec = SubjectSummaryInputSpec
    output_spec = SubjectSummaryOutputSpec

    def _run_interface(self, runtime):
        if self.inputs.subject_id:
            self._recon_id = f'sub-{self.inputs.subject_id}'
            if self.inputs.session_id:
                self._recon_id += f'_ses-{self.inputs.session_id}'
            self._results['subject_id'] = self._recon_id
        return super()._run_interface(runtime)

    def _generate_segment(self):
        BIDS_NAME = re.compile(
            r'^(.*\/)?'
            '(?P<subject_id>sub-[a-zA-Z0-9]+)'
            '(_(?P<session_id>ses-[a-zA-Z0-9]+))?'
            '(_(?P<task_id>task-[a-zA-Z0-9]+))?'
            '(_(?P<acq_id>acq-[a-zA-Z0-9]+))?'
            '(_(?P<rec_id>rec-[a-zA-Z0-9]+))?'
            '(_(?P<run_id>run-[a-zA-Z0-9]+))?'
        )

        recon_method = self.inputs.recon_method

        statuses = {'no': 'Not run', 'todo': 'Run by NiBabies', 'done': 'Pre-existing directory'}
        if not self.inputs.subjects_dir:
            recon_status = statuses['no']
        else:
            if recon_method == 'freesurfer':
                from smriprep.interfaces.freesurfer import ReconAll

                recon = ReconAll(
                    subjects_dir=self.inputs.subjects_dir,
                    subject_id=self._recon_id,
                    T1_files=self.inputs.t1w,
                    flags='-noskullstrip',
                )
                recon_status = (
                    statuses['done'] if recon.cmdline.startswith('echo') else statuses['todo']
                )

            elif recon_method == 'infantfs':
                from niworkflows.utils.connections import pop_file

                from nibabies.interfaces.freesurfer import InfantReconAll

                recon = InfantReconAll(
                    subjects_dir=self.inputs.subjects_dir,
                    subject_id=self._recon_id,
                    t1_file=pop_file(self.inputs.t1w),
                )
                recon_status = (
                    statuses['done'] if recon.cmdline.startswith('echo') else statuses['todo']
                )

            elif recon_method == 'mcribs':
                # Use fingerprint logfile to identify "MCRIBS" runs vs FreeSurfer
                fingerprint = (
                    Path(self.inputs.subjects_dir) / self._recon_id / 'scripts' / 'mcribs.log'
                )
                recon_status = statuses['done'] if fingerprint.exists() else statuses['todo']

        num_t1w, num_t2w = 0, 0
        if self.inputs.t1w:
            num_t1w = len(self.inputs.t1w)
        if self.inputs.t2w:
            num_t2w = len(self.inputs.t2w)

        # Add list of tasks with number of runs
        bold_series = self.inputs.bold if isdefined(self.inputs.bold) else []
        bold_series = [s[0] if isinstance(s, list) else s for s in bold_series]

        counts = Counter(
            BIDS_NAME.search(series).groupdict()['task_id'][5:] for series in bold_series
        )

        tasks = ''
        if counts:
            header = '\t\t<ul class="elem-desc">'
            footer = '\t\t</ul>'
            lines = [
                '\t\t\t<li>Task: {task_id} ({n_runs:d} run{s})</li>'.format(
                    task_id=task_id, n_runs=n_runs, s='' if n_runs == 1 else 's'
                )
                for task_id, n_runs in sorted(counts.items())
            ]
            tasks = '\n'.join([header] + lines + [footer])

        return SUBJECT_TEMPLATE.format(
            subject_id=self.inputs.subject_id,
            session_id=self.inputs.session_id,
            age=self.inputs.age,
            num_t1w=num_t1w,
            num_t2w=num_t2w,
            anat_ref=self.inputs.anatomical_reference,
            n_bold=len(bold_series),
            tasks=tasks,
            std_spaces=', '.join(self.inputs.std_spaces),
            nstd_spaces=', '.join(self.inputs.nstd_spaces),
            recon_method=recon_method,
            recon_status=recon_status,
        )


class FunctionalSummaryInputSpec(TraitedSpec):
    slice_timing = traits.Enum(False, True, 'TooShort', desc='Slice timing correction used')
    distortion_correction = traits.Str(
        desc='Susceptibility distortion correction method', mandatory=True
    )
    pe_direction = traits.Enum(
        None,
        'i',
        'i-',
        'j',
        'j-',
        'k',
        'k-',
        mandatory=True,
        desc='Phase-encoding direction detected',
    )
    registration = traits.Enum(
        'FSL', 'FreeSurfer', mandatory=True, desc='Functional/anatomical registration method'
    )
    fallback = traits.Bool(desc='Boundary-based registration rejected')
    registration_dof = traits.Enum(
        6, 9, 12, desc='Registration degrees of freedom', mandatory=True
    )
    registration_init = traits.Enum(
        't1w',
        't2w',
        'header',
        mandatory=True,
        desc='Whether to initialize registration with the "header"'
        ' or by centering the volumes ("t1w" or "t2w")',
    )
    tr = traits.Float(desc='Repetition time', mandatory=True)
    dummy_scans = traits.Either(traits.Int(), None, desc='number of dummy scans specified by user')
    algo_dummy_scans = traits.Int(desc='number of dummy scans determined by algorithm')
    echo_idx = InputMultiObject(traits.Str, usedefault=True, desc='BIDS echo identifiers')
    orientation = traits.Str(mandatory=True, desc='Orientation of the voxel axes')


class FunctionalSummary(SummaryInterface):
    input_spec = FunctionalSummaryInputSpec

    def _generate_segment(self):
        dof = self.inputs.registration_dof
        if isdefined(self.inputs.slice_timing):
            stc = {
                True: 'Applied',
                False: 'Not applied',
                'TooShort': 'Skipped (too few volumes)',
            }[self.inputs.slice_timing]
        else:
            stc = 'n/a'
        # TODO: Add a note about registration_init below?
        reg = {
            'FSL': [
                'FSL <code>flirt</code> with boundary-based registration'
                f' (BBR) metric - {dof} dof',
                'FSL <code>flirt</code> rigid registration - 6 dof',
            ],
            'FreeSurfer': [
                'FreeSurfer <code>bbregister</code> '
                f'(boundary-based registration, BBR) - {dof} dof',
                f'FreeSurfer <code>mri_coreg</code> - {dof} dof',
            ],
        }[self.inputs.registration][self.inputs.fallback]

        pedir = get_world_pedir(self.inputs.orientation, self.inputs.pe_direction)

        dummy_scan_tmp = '{n_dum}'
        if self.inputs.dummy_scans == self.inputs.algo_dummy_scans:
            dummy_scan_msg = ' '.join(
                [dummy_scan_tmp, '(Confirmed: {n_alg} automatically detected)']
            ).format(n_dum=self.inputs.dummy_scans, n_alg=self.inputs.algo_dummy_scans)
        # the number of dummy scans was specified by the user and
        # it is not equal to the number detected by the algorithm
        elif self.inputs.dummy_scans is not None:
            dummy_scan_msg = ' '.join(
                [dummy_scan_tmp, '(Warning: {n_alg} automatically detected)']
            ).format(n_dum=self.inputs.dummy_scans, n_alg=self.inputs.algo_dummy_scans)
        # the number of dummy scans was not specified by the user
        else:
            dummy_scan_msg = dummy_scan_tmp.format(n_dum=self.inputs.algo_dummy_scans)

        multiecho = 'Single-echo EPI sequence.'
        n_echos = len(self.inputs.echo_idx)
        if n_echos == 1:
            multiecho = (
                f'Multi-echo EPI sequence: only echo {self.inputs.echo_idx[0]} processed '
                'in single-echo mode.'
            )
        if n_echos > 2:
            multiecho = f'Multi-echo EPI sequence: {n_echos} echoes.'

        return FUNCTIONAL_TEMPLATE.format(
            pedir=pedir,
            stc=stc,
            sdc=self.inputs.distortion_correction,
            registration=reg,
            tr=self.inputs.tr,
            dummy_scan_desc=dummy_scan_msg,
            multiecho=multiecho,
            ornt=self.inputs.orientation,
        )


class AboutSummaryInputSpec(BaseInterfaceInputSpec):
    version = Str(desc='NiBabies version')
    command = Str(desc='NiBabies command')
    # Date not included - update timestamp only if version or command changes


class AboutSummary(SummaryInterface):
    input_spec = AboutSummaryInputSpec

    def _generate_segment(self):
        return ABOUT_TEMPLATE.format(
            version=self.inputs.version,
            command=self.inputs.command,
            date=time.strftime('%Y-%m-%d %H:%M:%S %z'),
        )


class LabeledHistogramInputSpec(nrb._SVGReportCapableInputSpec):
    in_file = traits.File(exists=True, mandatory=True, desc='Image containing values to plot')
    label_file = traits.File(
        exists=True,
        desc='Mask or label image where non-zero values will be used to extract data from in_file',
    )
    mapping = traits.Dict(desc='Map integer label values onto names of voxels')
    xlabel = traits.Str('voxels', usedefault=True, desc='Description of values plotted')


class LabeledHistogram(nrb.ReportingInterface):
    input_spec = LabeledHistogramInputSpec

    def _generate_report(self):
        import nibabel as nb
        import numpy as np
        import seaborn as sns
        from matplotlib import pyplot as plt
        from nilearn.image import resample_to_img

        report_file = self._out_report
        img = nb.load(self.inputs.in_file)
        data = img.get_fdata(dtype=np.float32)

        if self.inputs.label_file:
            label_img = nb.load(self.inputs.label_file)
            if label_img.shape != img.shape[:3] or not np.allclose(label_img.affine, img.affine):
                label_img = resample_to_img(label_img, img, interpolation='nearest')
            labels = np.uint16(label_img.dataobj)
        else:
            labels = np.uint8(data > 0)

        uniq_labels = np.unique(labels[labels > 0])
        label_map = self.inputs.mapping or {label: label for label in uniq_labels}

        rois = {label_map.get(label, label): data[labels == label] for label in label_map}
        with sns.axes_style('whitegrid'):
            fig = sns.histplot(rois, bins=50)
            fig.set_xlabel(self.inputs.xlabel)
        plt.savefig(report_file)
        plt.close()


def get_world_pedir(ornt, pe_direction):
    """Return world direction of phase encoding"""
    # TODO: Move to niworkflows
    axes = (('Right', 'Left'), ('Anterior', 'Posterior'), ('Superior', 'Inferior'))
    ax_idcs = {'i': 0, 'j': 1, 'k': 2}

    if pe_direction is not None:
        axcode = ornt[ax_idcs[pe_direction[0]]]
        inv = pe_direction[1:] == '-'

        for ax in axes:
            for flip in (ax, ax[::-1]):
                if flip[not inv].startswith(axcode):
                    return '-'.join(flip)
    LOGGER.warning(
        'Cannot determine world direction of phase encoding. '
        f'Orientation: {ornt}; PE dir: {pe_direction}'
    )
    return 'Could not be determined - assuming Anterior-Posterior'


class _CiftiSurfacesPlotInputSpec(BaseInterfaceInputSpec):
    in_file = File(exists=True, mandatory=True, desc='CIFTI dense timeseries (.dtseries.nii)')
    surface_type = traits.Enum(
        'inflated',
        'midthickness',
        'veryinflated',
        usedefault=True,
        desc='inflation level of the fsLR mesh used for rendering',
    )


class _CiftiSurfacesPlotOutputSpec(TraitedSpec):
    out_report = File(exists=True, desc='the output SVG reportlet')


class CiftiSurfacesPlot(SimpleInterface):
    """
    Render the mean BOLD of a CIFTI dense timeseries on the fsLR surfaces (QC).

    TODO: port this interface to ``nireports``
    """

    input_spec = _CiftiSurfacesPlotInputSpec
    output_spec = _CiftiSurfacesPlotOutputSpec

    def _run_interface(self, runtime):
        from nireports.reportlets.surface import cifti_surfaces_plot

        out_report = str(Path(runtime.cwd) / 'cifti_surfaces.svg')
        cifti_surfaces_plot(
            self.inputs.in_file,
            surface_type=self.inputs.surface_type,
            output_file=out_report,
        )
        self._results['out_report'] = out_report
        return runtime


class _SubcorticalAlignmentReportInputSpec(BaseInterfaceInputSpec):
    template = File(exists=True, mandatory=True, desc='standard-space anatomical underlay (T1w)')
    bold = File(exists=True, mandatory=True, desc='BOLD reference in standard space')
    labels = File(exists=True, mandatory=True, desc='subcortical atlas label volume')
    cuts = traits.Int(7, usedefault=True, desc='number of cuts per axis')


class _SubcorticalAlignmentReportOutputSpec(TraitedSpec):
    out_report = File(exists=True, desc='before/after flicker SVG reportlet')


class SubcorticalAlignmentReport(SimpleInterface):
    """Flickering before/after of the subcortical atlas parcels.

    The parcels (colorized by structure) are overlaid first on the standard template and
    then on the BOLD reference resampled to the same space, so the flicker shows whether
    the normalized functional data sits under the matching parcels. TODO: upstream to
    nireports.
    """

    input_spec = _SubcorticalAlignmentReportInputSpec
    output_spec = _SubcorticalAlignmentReportOutputSpec

    def _run_interface(self, runtime):
        self._results['out_report'] = _subcortical_alignment_report(
            self.inputs.template,
            self.inputs.bold,
            self.inputs.labels,
            cuts=self.inputs.cuts,
            out_file=str(Path(runtime.cwd) / 'subcortical_alignment.svg'),
        )
        return runtime


def _subcortical_alignment_report(template, bold, labels, cuts=7, out_file=None):
    from uuid import uuid4

    import nibabel as nb
    import numpy as np
    from nilearn.plotting import plot_roi
    from niworkflows.viz.utils import compose_view, cuts_from_bbox, extract_svg
    from svgutils.transform import fromstring

    lab_img = nb.load(labels)
    lab = np.asanyarray(lab_img.dataobj)
    # Remap the sparse, uneven label values to consecutive integers for distinct colors
    structures = sorted({int(v) for v in np.unique(lab)} - {0})
    lut = np.zeros(int(structures[-1]) + 1 if structures else 1, dtype='int16')
    for new, old in enumerate(structures, start=1):
        lut[old] = new
    labels_c = nb.Nifti1Image(lut[lab], lab_img.affine, lab_img.header)
    vmax = max(len(structures), 1)

    cut_coords = cuts_from_bbox(lab_img, cuts=cuts)

    def _panels(bg, div_id, title):
        panels = []
        for i, mode in enumerate(('z', 'x', 'y')):
            display = plot_roi(
                labels_c,
                bg_img=bg,
                display_mode=mode,
                cut_coords=cut_coords[mode],
                cmap='tab20',
                vmin=0,
                vmax=vmax,
                alpha=0.7,
                annotate=False,
                draw_cross=False,
                colorbar=False,
                black_bg=True,
            )
            if i == 0:
                display.title(title, size=10)
            svg = extract_svg(display)
            display.close()
            svg = svg.replace('figure_1', f'{div_id}-{mode}-{uuid4()}', 1)
            panels.append(fromstring(svg))
        return panels

    out_file = str(Path(out_file or 'subcortical_alignment.svg').absolute())
    compose_view(
        _panels(template, 'before', 'template'),
        _panels(bold, 'after', 'BOLD (MNI152NLin6Asym)'),
        out_file=out_file,
    )
    return out_file

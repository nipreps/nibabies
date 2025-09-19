#!/usr/bin/env python
"""
This is a Python wrapper to facilitate running NiBabies through Docker or Singularity services.
Docker or Singularity must be installed and running.
To test installation, you can use one of the following:

```bash
docker info
singularity version
```

Please report any feedback to our GitHub repository
(https://github.com/nipreps/nibabies) and do not
forget to credit all the authors of service that NiBabies
uses (https://fmriprep.readthedocs.io/en/latest/citing.html).
"""
import os
import re
import subprocess
import sys

try:
    from ._version import __version__
except ImportError:
    __version__ = '0+unknown'
__copyright__ = 'Copyright 2023, The NiPreps Developers'
__bugreports__ = 'https://github.com/nipreps/nibabies/issues'

MISSING = """
Image '{}' is missing
Would you like to download? [Y/n] """
PKG_PATH = '/opt/conda/envs/nibabies/lib/python3.12/site-packages'
TF_TEMPLATES = (
    'MNI152Lin',
    'MNI152NLin2009cAsym',
    'MNI152NLin6Asym',
    'MNI152NLin6Sym',
    'MNIInfant',
    'MNIPediatricAsym',
    'NKI',
    'OASIS30ANTs',
    'PNC',
    'UNCInfant',
    'fsLR',
    'fsaverage',
    'fsaverage5',
    'fsaverage6',
)
NONSTANDARD_REFERENCES = ('anat', 'T1w', 'run', 'func', 'sbref', 'fsnative')

# Monkey-patch Py2 subprocess
if not hasattr(subprocess, 'DEVNULL'):
    subprocess.DEVNULL = -3

if not hasattr(subprocess, 'run'):
    # Reimplement minimal functionality for usage in this file
    def _run(args, stdout=None, stderr=None):
        from collections import namedtuple

        result = namedtuple('CompletedProcess', 'stdout stderr returncode')

        devnull = None
        if subprocess.DEVNULL in (stdout, stderr):
            devnull = open(os.devnull, 'r+')
            if stdout == subprocess.DEVNULL:
                stdout = devnull
            if stderr == subprocess.DEVNULL:
                stderr = devnull

        proc = subprocess.Popen(args, stdout=stdout, stderr=stderr)
        stdout, stderr = proc.communicate()
        res = result(stdout, stderr, proc.returncode)

        if devnull is not None:
            devnull.close()

        return res

    subprocess.run = _run


# De-fang Python 2's input - we don't eval user input
try:
    input = raw_input
except NameError:
    pass


# The helper class to facilate Docker / Singularity nuiances
class ContainerManager:
    def __init__(self, service, image=None):
        """
        Inputs
        ------
        service : str
            container service to use; either docker or singularity
        """
        self.service = service
        self.image = image
        self.command = [service, 'run']
        self.mounts = []
        if service == 'docker':
            self.add_cmd('--rm')
        elif service == 'singularity':
            self.add_cmd('--cleanenv')

    def add_cmd(self, cmd):
        """Add single or multiple commands to final container run call"""
        if isinstance(cmd, str):
            self.command.append(cmd)
        elif isinstance(cmd, (list, tuple)):
            self.command += list(cmd)

    def add_mount(self, src, dst, read_only=True):
        """
        Generate bind mount string and add to ``mounts`` attribute.

        Inputs
        ------
        src : absolute local path
        dst : absolute container path
        read_only : disable writing to bound path
        """
        self.mounts.append('{0}:{1}{2}'.format(src, dst, ':ro' if read_only else ''))

    def check_install(self):
        """Verify that the service is installed and the user has permission to
        run images.

        Returns
        -------
        -1  Docker/Singularity can't be found
        0  Docker found, but user can't connect to daemon
        1  Test run OK
        """
        try:
            ret = subprocess.run(
                [self.service, 'version'],
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
            )
        except OSError as e:
            from errno import ENOENT

            if e.errno == ENOENT:
                return -1
            raise e
        if ret.stderr.startswith(b'Cannot connect to the Docker daemon.'):
            return 0
        return 1

    def check_image(self, image):
        """Check whether image is present on local system"""
        if self.service == 'docker':
            ret = subprocess.run(['docker', 'images', '-q', image], stdout=subprocess.PIPE)
            return bool(ret.stdout)
        elif self.service == 'singularity':
            # check if the image file exists
            return os.path.exists(os.path.abspath(image))
        raise NotImplementedError

    def check_memory(self, image):
        """Check total memory from within an image"""
        if self.service == 'docker':
            cmd = ['docker', 'run', '--rm', '--entrypoint=free', image, '-m']
        elif self.service == 'singularity':
            cmd = ['singularity', 'exec', image, 'free', '-m']
        ret = subprocess.run(cmd, stdout=subprocess.PIPE)
        if ret.returncode:
            return -1
        mem = [
            line.decode().split()[1]
            for line in ret.stdout.splitlines()
            if line.startswith(b'Mem:')
        ][0]
        return int(mem)

    def set_version(self):
        if self.service == 'docker':
            ret = subprocess.run(
                ['docker', 'version', '--format', '{{.Server.Version}}'], stdout=subprocess.PIPE
            )
        elif self.service == 'singularity':
            ret = subprocess.run(['singularity', 'version'], stdout=subprocess.PIPE)
        version = ret.stdout.decode('ascii').strip()
        version_env = '%s_VERSION_8395080871' % self.service.upper()
        self.add_envvar((version_env, version))

    def add_envvar(self, envtuple):
        """Set an environment variable

        Inputs
        ------
        envtuple : tuple in the form of ("ENV_VAR", "value")
        """
        if self.service == 'docker':
            env = '='.join(envtuple)
            self.add_cmd(['-e', env])
        elif self.service == 'singularity':
            # singularity will transfer over environment variables
            # with the prefix: SINGULARITYENV_
            envvar, value = envtuple
            envvar = 'SINGULARITYENV_' + envvar
            os.environ[envvar] = value

    def finalize_container_cmd(self):
        """Add bindings to final command, and finish with image"""
        if self.service == 'docker':
            mflag = '-v'
        elif self.service == 'singularity':
            mflag = '-B'

        for mount in self.mounts:
            self.add_cmd((mflag, mount))
        self.add_cmd(self.image)


def merge_help(wrapper_help, target_help):
    def _get_posargs(usage):
        """
        Extract positional arguments from usage string.

        This function can be used by both native fmriprep (`fmriprep -h`)
        and the docker wrapper (`fmriprep-docker -h`).
        """
        posargs = []
        for targ in usage.split('\n')[-3:]:
            line = targ.lstrip()
            if line.startswith('usage'):
                continue
            if line[0].isalnum() or line[0] == '{':
                posargs.append(line)
            elif line[0] == '[' and (line[1].isalnum() or line[1] == '{'):
                posargs.append(line)
        return ' '.join(posargs)

    # Matches all flags with up to one nested square bracket
    opt_re = re.compile(r'(\[--?[\w-]+(?:[^\[\]]+(?:\[[^\[\]]+\])?)?\])')
    # Matches flag name only
    flag_re = re.compile(r'\[--?([\w-]+)[ \]]')

    # Normalize to Unix-style line breaks
    w_help = wrapper_help.rstrip().replace('\r', '')
    t_help = target_help.rstrip().replace('\r', '')

    w_usage, w_details = w_help.split('\n\n', 1)
    w_groups = w_details.split('\n\n')
    t_usage, t_details = t_help.split('\n\n', 1)
    t_groups = t_details.split('\n\n')

    w_posargs = _get_posargs(w_usage)
    t_posargs = _get_posargs(t_usage)

    w_options = opt_re.findall(w_usage)
    w_flags = sum(map(flag_re.findall, w_options), [])
    t_options = opt_re.findall(t_usage)
    t_flags = sum(map(flag_re.findall, t_options), [])

    # The following code makes this assumption
    # assert w_flags[:2] == ["h", "version"]
    # assert w_posargs.replace("]", "").replace("[", "") == t_posargs

    # Make sure we're not clobbering options we don't mean to
    overlap = set(w_flags).intersection(t_flags)
    expected_overlap = {
        'bids-database-dir',
        'bids-filter-file',
        'derivatives',
        'deriv-filter-file',
        'fs-license-file',
        'fs-subjects-dir',
        'config-file',
        'segmentation-atlases-dir',
        'h',
        'use-plugin',
        'version',
        'w',
    }

    assert overlap == expected_overlap, 'Clobbering options: {}'.format(
        ', '.join(overlap - expected_overlap)
    )

    sections = []

    # Construct usage
    start = w_usage[: w_usage.index(' [')]
    indent = ' ' * len(start)
    new_options = sum(
        (
            w_options[:2],
            [opt for opt, flag in zip(t_options, t_flags, strict=False) if flag not in overlap],
            w_options[2:],
        ),
        [],
    )
    opt_line_length = 79 - len(start)
    length = 0
    opt_lines = [start]
    for opt in new_options:
        opt = ' ' + opt
        olen = len(opt)
        if length + olen <= opt_line_length:
            opt_lines[-1] += opt
            length += olen
        else:
            opt_lines.append(indent + opt)
            length = olen
    opt_lines.append(indent + ' ' + t_posargs)
    sections.append('\n'.join(opt_lines))

    # Use target description and positional args
    sections.extend(t_groups[:2])

    for line in t_groups[2].split('\n')[1:]:
        content = line.lstrip().split(',', 1)[0]
        if content[1:] not in overlap:
            w_groups[2] += '\n' + line

    sections.append(w_groups[2])

    # All remaining sections, show target then wrapper (skipping duplicates)
    sections.extend(t_groups[3:] + w_groups[6:])
    return '\n\n'.join(sections)


def is_in_directory(filepath, directory):
    return os.path.realpath(filepath).startswith(os.path.realpath(directory) + os.sep)


def get_parser():
    """Defines the command line interface of the wrapper"""
    import argparse
    from functools import partial

    class ToDict(argparse.Action):
        def __call__(self, parser, namespace, values, option_string=None):
            d = {}
            for kv in values:
                k, v = kv.split('=')
                d[k] = os.path.abspath(v)
            setattr(namespace, self.dest, d)

    def _is_file(path, parser):
        """Ensure a given path exists and it is a file."""
        path = os.path.abspath(path)
        if not os.path.isfile(path):
            raise parser.error('Path should point to a file (or symlink of file): <%s>.' % path)
        return path

    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        add_help=False,
    )

    IsFile = partial(_is_file, parser=parser)

    # require users to specify container service
    parser.add_argument('service', nargs='?', choices=('docker', 'singularity'))
    # Standard NiBabies arguments
    parser.add_argument('bids_dir', nargs='?', type=os.path.abspath, default='')
    parser.add_argument('output_dir', nargs='?', type=os.path.abspath, default='')
    parser.add_argument(
        'analysis_level', nargs='?', choices=['participant'], default='participant'
    )

    parser.add_argument(
        '-h', '--help', action='store_true', help='show this help message and exit'
    )
    parser.add_argument(
        '--version', action='store_true', help="show program's version number and exit"
    )

    # Allow alternative images (semi-developer)
    parser.add_argument(
        '-i',
        '--image',
        metavar='IMG',
        type=str,
        default='nipreps/nibabies:%s' % __version__,
        help='image name',
    )

    # Options for mapping files and directories into container
    # Update `expected_overlap` variable in merge_help() when adding to this
    g_wrap = parser.add_argument_group(
        'Wrapper options',
        'Standard options that require mapping files into the container',
    )
    g_wrap.add_argument(
        '-w',
        '--work-dir',
        action='store',
        type=os.path.abspath,
        help='path where intermediate results should be stored',
    )
    g_wrap.add_argument(
        '--output-spaces',
        nargs='*',
        help="""\
Standard and non-standard spaces to resample anatomical and functional images to. \
Standard spaces may be specified by the form \
``<TEMPLATE>[:res-<resolution>][:cohort-<label>][...]``, where ``<TEMPLATE>`` is \
a keyword (valid keywords: %s) or path pointing to a user-supplied template, and \
may be followed by optional, colon-separated parameters. \
Non-standard spaces (valid keywords: %s) imply specific orientations and sampling \
grids. \
Important to note, the ``res-*`` modifier does not define the resolution used for \
the spatial normalization."""
        % (
            ', '.join('"%s"' % s for s in TF_TEMPLATES),
            ', '.join(NONSTANDARD_REFERENCES),
        ),
    )

    g_wrap.add_argument(
        '--fs-license-file',
        metavar='PATH',
        type=IsFile,
        default=os.getenv('FS_LICENSE', None),
        help='Path to FreeSurfer license key file. Get it (for free) by registering'
        ' at https://surfer.nmr.mgh.harvard.edu/registration.html',
    )
    g_wrap.add_argument(
        '--fs-subjects-dir',
        metavar='PATH',
        type=os.path.abspath,
        help='Path to existing Infant FreeSurfer subjects directory to reuse. ',
    )
    g_wrap.add_argument(
        '--config-file',
        metavar='PATH',
        type=os.path.abspath,
        help='Use pre-generated configuration file. Values in file will be overridden '
        'by command-line arguments.',
    )
    g_wrap.add_argument(
        '--use-plugin',
        metavar='PATH',
        action='store',
        default=None,
        type=os.path.abspath,
        help='nipype plugin configuration file',
    )
    g_wrap.add_argument(
        '--bids-database-dir',
        metavar='PATH',
        type=os.path.abspath,
        help='Path to an existing PyBIDS database folder, for faster indexing '
        '(especially useful for large datasets).',
    )
    g_wrap.add_argument(
        '--segmentation-atlases-dir',
        metavar='PATH',
        type=os.path.abspath,
        help='Directory containing prelabeled segmentations to use for JointLabelFusion.',
    )
    g_wrap.add_argument(
        '-d',
        '--derivatives',
        nargs='+',
        metavar='PATH',
        action=ToDict,
        help='Search PATH(s) for pre-computed derivatives.',
    )
    g_wrap.add_argument(
        '--bids-filter-file',
        metavar='PATH',
        type=os.path.abspath,
        help='Filter file',
    )
    g_wrap.add_argument(
        '--deriv-filter-file',
        metavar='PATH',
        type=os.path.abspath,
        help='Filter file',
    )

    # Developer patch/shell options
    g_dev = parser.add_argument_group(
        'Developer options', 'Tools for testing and debugging nibabies'
    )
    g_dev.add_argument(
        '--patch',
        nargs='+',
        metavar='PACKAGE=PATH',
        action=ToDict,
        help='local repository to use within container',
    )
    g_dev.add_argument(
        '--shell',
        action='store_true',
        help='open shell in image instead of running nibabies',
    )
    g_dev.add_argument(
        '--config',
        metavar='PATH',
        action='store',
        type=os.path.abspath,
        help='Use custom nipype.cfg file',
    )
    g_dev.add_argument(
        '-e',
        '--env',
        action='append',
        nargs=2,
        metavar=('ENV_VAR', 'value'),
        help='Set custom environment variable within container',
    )
    g_dev.add_argument(
        '-u',
        '--user',
        action='store',
        help='Run container as a given user/uid. Additionally, group/gid can be'
        'assigned, (i.e., --user <UID>:<GID>)',
    )
    g_dev.add_argument(
        '--network',
        action='store',
        help='Run container with a different network driver '
        '("none" to simulate no internet connection)',
    )
    g_dev.add_argument('--no-tty', action='store_true', help='Run docker without TTY flag -it')

    return parser


def main():
    """Entry point"""

    parser = get_parser()
    # Capture additional arguments to pass inside container
    opts, unknown_args = parser.parse_known_args()

    if opts.version:
        print('nibabies wrapper {!s}'.format(__version__))
        return

    # Set help if no directories set
    if opts.help or not all((opts.service, opts.bids_dir, opts.output_dir)):
        parser.print_help()
        return

    container = ContainerManager(opts.service, image=opts.image)
    check = container.check_install()
    if check < 1:
        if check == -1:
            print(
                'nibabies: Could not find %s command... Is it installed?' % opts.service,
            )
        else:
            print(
                "nibabies: Make sure you have permission to run '%s'" % opts.service,
            )
        return 1

    if not container.check_image(opts.image):
        resp = 'Y'
        if opts.service == 'singularity':
            print('Singularity image must already exist locally.')
            return 1
        try:
            resp = input(MISSING.format(opts.image))
        except KeyboardInterrupt:
            print()
            return 1
        if resp not in ('y', 'Y', ''):
            return 0
        print('Downloading. This may take a while...')

    # Warn on low memory allocation
    mem_total = container.check_memory(opts.image)
    if mem_total == -1:
        print(
            'Could not detect memory capacity of Docker container.\n'
            'Do you have permission to run docker?'
        )
        return 1
    if '--reports-only' not in unknown_args and mem_total < 8000:
        print(
            'Warning: <8GB of RAM is available within your environment.\n'
            'Some parts of nibabies may fail to complete.'
        )
        if '--mem_mb' not in unknown_args:
            resp = 'N'
            try:
                resp = input('Continue anyway? [y/N]')
            except KeyboardInterrupt:
                print()
                return 1
            if resp not in ('y', 'Y', ''):
                return 0

    container.set_version()

    if opts.service == 'docker':
        if not opts.no_tty:
            container.add_cmd('-it')
        if opts.user:
            container.add_cmd(('-u', opts.user))
        if opts.network:
            container.add_cmd('--network=%s' % opts.network)

    # Patch working repositories into installed package directories
    if opts.patch:
        for pkg, repo_path in opts.patch.items():
            container.add_mount(repo_path, os.sep.join((PKG_PATH, pkg)))

    if opts.env:
        for envvar in opts.env:
            container.add_envvar(tuple(envvar))

    if opts.fs_license_file:
        container.add_mount(opts.fs_license_file, '/opt/freesurfer/license.txt')

    main_args = []
    if opts.bids_dir:
        container.add_mount(opts.bids_dir, '/data')
        main_args.append('/data')
    if opts.output_dir:
        if not os.path.exists(opts.output_dir):
            # create it before the container does
            os.makedirs(opts.output_dir)
        container.add_mount(opts.output_dir, '/out', read_only=False)
        main_args.append('/out')
    main_args.append(opts.analysis_level)

    if opts.fs_subjects_dir:
        container.add_mount(opts.fs_subjects_dir, '/opt/subjects', read_only=False)
        unknown_args.extend(['--fs-subjects-dir', '/opt/subjects'])

    if opts.config_file:
        container.add_mount(opts.config_file, '/tmp/config.toml', read_only=False)
        unknown_args.extend(['--config-file', '/tmp/config.toml'])

    if opts.segmentation_atlases_dir:
        container.add_mount(opts.segmentation_atlases_dir, '/opt/segmentations')
        unknown_args.extend(['--segmentation-atlases-dir', '/opt/segmentations'])
    if opts.bids_filter_file:
        container.add_mount(opts.bids_filter_file, '/opt/bids_filters.json')
        unknown_args.extend(['--bids-filter-file', '/opt/bids_filters.json'])
    if opts.deriv_filter_file:
        container.add_mount(opts.deriv_filter_file, '/opt/derivative_filters.json')
        unknown_args.extend(['--deriv-filter-file', '/opt/derivative_filters.json'])
    # Patch derivatives for searching
    if opts.derivatives:
        deriv_args = ['--derivatives']
        for deriv, deriv_path in opts.derivatives.items():
            deriv_target = '/deriv/%s' % deriv
            container.add_mount(deriv_path, deriv_target)
            deriv_args.append('='.join([deriv, deriv_target]))
        unknown_args.extend(deriv_args)

    # Check that work_dir is not a child of bids_dir
    if opts.work_dir and opts.bids_dir:
        if is_in_directory(opts.work_dir, opts.bids_dir):
            print(
                'The selected working directory is a subdirectory of the input BIDS folder. '
                'Please modify the output path.'
            )
            return 1

        if not os.path.exists(opts.work_dir):
            # create it before the container does
            os.makedirs(opts.work_dir)
        container.add_mount(opts.work_dir, '/scratch', read_only=False)
        unknown_args.extend(['-w', '/scratch'])

    if opts.config:
        container.add_mount('opts.config', '/home/fmriprep/.nipype/nipype.cfg')

    if opts.use_plugin:
        container.add_mount(opts.use_plugin, '/tmp/plugin.yml')
        unknown_args.extend(['--use-plugin', '/tmp/plugin.yml'])

    if opts.bids_database_dir:
        container.add_mount(opts.bids_database_dir, '/tmp/bids_db', read_only=False)
        unknown_args.extend(['--bids-database-dir', '/tmp/bids_db'])

    if opts.output_spaces:
        spaces = []
        for space in opts.output_spaces:
            if space.split(':')[0] not in (TF_TEMPLATES + NONSTANDARD_REFERENCES):
                tpl = os.path.basename(space)
                if not tpl.startswith('tpl-'):
                    raise RuntimeError('Custom template %s requires a `tpl-` prefix' % tpl)
                target = '/home/fmriprep/.cache/templateflow/' + tpl
                container.add_mount(os.path.abspath(space), target)
                spaces.append(tpl[4:])
            else:
                spaces.append(space)
        unknown_args.extend(['--output-spaces'] + spaces)

    if opts.shell:
        if opts.service == 'docker':
            container.add_cmd('--entrypoint=bash')
        elif opts.service == 'singularity':
            # replace default "run" command
            container.command[1] = 'shell'

    container.image = opts.image
    # after this, all call to ``container.add_cmd``
    # will be for nibabies arguments
    container.finalize_container_cmd()

    # Override help and version to describe underlying program
    # Respects '-i' flag, so will retrieve information from any image
    # if opts.help:
    #     container.add_cmd("-h")
    #     targethelp = subprocess.check_output(container.command).decode()
    #     print(merge_help(parser.format_help(), targethelp))
    #     return 0
    # elif opts.version:
    #     # Get version to be run and exit
    #     container.add_cmd("--version")
    #     ret = subprocess.run(container.command)
    #     return ret.returncode

    if not opts.shell:
        container.add_cmd(main_args)
        container.add_cmd(unknown_args)

    print('RUNNING: ' + ' '.join(container.command))
    ret = subprocess.run(container.command)
    if ret.returncode:
        print('nibabies: Please report errors to %s' % __bugreports__)
    return ret.returncode


if __name__ == '__main__':
    if '__main__.py' in sys.argv[0]:
        from . import __name__ as module

        sys.argv[0] = '%s -m %s' % (sys.executable, module)
    sys.exit(main())

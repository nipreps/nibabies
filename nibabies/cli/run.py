#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""NiBabies runner."""
from .. import config

EXITCODE: int = -1


def main():
    """Entry point."""
    import atexit
    import gc
    import os
    import sys
    from pathlib import Path

    from ..utils.bids import write_bidsignore, write_derivative_description
    from .parser import parse_args
    from .workflow import build_boilerplate, build_workflow

    _cwd = os.getcwd()
    # Revert OMP_NUM_THREADS + other runtime set environment variables
    atexit.register(config.restore_env)

    parse_args()

    # collect and submit telemetry information
    # if `--notrack` is specified, nothing is done.
    global EXITCODE
    if not config.execution.notrack:
        from nibabies.utils.telemetry import setup_migas

        setup_migas(init=True)
        atexit.register(migas_exit)

    if "participant" in config.workflow.analysis_level:
        _pool = None
        if config.nipype.plugin == "MultiProc":
            import multiprocessing as mp
            from concurrent.futures import ProcessPoolExecutor
            from contextlib import suppress

            # should drastically reduce VMS
            # see https://github.com/nipreps/mriqc/pull/984 for more details
            os.environ["OMP_NUM_THREADS"] = "1"

            with suppress(RuntimeError):
                mp.set_start_method("fork")
            gc.collect()

            _pool = ProcessPoolExecutor(
                max_workers=config.nipype.nprocs,
                initializer=config._process_initializer,
                initargs=(_cwd, config.nipype.omp_nthreads),
            )

        config_file = config.execution.work_dir / config.execution.run_uuid / "config.toml"
        config_file.parent.mkdir(exist_ok=True, parents=True)
        config.to_filename(config_file)

        # build the workflow within the same process
        # it still needs to be saved / loaded to be properly initialized
        retval = build_workflow(config_file)
        EXITCODE = retval['return_code']
        nibabies_wf = retval['workflow']

        # exit conditions:
        # - no workflow (--reports-only)
        # - retcode is not 0
        # - boilerplate only

        if nibabies_wf is None and not config.execution.reports_only:
            sys.exit(EXITCODE)

        if config.execution.write_graph:
            nibabies_wf.write_graph(graph2use="colored", format="svg", simple_form=True)

        if EXITCODE != 0:
            sys.exit(EXITCODE)

        # generate boilerplate
        build_boilerplate(nibabies_wf)
        if config.execution.boilerplate_only:
            sys.exit(EXITCODE)

        gc.collect()

        config.loggers.workflow.log(
            15,
            "\n".join(["nibabies config:"] + ["\t\t%s" % s for s in config.dumps().splitlines()]),
        )
        config.loggers.workflow.log(25, "nibabies started!")

        # Hack MultiProc's pool to reduce VMS
        _plugin = config.nipype.get_plugin()
        if _pool:
            from nipype.pipeline.plugins.multiproc import MultiProcPlugin

            multiproc = MultiProcPlugin(plugin_args=config.nipype.plugin_args)
            multiproc.pool = _pool
            _plugin = {"plugin": multiproc}

        gc.collect()
        try:
            nibabies_wf.run(**_plugin)
        except Exception as e:
            config.loggers.workflow.critical("nibabies failed: %s", e)
            EXITCODE = 1
            raise
        else:
            config.loggers.workflow.log(25, "nibabies finished successfully!")
            # Bother users with the boilerplate only iff the workflow went okay.
            boiler_file = config.execution.nibabies_dir / "logs" / "CITATION.md"
            if boiler_file.exists():
                if config.environment.exec_env in (
                    "singularity",
                    "docker",
                    "nibabies-docker",
                ):
                    boiler_file = Path("<OUTPUT_PATH>") / boiler_file.relative_to(
                        config.execution.output_dir
                    )
                config.loggers.workflow.log(
                    25,
                    "Works derived from this nibabies execution should include the "
                    f"boilerplate text found in {boiler_file}.",
                )

            if config.workflow.run_reconall:
                from niworkflows.utils.misc import _copy_any
                from templateflow import api

                dseg_tsv = str(api.get("fsaverage", suffix="dseg", extension=[".tsv"]))
                _copy_any(dseg_tsv, str(config.execution.nibabies_dir / "desc-aseg_dseg.tsv"))
                _copy_any(dseg_tsv, str(config.execution.nibabies_dir / "desc-aparcaseg_dseg.tsv"))
        # errno = 0
        finally:
            from pkg_resources import resource_filename as pkgrf

            from ..reports.core import generate_reports

            # Generate reports phase
            generate_reports(
                config.execution.participant_label,
                config.execution.session_id,
                config.execution.nibabies_dir,
                config.execution.run_uuid,
                config=pkgrf("nibabies", "data/reports-spec.yml"),
                packagename="nibabies",
            )
            write_derivative_description(config.execution.bids_dir, config.execution.nibabies_dir)
            write_bidsignore(config.execution.nibabies_dir)


def migas_exit() -> None:
    """
    Send a final crumb to the migas server signaling if the run successfully completed
    This function should be registered with `atexit` to run at termination.
    """
    import sys

    from nibabies.utils.telemetry import send_breadcrumb

    global EXITCODE
    migas_kwargs = {'status': 'C'}
    # `sys` will not have these attributes unless an error has been handled
    if hasattr(sys, 'last_type'):
        migas_kwargs = {
            'status': 'F',
            'status_desc': 'Finished with error(s)',
            'error_type': sys.last_type,
            'error_desc': sys.last_value,
        }
    elif EXITCODE != 0:
        migas_kwargs.update(
            {
                'status': 'F',
                'status_desc': f'Completed with exitcode {EXITCODE}',
            }
        )
    else:
        migas_kwargs['status_desc'] = 'Success'

    send_breadcrumb(**migas_kwargs)


if __name__ == "__main__":
    raise RuntimeError(
        "Please `pip install` this and run via the commandline interfaces, `nibabies <command>`"
    )

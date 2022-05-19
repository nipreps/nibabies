from itertools import product
from pathlib import Path

from niworkflows.reports.core import Report as _Report
from pkg_resources import resource_filename as pkgrf


class Report(_Report):
    # niworkflows patch to preserve `out_filename` even if subject_id is present
    def __init__(
        self,
        out_dir,
        run_uuid,
        config=None,
        out_filename=None,
        packagename=None,
        reportlets_dir=None,
        subject_id=None,
    ):
        self.root = Path(reportlets_dir or out_dir)

        # Initialize structuring elements
        self.sections = []
        self.errors = []
        self.out_dir = Path(out_dir)
        self.run_uuid = run_uuid
        self.packagename = packagename
        self.subject_id = subject_id
        if subject_id is not None:
            self.subject_id = subject_id[4:] if subject_id.startswith("sub-") else subject_id
            # ensure set output filename is preserved
            if not out_filename:
                out_filename = f"sub-{self.subject_id}.html"

        self.out_filename = out_filename or "report.html"

        # Default template from niworkflows
        self.template_path = Path(pkgrf("niworkflows", "reports/report.tpl"))
        self._load_config(Path(config or pkgrf("niworkflows", "reports/default.yml")))
        assert self.template_path.exists()

    # TODO: Upstream ``Report._load_config`` to niworkflows
    def _load_config(self, config):
        from yaml import safe_load as load

        settings = load(config.read_text())
        self.packagename = self.packagename or settings.get("package", None)

        # Removed from here: Appending self.packagename to self.root and self.out_dir
        # In this version, pass reportlets_dir and out_dir with nibabies in the path.

        if self.subject_id is not None:
            self.root = self.root / "sub-{}".format(self.subject_id)

        if "template_path" in settings:
            self.template_path = config.parent / settings["template_path"]

        self.index(settings["sections"])


#
# The following are the interface used directly by NiBabies
#


def run_reports(
    out_dir,
    subject_label,
    run_uuid,
    config=None,
    out_filename='report.html',
    reportlets_dir=None,
    packagename=None,
):
    """
    Run the reports.
    """
    return Report(
        out_dir,
        run_uuid,
        config=config,
        out_filename=out_filename,
        subject_id=subject_label,
        packagename=packagename,
        reportlets_dir=reportlets_dir,
    ).generate_report()


def generate_reports(
    subject_list,
    sessions_list,
    output_dir,
    run_uuid,
    config=None,
    work_dir=None,
    packagename=None,
):
    """Execute run_reports on a list of subjects."""
    reportlets_dir = None
    if work_dir is not None:
        reportlets_dir = Path(work_dir) / "reportlets"

    if sessions_list is None:
        sessions_list = [None]

    report_errors = []
    for subject_label, session in product(subject_list, sessions_list):
        html_report = f"sub-{subject_label}"
        if session:
            html_report += f"_ses-{session}"
        html_report += ".html"
        report_errors.append(
            run_reports(
                output_dir,
                subject_label,
                run_uuid,
                config=config,
                out_filename=html_report,
                packagename=packagename,
                reportlets_dir=reportlets_dir,
            )
        )

    errno = sum(report_errors)
    if errno:
        import logging

        logger = logging.getLogger("cli")
        error_list = ", ".join(
            "%s (%d)" % (subid, err) for subid, err in zip(subject_list, report_errors) if err
        )
        logger.error(
            "Preprocessing did not finish successfully. Errors occurred while processing "
            "data from participants: %s. Check the HTML reports for details.",
            error_list,
        )
    return errno

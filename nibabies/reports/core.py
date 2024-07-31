from pathlib import Path

from nireports.assembler.report import Report

from nibabies.data import load as load_data


def run_reports(
    out_dir,
    subject,
    run_uuid,
    session=None,
    out_filename=None,
    reportlets_dir=None,
    packagename=None,
):
    """
    Run the reports.
    """
    return Report(
        out_dir,
        run_uuid,
        subject=subject,
        session=session,
        bootstrap_file=load_data.readable('reports-spec.yml'),
        reportlets_dir=reportlets_dir,
    ).generate_report()


def generate_reports(
    sub_ses_list,
    output_dir,
    run_uuid,
    work_dir=None,
    packagename=None,
):
    """Execute run_reports on a list of subjects."""
    reportlets_dir = None
    if work_dir is not None:
        reportlets_dir = Path(work_dir) / 'reportlets'

    report_errors = []
    for subject, session in sub_ses_list:
        report_errors.append(
            run_reports(
                output_dir,
                subject,
                run_uuid,
                session=session,
                packagename=packagename,
                reportlets_dir=reportlets_dir,
            )
        )

    errno = sum(report_errors)
    if errno:
        import logging

        logger = logging.getLogger('cli')
        error_list = ', '.join(
            '%s (%d)' % (subid, err)
            for subid, err in zip(sub_ses_list, report_errors, strict=False)
            if err
        )
        logger.error(
            'Preprocessing did not finish successfully. Errors occurred while processing '
            'data from participants: %s. Check the HTML reports for details.',
            error_list,
        )
    return errno

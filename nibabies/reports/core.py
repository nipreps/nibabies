from pathlib import Path

from nireports.assembler.report import Report

from nibabies.data import load as load_data


def run_reports(
    out_dir,
    subject,
    run_uuid,
    session=None,
    bootstrap_file=None,
    out_filename='report.html',
    reportlets_dir=None,
):
    """
    Run the reports.
    """
    return Report(
        out_dir,
        run_uuid,
        subject=subject,
        session=session,
        bootstrap_file=load_data('reports-spec.yml'),
        reportlets_dir=reportlets_dir,
        out_filename=out_filename,
    ).generate_report()


def generate_reports(
    sub_ses_list,
    output_dir,
    run_uuid,
    *,
    work_dir=None,
    bootstrap_file=None,
    config_hash=None,
):
    """Execute run_reports on a list of subjects."""
    reportlets_dir = None
    if work_dir is not None:
        reportlets_dir = Path(work_dir) / 'reportlets'

    report_errors = []
    for subject, session in sub_ses_list:
        # Determine the output filename
        html_report = f'sub-{subject}'
        if session is not None:
            html_report += f'_ses-{session}'
        if config_hash is not None:
            html_report += f'_{config_hash}'
        html_report += '.html'

        report_errors.append(
            run_reports(
                output_dir,
                subject,
                run_uuid,
                bootstrap_file=bootstrap_file,
                reportlets_dir=reportlets_dir,
                out_filename=html_report,
                session=session,
            )
        )

    errno = sum(report_errors)
    if errno:
        import logging

        logger = logging.getLogger('cli')
        error_list = ', '.join(
            f'{subid} ({err})'
            for subid, err in zip(sub_ses_list, report_errors, strict=False)
            if err
        )
        logger.error(
            'Preprocessing did not finish successfully. Errors occurred while processing '
            'data from participants: %s. Check the HTML reports for details.',
            error_list,
        )
    return errno

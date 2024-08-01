from nibabel.optpkg import optional_package

from .. import __version__, config

migas = optional_package('migas')[0]


def setup_migas(init_ping: bool = True, exit_ping: bool = True) -> None:
    """
    Prepare the migas python client to communicate with a migas server.
    If ``init`` is ``True``, send an initial breadcrumb.
    """
    # generate session UUID from generated run UUID
    session_id = None
    if config.execution.run_uuid:
        session_id = config.execution.run_uuid.split('_', 1)[-1]

    migas.setup(session_id=session_id)
    if init_ping:
        # send initial status ping
        send_crumb(status='R', status_desc='workflow start')
    if exit_ping:
        from migas.error.nipype import node_execution_error

        migas.track_exit(
            'nipreps/nibabies',
            __version__,
            {'NodeExecutionError': node_execution_error},
        )


def send_crumb(**kwargs) -> dict:
    """
    Communicate with the migas telemetry server. This requires `migas.setup()` to be called.
    """
    return migas.add_breadcrumb('nipreps/nibabies', __version__, **kwargs)

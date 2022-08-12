import migas

from .. import __version__, config


def setup_migas(init_ping: bool = True) -> None:
    """
    Prepare the migas python client to communicate with a migas server.
    If `init_ping` is `True`, send an initial breadcrumb.
    """
    # generate session UUID from generated run UUID
    session_id = None
    if config.execution.run_uuid:
        session_id = config.execution.run_uuid.split('_', 1)[-1]

    migas.setup(session_id=session_id)
    if init_ping:
        # send initial status ping
        ping_migas(status='pending')


def ping_migas(*, status: str) -> dict:
    """
    Communicate with the migas telemetry server. This requires `migas.setup()` to be called.
    """
    res = migas.add_project("nipreps/nibabies", __version__, status=status)
    return res

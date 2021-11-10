from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

# Silence PyBIDS warning for extension entity behavior
# Can be removed once minimum PyBIDS dependency hits 0.14
try:
    import bids
    from packaging.version import Version

    if Version(bids.__version__) < Version("0.14.0"):
        bids.config.set_option("extension_initial_dot", True)
except (ImportError, ValueError):
    pass
else:
    del bids

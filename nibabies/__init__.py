from ._version import get_versions

__version__ = get_versions()["version"]
del get_versions

# Silence PyBIDS warning for extension entity behavior
# Can be removed once minimum PyBIDS dependency hits 0.14
try:
    import bids

    bids.config.set_option("extension_initial_dot", True)
except (ImportError, ValueError):
    pass
else:
    del bids

from niworkflows.interfaces.bids import DerivativesDataSink as _DDS


# TODO: Set default as default in niworkflows
class DerivativesDataSink(_DDS):
    out_path_base = ""

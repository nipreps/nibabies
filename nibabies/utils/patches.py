"Patch functions to allow 'flexible' connections"


def set_tf_resolution(template, sloppy):
    "'Dynamic' workaround to apply tweaks based on templates"
    if template in ('UNCInfant', 'MNIInfant'):
        from nipype.interfaces.base import Undefined
        return Undefined
    return sloppy + 1


def set_reg_resolution(template):
    if template in ('UNCInfant', 'MNIInfant'):
        return None
    from nipype.interfaces.base import Undefined
    return Undefined

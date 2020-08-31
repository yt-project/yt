import numpy as np
import traitlets


class YTPositionTrait(traitlets.TraitType):
    default_value = None
    info_text = "A position in code_length"

    def validate(self, obj, value):
        if isinstance(value, (tuple, list)):
            value = np.array(value)
        if hasattr(value, "in_units"):
            value = value.in_units("unitary").d
        if not isinstance(value, np.ndarray):
            self.error(obj, value)
        return value.astype("f4")


def ndarray_shape(*dimensions):
    # http://traittypes.readthedocs.io/en/latest/api_documentation.html
    def validator(trait, value):
        if value.shape != dimensions:
            raise traitlets.TraitError(
                "Expected an of shape %s and got and array with shape %s"
                % (dimensions, value.shape)
            )
        else:
            return value

    return validator


def ndarray_ro():
    def validator(trait, value):
        if value.flags["WRITEABLE"]:
            value = value.copy()
            value.flags["WRITEABLE"] = False
        return value

    return validator

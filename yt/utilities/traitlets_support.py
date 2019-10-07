import traitlets
import traittypes
import numpy as np
import yt.units.dimensions as dims

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

class YTDimensionfulTrait(traitlets.TraitType):
    dimensions = None
    default_value = None
    info_text = "A length"

    def __init__(self, default_value=None, allow_none=False, dimensions=dims.length, **kwargs):    
        self.dimensions = dimensions
        super().__init__(default_value=default_value, allow_none=allow_none, **kwargs)

    def validate(self, obj, value):
        # We assume that obj has a `ds` attribute
        if obj.ds is None:
            return value
        if isinstance(value, (tuple, list)):
            v1 = value
            value = obj.ds.quan(value[0], value[1])

        if not isinstance(value, np.ndarray):
            self.error(obj, value)
        if value.units.dimensions != self.dimensions:
            raise traitlets.TraitError('Expected dimensions of %s, got %s' % (self.dimensions, value.units.dimensions))
        return value

def ndarray_shape(*dimensions):
    # http://traittypes.readthedocs.io/en/latest/api_documentation.html
    def validator(trait, value):
        if value.shape != dimensions:
            raise traitlets.TraitError('Expected an of shape %s and got and array with shape %s' % (dimensions, value.shape))
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



import numpy as np

from yt import YTQuantity, YTArray
from yt.data_objects.static_output import Dataset
from yt.utilities.exceptions import YTInvalidArgumentType

def input_validator(_type_name):
    def validator(func):
        def wrapper(*args, **kwargs):
            bad_input = False
            data_category_1 = (list, np.ndarray, YTArray, YTQuantity)
            data_category_2 = (int, float, tuple, YTQuantity)

            if _type_name == 'disk':
                kwargs_list = ['fields', 'ds', 'field_parameters', 'data_source']

                if not isinstance(args[1], data_category_1) \
                        or not isinstance(args[2], data_category_1) \
                        or not isinstance(args[3], data_category_2) \
                        or not isinstance(args[4], data_category_2):
                    bad_input = True

                if len(args[1]) != 3 or len(args[2]) != 3 \
                    or (isinstance(args[3], YTQuantity) and len(args[3]) != 1) \
                    or (isinstance(args[4], YTQuantity) and len(args[4]) != 1):
                    bad_input = True

                from yt.data_objects.selection_data_containers import YTRegion
                for k in kwargs :
                    if (k == kwargs_list[0] and
                        not isinstance(kwargs[k], (list, np.ndarray))) \
                        or (k == kwargs_list[1] and not isinstance(kwargs[k], Dataset)) \
                        or (k == kwargs_list[2] and
                            not isinstance(kwargs[k], dict)) \
                        or (k == kwargs_list[3]
                            and not isinstance(kwargs[k], YTRegion)):
                        bad_input = True
                    if bad_input:
                        break

            if not bad_input:
                func(*args, **kwargs)
            else:
                raise YTInvalidArgumentType(_type_name)
        return wrapper
    return validator

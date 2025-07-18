# avoid deprecation warnings in numpy >= 2.0

import numpy as np

if hasattr(np, "trapezoid"):
    # np.trapz is deprecated in numpy 2.0 in favor of np.trapezoid
    trapezoid = np.trapezoid
else:
    trapezoid = np.trapz  # type: ignore[assignment] # noqa: NPY201

import re
from typing import Tuple

import matplotlib


def version_tuple(version: str) -> Tuple[int, ...]:
    elems = version.split(".")
    if len(elems) > 3:
        elems = elems[:3]

    if not elems[-1].isnumeric():
        # allow alpha/beta/release candidate versions
        match = re.search(r"^\d+", elems[-1])
        if match is None:
            elems.pop()
        else:
            elems[-1] = match.group()

    return tuple(int(_) for _ in elems)


MPL_VERSION = version_tuple(matplotlib.__version__)

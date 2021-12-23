import sys
from typing import List, Optional, Tuple

if sys.version_info >= (3, 10):
    from typing import TypeAlias
else:
    from typing_extensions import TypeAlias


FieldDescT: TypeAlias = Tuple[str, Tuple[str, List[str], Optional[str]]]
KnownFieldsT: TypeAlias = Tuple[FieldDescT, ...]

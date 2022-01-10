from typing import List, Optional, Tuple

FieldDescT = Tuple[str, Tuple[str, List[str], Optional[str]]]
KnownFieldsT = Tuple[FieldDescT, ...]

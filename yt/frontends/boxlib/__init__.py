"""

import sys

from .. import amrex

sys.module[__name__] = __import__('amrex')

sys.module[__name__.data_structures] = __import__('amrex.data_structures')

sys.module[__name__.tests] = __import__('amrex.tests')

sys.module[__name__.fields] = __import__('amrex.fields')

sys.module[__name__.io] = __import__('amrex.io')

"""
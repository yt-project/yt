import os
import sys
from difflib import unified_diff

from nose.plugins import Plugin

from yt.utilities.exceptions import YTAmbiguousFieldName

ROOT = os.path.abspath(__file__)

diff_file = open("diff.txt", mode="w")


class AmbiguousResolvePlugin(Plugin):
    """Plugin that takes no command-line arguments"""

    name = "auto-update-ambiguous-fields"
    enabled = True
    _my_version = None

    def configure(self, options, conf):
        pass

    def options(self, parser, env=None):
        pass

    def addError(self, test, err: sys.exc_info):
        try:
            test_path = test.context.__file__
        except Exception:
            test_path = os.path.join(ROOT, test.context.__module__.replace(".", "/"))
        t, v, tb = err
        if t is not YTAmbiguousFieldName:
            return

        import traceback

        # print('TYPE:', t)
        # print('VALUE:', v)
        # print('TRACEBACK:', tb)

        ambiguous_fname = v.fname
        ok = False
        for ft, _ in traceback.walk_tb(tb):
            # line = (
            #     open(ft.f_code.co_filename).readlines()[ft.f_lineno].replace("\n", "")
            # )
            # print(f"{ft.f_code.co_filename}:{ft.f_code.co_firstlineno} {line}")
            if test_path == ft.f_code.co_filename:
                ft_err = ft
                ok = True

        if not ok:
            return

        # Now, ft contains the current frame,
        # need to correct the test!
        # print(f"error from {ft_err.f_code.co_filename}:{ft_err.f_code.co_firstlineno}:{ft_err.f_lineno}")

        with open(test_path) as f:
            lines = f.readlines()

        lineno = ft_err.f_lineno - 1

        suggested_ftype = " OR ".join(v.possible_ftypes)

        corrected = lines.copy()
        corrected[lineno] = (
            corrected[lineno]
            .replace(
                '"%s"' % ambiguous_fname, f'("{suggested_ftype}", "{ambiguous_fname}")'
            )
            .replace(
                "'%s'" % ambiguous_fname, f"('{suggested_ftype}', '{ambiguous_fname}')"
            )
        )

        rel_path = os.path.relpath(test_path)
        diff_file.writelines(
            unified_diff(lines, corrected, fromfile=rel_path, tofile=rel_path)
        )
        diff_file.flush()


# run(
#     argv=["nosetests", "yt"],
#     plugins=[AmbiguousResolvePlugin()]
# )

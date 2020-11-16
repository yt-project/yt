import os
import sys
import traceback
from difflib import unified_diff

from nose.plugins import Plugin
from nose.plugins.plugintest import run

from yt.utilities.exceptions import YTAmbiguousFieldName

ROOT = os.path.join(os.path.abspath(__file__), "..")

diff_file = open("diff.txt", mode="w")


class AmbiguousResolvePlugin(Plugin):
    """Plugin that takes no command-line arguments"""

    name = "ambiguous-solver"
    enabled = True

    def configure(self, options, conf):
        pass

    def options(self, parser, env={}):
        pass

    def addError(self, test, err: sys.exc_info):
        try:
            test_path = test.context.__file__
        except AttributeError:
            test_path = os.path.join(ROOT, test.context.__module__.replace(".", "/"))
        t, v, tb = err
        if t is not YTAmbiguousFieldName:
            return

        ambiguous_fname = v.fname
        ok = False
        for ft, _ in traceback.walk_tb(tb):
            # line = (
            #     open(ft.f_code.co_filename, "r")
            #     .readlines()[ft.f_lineno]
            #     .replace("\n", "")
            # )
            # print(f"{ft.f_code.co_filename}:{ft.f_code.co_firstlineno} {line}")
            if test_path == ft.f_code.co_filename:
                ft_err = ft
                ok = True

        if not ok:
            return

        # Now, ft contains the current frame, need to correct the test!
        with open(test_path, "r") as f:
            lines = f.readlines()

        lineno = ft_err.f_lineno - 1

        suggested_ftype = v.possible_ftypes[-1]
        corrected = lines.copy()
        corrected[lineno] = (
            corrected[lineno]
            .replace(
                f'"{ambiguous_fname}"', f'("{suggested_ftype}", "{ambiguous_fname}")',
            )
            .replace(
                f"'{ambiguous_fname}'", f"('{suggested_ftype}', '{ambiguous_fname}')",
            )
        )

        rel_path = os.path.relpath(test_path, ROOT)
        diff_file.writelines(
            unified_diff(lines, corrected, fromfile=rel_path, tofile=rel_path)
        )
        diff_file.flush()


run(argv=["nosetests", "yt"], plugins=[AmbiguousResolvePlugin()])

import subprocess
import sys

from docutils.parsers.rst import Directive


def setup(app):
    app.add_directive("yt_showfields", ShowFields)
    setup.app = app
    setup.config = app.config
    setup.confdir = app.confdir

    retdict = dict(version="1.0", parallel_read_safe=True, parallel_write_safe=True)

    return retdict


class ShowFields(Directive):
    required_arguments = 0
    optional_arguments = 0
    parallel_read_safe = True
    parallel_write_safe = True

    def run(self):
        rst_file = self.state_machine.document.attributes["source"]
        lines = subprocess.check_output(
            [sys.executable, "./helper_scripts/show_fields.py"]
        )
        lines = lines.decode("utf8")
        lines = lines.split("\n")
        self.state_machine.insert_input(lines, rst_file)
        return []

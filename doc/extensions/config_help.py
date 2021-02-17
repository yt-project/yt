import re
import subprocess

from docutils import statemachine
from docutils.parsers.rst import Directive


def setup(app):
    app.add_directive("config_help", GetConfigHelp)
    setup.app = app
    setup.config = app.config
    setup.confdir = app.confdir

    retdict = dict(version="1.0", parallel_read_safe=True, parallel_write_safe=True)

    return retdict


class GetConfigHelp(Directive):
    required_arguments = 1
    optional_arguments = 0
    final_argument_whitespace = True

    def run(self):
        rst_file = self.state_machine.document.attributes["source"]
        data = (
            subprocess.check_output(self.arguments[0].split(" ") + ["-h"])
            .decode("utf8")
            .split("\n")
        )
        ind = next(
            (i for i, val in enumerate(data) if re.match(r"\s{0,3}\{.*\}\s*$", val))
        )
        lines = [".. code-block:: none", ""] + data[ind + 1 :]
        self.state_machine.insert_input(
            statemachine.string2lines("\n".join(lines)), rst_file
        )
        return []

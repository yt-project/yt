import tempfile
import os
import glob
import base64
import shutil
import subprocess
from sphinx.util.compat import Directive
from docutils import nodes


class PythonScriptDirective(Directive):
    """Execute an inline python script and display images.

    This uses exec to execute an inline python script, copies
    any images produced by the script, and embeds them in the document
    along with the script.

    """
    required_arguments = 0
    optional_arguments = 0
    has_content = True

    def run(self):
        cwd = os.getcwd()
        tmpdir = tempfile.mkdtemp()
        os.chdir(tmpdir)

        # Construct script from cell content
        content = "\n".join(self.content)
        with open("temp.py", "w") as f:
            f.write(content)

        # Use sphinx logger?
        print ""
        print content
        print ""

        subprocess.call(['python', 'temp.py'])

        text = ''
        for im in sorted(glob.glob("*.png")):
            text += get_image_tag(im)

        code = content

        literal = nodes.literal_block(code, code)
        literal['language'] = 'python'

        attributes = {'format': 'html'}
        img_node = nodes.raw('', text, **attributes)

        # clean up
        os.chdir(cwd)
        shutil.rmtree(tmpdir, True)

        return [literal, img_node]


def setup(app):
    app.add_directive('python-script', PythonScriptDirective)
    setup.app = app
    setup.config = app.config
    setup.confdir = app.confdir

    retdict = dict(
        version='0.1',
        parallel_read_safe=True,
        parallel_write_safe=True
    )

    return retdict


def get_image_tag(filename):
    with open(filename, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read())
        return '<img src="data:image/png;base64,%s" width="600"><br>' \
            % encoded_string

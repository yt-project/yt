import errno
import glob
import os
import shutil
import subprocess
import tempfile
import time
import uuid

from docutils import nodes
from docutils.parsers.rst import Directive


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

        rst_file = self.state_machine.document.attributes["source"]
        rst_dir = os.path.abspath(os.path.dirname(rst_file))

        image_dir, image_rel_dir = make_image_dir(setup, rst_dir)

        # Construct script from cell content
        content = "\n".join(self.content)
        with open("temp.py", "w") as f:
            f.write(content)

        # Use sphinx logger?
        uid = uuid.uuid4().hex[:8]
        print("")
        print(f">> Contents of the script: {uid}")
        print(content)
        print("")

        start = time.time()
        subprocess.call(["python", "temp.py"])
        print(f">> The execution of the script {uid} took {time.time() - start:f} s")
        text = ""
        for im in sorted(glob.glob("*.png")):
            text += get_image_tag(im, image_dir, image_rel_dir)

        code = content

        literal = nodes.literal_block(code, code)
        literal["language"] = "python"

        attributes = {"format": "html"}
        img_node = nodes.raw("", text, **attributes)

        # clean up
        os.chdir(cwd)
        shutil.rmtree(tmpdir, True)

        return [literal, img_node]


def setup(app):
    app.add_directive("python-script", PythonScriptDirective)
    setup.app = app
    setup.config = app.config
    setup.confdir = app.confdir

    retdict = dict(version="0.1", parallel_read_safe=True, parallel_write_safe=True)

    return retdict


def get_image_tag(filename, image_dir, image_rel_dir):
    my_uuid = uuid.uuid4().hex
    shutil.move(filename, image_dir + os.path.sep + my_uuid + filename)
    relative_filename = image_rel_dir + os.path.sep + my_uuid + filename
    return f'<img src="{relative_filename}" width="600"><br>'


def make_image_dir(setup, rst_dir):
    image_dir = setup.app.builder.outdir + os.path.sep + "_images"
    rel_dir = os.path.relpath(setup.confdir, rst_dir)
    image_rel_dir = rel_dir + os.path.sep + "_images"
    thread_safe_mkdir(image_dir)
    return image_dir, image_rel_dir


def thread_safe_mkdir(dirname):
    try:
        os.makedirs(dirname)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise

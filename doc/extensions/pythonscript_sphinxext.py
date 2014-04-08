from sphinx.util.compat import Directive
from subprocess import Popen,PIPE
from docutils.parsers.rst import directives
from docutils import nodes
import os, glob, base64

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
        # Construct script from cell content
        content = "\n".join(self.content)
        with open("temp.py", "w") as f:
            f.write(content)

        # Use sphinx logger?
        print ""
        print content
        print ""

        codeproc = Popen(['python', 'temp.py'], stdout=PIPE)
        out = codeproc.stdout.read()

        images = sorted(glob.glob("*.png"))
        fns = []
        text = ''
        for im in images:
            text += get_image_tag(im)
            os.remove(im)
            
        os.remove("temp.py")

        code = content

        literal = nodes.literal_block(code,code)
        literal['language'] = 'python'

        attributes = {'format': 'html'}
        img_node = nodes.raw('', text, **attributes)
        
        return [literal, img_node]

def setup(app):
    app.add_directive('python-script', PythonScriptDirective)
    setup.app = app
    setup.config = app.config
    setup.confdir = app.confdir

def get_image_tag(filename):
    with open(filename, "rb") as image_file:
        encoded_string = base64.b64encode(image_file.read())
        return '<img src="data:image/png;base64,%s" width="600"><br>' % encoded_string

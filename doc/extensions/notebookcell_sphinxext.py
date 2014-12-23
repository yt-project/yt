import os
import shutil
import io
import tempfile
from sphinx.util.compat import Directive
from docutils.parsers.rst import directives
from IPython.nbformat import current
from notebook_sphinxext import \
    notebook_node, visit_notebook_node, depart_notebook_node, \
    evaluate_notebook, make_image_dir, write_notebook_output


class NotebookCellDirective(Directive):
    """Insert an evaluated notebook cell into a document

    This uses runipy and nbconvert to transform an inline python
    script into html suitable for embedding in a Sphinx document.
    """
    required_arguments = 0
    optional_arguments = 1
    has_content = True
    option_spec = {'skip_exceptions': directives.flag}

    def run(self):
        # check if raw html is supported
        if not self.state.document.settings.raw_enabled:
            raise self.warning('"%s" directive disabled.' % self.name)

        cwd = os.getcwd()
        tmpdir = tempfile.mkdtemp()
        os.chdir(tmpdir)

        rst_file = self.state_machine.document.attributes['source']
        rst_dir = os.path.abspath(os.path.dirname(rst_file))

        image_dir, image_rel_dir = make_image_dir(setup, rst_dir)

        # Construct notebook from cell content
        content = "\n".join(self.content)
        with open("temp.py", "w") as f:
            f.write(content)

        convert_to_ipynb('temp.py', 'temp.ipynb')

        skip_exceptions = 'skip_exceptions' in self.options

        evaluated_text, resources = evaluate_notebook(
            'temp.ipynb', skip_exceptions=skip_exceptions)

        evaluated_text = write_notebook_output(
            resources, image_dir, image_rel_dir, evaluated_text)

        # create notebook node
        attributes = {'format': 'html', 'source': 'nb_path'}
        nb_node = notebook_node('', evaluated_text, **attributes)
        (nb_node.source, nb_node.line) = \
            self.state_machine.get_source_and_line(self.lineno)

        # clean up
        os.chdir(cwd)
        shutil.rmtree(tmpdir, True)

        return [nb_node]

def setup(app):
    setup.app = app
    setup.config = app.config
    setup.confdir = app.confdir

    app.add_node(notebook_node,
                 html=(visit_notebook_node, depart_notebook_node))

    app.add_directive('notebook-cell', NotebookCellDirective)

    retdict = dict(
        version='0.1',
        parallel_read_safe=True,
        parallel_write_safe=True
    )

    return retdict

def convert_to_ipynb(py_file, ipynb_file):
    with io.open(py_file, 'r', encoding='utf-8') as f:
        notebook = current.reads(f.read(), format='py')
    with io.open(ipynb_file, 'w', encoding='utf-8') as f:
        current.write(notebook, f, format='ipynb')

import errno
import os
import shutil
import string
import re
import tempfile
import uuid
from sphinx.util.compat import Directive
from docutils import nodes
from docutils.parsers.rst import directives
from IPython.config import Config
from IPython.nbconvert import html, python
from IPython.nbformat import current as nbformat
from runipy.notebook_runner import NotebookRunner, NotebookError

class NotebookDirective(Directive):
    """Insert an evaluated notebook into a document

    This uses runipy and nbconvert to transform a path to an unevaluated notebook
    into html suitable for embedding in a Sphinx document.
    """
    required_arguments = 1
    optional_arguments = 1
    option_spec = {'skip_exceptions': directives.flag}
    final_argument_whitespace = True

    def run(self): # check if there are spaces in the notebook name
        nb_path = self.arguments[0]
        if ' ' in nb_path: raise ValueError(
            "Due to issues with docutils stripping spaces from links, white "
            "space is not allowed in notebook filenames '{0}'".format(nb_path))
        # check if raw html is supported
        if not self.state.document.settings.raw_enabled:
            raise self.warning('"%s" directive disabled.' % self.name)

        cwd = os.getcwd()
        tmpdir = tempfile.mkdtemp()
        os.chdir(tmpdir)

        # get path to notebook
        nb_filename = self.arguments[0]
        nb_basename = os.path.basename(nb_filename)
        rst_file = self.state_machine.document.attributes['source']
        rst_dir = os.path.abspath(os.path.dirname(rst_file))
        nb_abs_path = os.path.abspath(os.path.join(rst_dir, nb_filename))

        # Move files around.
        rel_dir = os.path.relpath(rst_dir, setup.confdir)
        dest_dir = os.path.join(setup.app.builder.outdir, rel_dir)
        dest_path = os.path.join(dest_dir, nb_basename)

        image_dir, image_rel_dir = make_image_dir(setup, rst_dir)

        # Ensure desination build directory exists
        thread_safe_mkdir(os.path.dirname(dest_path))

        # Copy unevaluated notebook
        shutil.copyfile(nb_abs_path, dest_path)

        # Construct paths to versions getting copied over
        dest_path_eval = string.replace(dest_path, '.ipynb', '_evaluated.ipynb')
        dest_path_script = string.replace(dest_path, '.ipynb', '.py')
        rel_path_eval = string.replace(nb_basename, '.ipynb', '_evaluated.ipynb')
        rel_path_script = string.replace(nb_basename, '.ipynb', '.py')

        # Create python script vesion
        script_text = nb_to_python(nb_abs_path)
        f = open(dest_path_script, 'w')
        f.write(script_text.encode('utf8'))
        f.close()

        skip_exceptions = 'skip_exceptions' in self.options

        ret = evaluate_notebook(
            nb_abs_path, dest_path_eval, skip_exceptions=skip_exceptions)

        try:
            evaluated_text, resources = ret
            evaluated_text = write_notebook_output(
                resources, image_dir, image_rel_dir, evaluated_text)
        except ValueError:
            # This happens when a notebook raises an unhandled exception
            evaluated_text = ret

        # Create link to notebook and script files
        link_rst = "(" + \
                   formatted_link(nb_basename) + "; " + \
                   formatted_link(rel_path_eval) + "; " + \
                   formatted_link(rel_path_script) + \
                   ")"

        self.state_machine.insert_input([link_rst], rst_file)

        # create notebook node
        attributes = {'format': 'html', 'source': 'nb_path'}
        nb_node = notebook_node('', evaluated_text, **attributes)
        (nb_node.source, nb_node.line) = \
            self.state_machine.get_source_and_line(self.lineno)

        # add dependency
        self.state.document.settings.record_dependencies.add(nb_abs_path)

        # clean up
        os.chdir(cwd)
        shutil.rmtree(tmpdir, True)

        return [nb_node]


class notebook_node(nodes.raw):
    pass

def nb_to_python(nb_path):
    """convert notebook to python script"""
    exporter = python.PythonExporter()
    output, resources = exporter.from_filename(nb_path)
    return output

def nb_to_html(nb_path):
    """convert notebook to html"""
    c = Config({'ExtractOutputPreprocessor':{'enabled':True}})

    exporter = html.HTMLExporter(template_file='full', config=c)
    notebook = nbformat.read(open(nb_path), 'json')
    output, resources = exporter.from_notebook_node(notebook)
    header = output.split('<head>', 1)[1].split('</head>',1)[0]
    body = output.split('<body>', 1)[1].split('</body>',1)[0]

    # http://imgur.com/eR9bMRH
    header = header.replace('<style', '<style scoped="scoped"')
    header = header.replace('body {\n  overflow: visible;\n  padding: 8px;\n}\n',
                            '')
    header = header.replace("code,pre{", "code{")

    # Filter out styles that conflict with the sphinx theme.
    filter_strings = [
        'navbar',
        'body{',
        'alert{',
        'uneditable-input{',
        'collapse{',
    ]

    filter_strings.extend(['h%s{' % (i+1) for i in range(6)])

    line_begin = [
        'pre{',
        'p{margin'
    ]

    filterfunc = lambda x: not any([s in x for s in filter_strings])
    header_lines = filter(filterfunc, header.split('\n'))

    filterfunc = lambda x: not any([x.startswith(s) for s in line_begin])
    header_lines = filter(filterfunc, header_lines)

    header = '\n'.join(header_lines)

    # concatenate raw html lines
    lines = ['<div class="ipynotebook">']
    lines.append(header)
    lines.append(body)
    lines.append('</div>')
    return '\n'.join(lines), resources

def evaluate_notebook(nb_path, dest_path=None, skip_exceptions=False):
    # Create evaluated version and save it to the dest path.
    notebook = nbformat.read(open(nb_path), 'json')
    nb_runner = NotebookRunner(notebook, pylab=False)
    try:
        nb_runner.run_notebook(skip_exceptions=skip_exceptions)
    except NotebookError as e:
        print('')
        print(e)
        # Return the traceback, filtering out ANSI color codes.
        # http://stackoverflow.com/questions/13506033/filtering-out-ansi-escape-sequences
        return "Notebook conversion failed with the " \
               "following traceback: \n%s" % \
            re.sub(r'\\033[\[\]]([0-9]{1,2}([;@][0-9]{0,2})*)*[mKP]?', '',
                   str(e))

    if dest_path is None:
        dest_path = 'temp_evaluated.ipynb'
    nbformat.write(nb_runner.nb, open(dest_path, 'w'), 'json')
    ret = nb_to_html(dest_path)
    if dest_path is 'temp_evaluated.ipynb':
        os.remove(dest_path)
    return ret

def formatted_link(path):
    return "`%s <%s>`__" % (os.path.basename(path), path)

def visit_notebook_node(self, node):
    self.visit_raw(node)

def depart_notebook_node(self, node):
    self.depart_raw(node)

def setup(app):
    setup.app = app
    setup.config = app.config
    setup.confdir = app.confdir

    app.add_node(notebook_node,
                 html=(visit_notebook_node, depart_notebook_node))

    app.add_directive('notebook', NotebookDirective)

    retdict = dict(
        version='0.1',
        parallel_read_safe=True,
        parallel_write_safe=True
    )

    return retdict

def make_image_dir(setup, rst_dir):
    image_dir = setup.app.builder.outdir + os.path.sep + '_images'
    rel_dir = os.path.relpath(setup.confdir, rst_dir)
    image_rel_dir = rel_dir + os.path.sep + '_images'
    thread_safe_mkdir(image_dir)
    return image_dir, image_rel_dir

def write_notebook_output(resources, image_dir, image_rel_dir, evaluated_text):
    my_uuid = uuid.uuid4().hex

    for output in resources['outputs']:
        new_name = image_dir + os.path.sep + my_uuid + output
        new_relative_name = image_rel_dir + os.path.sep + my_uuid + output
        evaluated_text = evaluated_text.replace(output, new_relative_name)
        with open(new_name, 'wb') as f:
            f.write(resources['outputs'][output])
    return evaluated_text

def thread_safe_mkdir(dirname):
    try:
        os.makedirs(dirname)
    except OSError as e:
        if e.errno != errno.EEXIST:
            raise
        pass

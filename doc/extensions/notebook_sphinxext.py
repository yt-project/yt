import os, shutil, string, glob
from sphinx.util.compat import Directive
from docutils import nodes
from docutils.parsers.rst import directives
from IPython.nbconvert import html, python
from runipy.notebook_runner import NotebookRunner

class NotebookDirective(Directive):
    """Insert an evaluated notebook into a document

    This uses runipy and nbconvert to transform a path to an unevaluated notebook
    into html suitable for embedding in a Sphinx document.
    """
    required_arguments = 1
    optional_arguments = 1
    option_spec = {'skip_exceptions' : directives.flag}

    def run(self):
        # check if raw html is supported
        if not self.state.document.settings.raw_enabled:
            raise self.warning('"%s" directive disabled.' % self.name)

        # get path to notebook
        source_dir = os.path.dirname(
            os.path.abspath(self.state.document.current_source))
        nb_basename = os.path.basename(self.arguments[0])
        rst_file = self.state_machine.document.attributes['source']
        rst_dir = os.path.abspath(os.path.dirname(rst_file))
        nb_abs_path = os.path.join(rst_dir, nb_basename)

        # Move files around.
        rel_dir = os.path.relpath(rst_dir, setup.confdir)
        rel_path = os.path.join(rel_dir, nb_basename)
        dest_dir = os.path.join(setup.app.builder.outdir, rel_dir)
        dest_path = os.path.join(dest_dir, nb_basename)

        if not os.path.exists(dest_dir):
            os.makedirs(dest_dir)

        # Copy unevaluated script
        try:
            shutil.copyfile(nb_abs_path, dest_path)
        except IOError:
            raise RuntimeError("Unable to copy notebook to build destination.")

        dest_path_eval = string.replace(dest_path, '.ipynb', '_evaluated.ipynb')
        dest_path_script = string.replace(dest_path, '.ipynb', '.py')
        rel_path_eval = string.replace(nb_basename, '.ipynb', '_evaluated.ipynb')
        rel_path_script = string.replace(nb_basename, '.ipynb', '.py')

        # Create python script vesion
        unevaluated_text = nb_to_html(nb_abs_path)
        script_text = nb_to_python(nb_abs_path)
        f = open(dest_path_script, 'w')
        f.write(script_text.encode('utf8'))
        f.close()

        skip_exceptions = 'skip_exceptions' in self.options

        try:
            evaluated_text = evaluate_notebook(nb_abs_path, dest_path_eval,
                                               skip_exceptions=skip_exceptions)
        except:
            # bail
            return []

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

        # clean up png files left behind by notebooks.
        png_files = glob.glob("*.png")
        fits_files = glob.glob("*.fits")
        h5_files = glob.glob("*.h5")
        for file in png_files:
            os.remove(file)

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
    exporter = html.HTMLExporter(template_file='full')
    output, resources = exporter.from_filename(nb_path)
    header = output.split('<head>', 1)[1].split('</head>',1)[0]
    body = output.split('<body>', 1)[1].split('</body>',1)[0]

    # http://imgur.com/eR9bMRH
    header = header.replace('<style', '<style scoped="scoped"')
    header = header.replace('body{background-color:#ffffff;}\n', '')
    header = header.replace('body{background-color:white;position:absolute;'
                            'left:0px;right:0px;top:0px;bottom:0px;'
                            'overflow:visible;}\n', '')
    header = header.replace('body{margin:0;'
                            'font-family:"Helvetica Neue",Helvetica,Arial,'
                            'sans-serif;font-size:13px;line-height:20px;'
                            'color:#000000;background-color:#ffffff;}', '')
    header = header.replace('\na{color:#0088cc;text-decoration:none;}', '')
    header = header.replace(
        'a:focus{color:#005580;text-decoration:underline;}', '')
    header = header.replace(
        '\nh1,h2,h3,h4,h5,h6{margin:10px 0;font-family:inherit;font-weight:bold;'
        'line-height:20px;color:inherit;text-rendering:optimizelegibility;}'
        'h1 small,h2 small,h3 small,h4 small,h5 small,'
        'h6 small{font-weight:normal;line-height:1;color:#999999;}'
        '\nh1,h2,h3{line-height:40px;}\nh1{font-size:35.75px;}'
        '\nh2{font-size:29.25px;}\nh3{font-size:22.75px;}'
        '\nh4{font-size:16.25px;}\nh5{font-size:13px;}'
        '\nh6{font-size:11.049999999999999px;}\nh1 small{font-size:22.75px;}'
        '\nh2 small{font-size:16.25px;}\nh3 small{font-size:13px;}'
        '\nh4 small{font-size:13px;}', '')
    header = header.replace('body {\n  overflow: visible;\n  padding: 8px;\n}\n', '')
    header = header.replace('background-color:#ffffff;', '', 1)

    # concatenate raw html lines
    lines = ['<div class="ipynotebook">']
    lines.append(header)
    lines.append(body)
    lines.append('</div>')
    return '\n'.join(lines)

def evaluate_notebook(nb_path, dest_path=None, skip_exceptions=False):
    # Create evaluated version and save it to the dest path.
    # Always use --pylab so figures appear inline
    # perhaps this is questionable?
    nb_runner = NotebookRunner(nb_in=nb_path, pylab=True)
    nb_runner.run_notebook(skip_exceptions=skip_exceptions)
    if dest_path is None:
        dest_path = 'temp_evaluated.ipynb'
    nb_runner.save_notebook(dest_path)
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

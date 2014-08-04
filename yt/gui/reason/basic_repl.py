"""
A read-eval-print-loop.  This code was released in the CherryPy project
initially, but has been heavily modified and again re-released in compliance
with the terms of its original license.



"""

#-----------------------------------------------------------------------------
# Copyright (c) 2013, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
#-----------------------------------------------------------------------------

import codeop
import inspect
import re
import json
import sys
import traceback
from yt.extern.six.moves import StringIO

class ProgrammaticREPL(object):
    stopped = False
    debug = False
    def __init__(self, locals=None):
        self.executed_cell_texts = []
        self.locals = {}
        if locals:
            self.locals.update(locals)
        self.buffer = []
        # Nominally at this point we could populate our namespace with widgets
        # or other useful functions.  We aren't really ready for that just yet.

    def evaluate_cell(self, input):
        result = []
        for line in input.split("\n"):
            r = self.push(line)
            if r is not None: result.append(r)
        # To catch if people leave a little bit at the end for us
        if r is None:
            r = self.push("\n")
            if r is not None: result.append(r)
        # Usually they have \n already
        return "".join(result)
    
    def push(self, line):
        """Push 'line' and return exec results (None if more input needed)."""
        if line == "help":
            return "Type help(object) for help about object."
        if line == "help()":
            return "You cannot call help() without an argument."
        
        self.buffer.append(line)
        source = "\n".join(self.buffer)
        
        try:
            code = codeop.compile_command(source, "<HTTP input>", 'single')
        except (OverflowError, SyntaxError, ValueError):
            self.buffer = []
            return traceback.format_exc()
        
        if code is None:
            # More lines needed.
            return None
        
        self.buffer = []
        return self.execute(code)
    
    def execute(self, code):
        """Execute the given code in self.locals and return any stdout/sterr."""
        out = StringIO()
        oldout = sys.stdout
        olderr = sys.stderr
        sys.stdout = sys.stderr = out
        try:
            try:
                exec code in self.locals
            except:
                result = traceback.format_exc()
            else:
                result = out.getvalue()
        finally:
            sys.stdout = oldout
            sys.stderr = olderr
        out.close()
        return result
    
    def dir(self, line):
        """Examine a partial line and provide attr list of final expr."""
        line = re.split(r"\s", line)[-1].strip()
        # Support lines like "thing.attr" as "thing.", because the browser
        # may not finish calculating the partial line until after the user
        # has clicked on a few more keys.
        line = ".".join(line.split(".")[:-1])
        try:
            result = eval("dir(%s)" % line, {}, self.locals)
        except:
            return []
        return result
    
    def doc(self, line):
        """Examine a partial line and provide sig+doc of final expr."""
        line = re.split(r"\s", line)[-1].strip()
        # Support lines like "func(text" as "func(", because the browser
        # may not finish calculating the partial line until after the user
        # has clicked on a few more keys.
        line = "(".join(line.split("(")[:-1])
        try:
            result = eval(line, {}, self.locals)
            try:
                if isinstance(result, type):
                    func = result.__init__
                else:
                    func = result
                args, varargs, varkw, defaults = inspect.getargspec(func)
            except TypeError:
                if callable(result):
                    doc = getattr(result, "__doc__", "") or ""
                    return "%s\n\n%s" % (line, doc)
                return None
        except:
            return None
        
        if args and args[0] == 'self':
            args.pop(0)
        missing = object()
        defaults = defaults or []
        defaults = ([missing] * (len(args) - len(defaults))) + list(defaults)
        arglist = []
        for a, d in zip(args, defaults):
            if d is missing:
                arglist.append(a)
            else:
                arglist.append("%s=%s" % (a, d))
        if varargs:
            arglist.append("*%s" % varargs)
        if varkw:
            arglist.append("**%s" % varkw)
        doc = getattr(result, "__doc__", "") or ""
        return "%s(%s)\n%s" % (line, ", ".join(arglist), doc)

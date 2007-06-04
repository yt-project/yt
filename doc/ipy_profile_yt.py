"""YT mode for IPython

Start ipython in shell mode by invoking "ipython -p yt"

"""

from IPython import ipapi, Prompts
import os,textwrap

import yt
import yt.lagos as lagos
import yt.raven as raven
import yt.fido as fido
import yt.enki as enki

# The import below effectively obsoletes your old-style ipythonrc[.ini],
# so consider yourself warned!

import ipy_defaults

def runBrowse(shell, args):
    #return shell
    retval=fido.runBrowse()
    if isinstance(retval, lagos.EnzoRun):
        print "Setting mr = %s" % (retval)
        print "Setting mh = %s" % (None)
        shell.api.user_ns['mr'] = retval
        shell.api.user_ns['mh'] = None
    elif isinstance(retval, lagos.EnzoHierarchy):
        print "Setting mr = %s" % (retval.run)
        print "Setting mh = %s" % (retval)
        shell.api.user_ns['mr'] = retval.run
        shell.api.user_ns['mh'] = retval
    elif isinstance(retval, lagos.EnzoParameterFile):
        hh = retval.promote()
        print "Setting mr = %s" % (retval.run)
        print "Setting mh = %s (promoted)" % (hh)
        shell.api.user_ns['mr'] = retval.run
        shell.api.user_ns['mh'] = hh
    else:
        print "Choose a real one next time!"
        return None

def runHippo(shell, args):
    if not os.getenv("DISPLAY"):
        print "DISPLAY not set; not starting hippo"
        return None
    if not shell.api.user_ns['mh'] and len(args) < 1:
        print "mh needs to be set, or a hierarchy has to be provided"
        return None
    if len(args) > 0:
        hh = args[0]
    else:
        hh = shell.api.user_ns['mh']
    app = None
    canvas = None
    if shell.api.user_ns.has_key("_hippoapp"):
        app = shell.api.user_ns['_hippoapp']
        canvas = -1
    hip = raven.EnzoHippo(hh, app=app, canvas=canvas)
    shell.api.user_ns['_hippoapp'] = hip.app
    shell.api.user_ns['hippos'].append(hip)
    shell.api.user_ns['ch'] = hip
    return hip

def logOff(shell, args):
    shell.api.ex("yt.logger.logging.disable(50)")

def logOn(shell, args):
    shell.api.ex("yt.logger.logging.disable(30)")

def logAll(shell, args):
    shell.api.ex("yt.logger.logging.disable(0)")

def addSlice(shell, args):
    pass
    

def myinputprompt(self, cont):
    #Prompts.PromptColors.set_active_scheme("Linux")
    ip=self.api
    colors = Prompts.PromptColors[""].colors
    count = str(len(ip.user_ns["_ih"]))
    run = ip.user_ns['mr']
    h = ip.user_ns['mh']
    if cont:
        return "%s%s%s: " % ( \
            colors.in_prompt2,
            "."*(4+len(str(run))+2+len(str(h))+len(count)+1),
            colors.normal)
    else:
        return "%sIn [%s%s%s|%s|%s%s%s]: %s" % ( \
            colors.in_prompt,
            colors.in_prompt2,
            ip.user_ns['mr'],
            colors.in_prompt,
            ip.user_ns['mh'],
            colors.out_number,
            count,
            colors.in_prompt,
            colors.normal)

def main():
    ip = ipapi.get()
    o = ip.options
    # autocall to "full" mode (smart mode is default, I like full mode)
    
    o.autocall = 2
    
    # Jason Orendorff's path class is handy to have in user namespace
    # if you are doing shell-like stuff
    try:
        ip.ex("from path import path" )
    except ImportError:
        pass
    
    ip.ex('import os')
    ip.ex("def up(): os.chdir('..')")
    ip.ex('import yt')
    ip.ex('import yt.lagos as lagos')
    ip.ex('import yt.raven as raven')
    ip.ex('import yt.fido as fido')
    ip.ex('import yt.enki as enki')

    ip.ex('ch = None')
    ip.ex('mh = None')
    ip.ex('mr = None')
    ip.ex('hippos = []')
        
    # Get pysh-like prompt for all profiles. 
    
    #o.prompt_in1 = prompt_format % run
    ip.set_hook("generate_prompt", myinputprompt)
    o.prompt_in2= '\C_Green|\C_LightGreen\D\C_Green> '
    o.prompt_out= '<\#> '
    
    from IPython import Release

    import sys
    # I like my banner minimal.
    #o.banner = "Py %s IPy %s\n" % (sys.version.split('\n')[0],Release.version)
    o.banner = "\nyt (http://www.slac.stanford.edu/~mturk/yt_doc/)\n" + \
               "  yt, lagos, fido, raven, enki available\n\n"

    ip.expose_magic("gr",runBrowse)
    ip.expose_magic("gh",runHippo)
    ip.expose_magic("logoff",logOff)
    ip.expose_magic("logon",logOn)
    ip.expose_magic("logall",logAll)
    
    # make 'd' an alias for ls -F
    
    ip.magic('alias d ls -F --color=auto')
    
    # Make available all system commands through "rehashing" immediately. 
    # You can comment these lines out to speed up startup on very slow 
    # machines, and to conserve a bit of memory. Note that pysh profile does this
    # automatically
    ip.IP.default_option('cd','-q')
    

    o.prompts_pad_left="1"
    # Remove all blank lines in between prompts, like a normal shell.
    o.separate_in="0"
    o.separate_out="0"
    o.separate_out2="0"
    
    # now alias all syscommands
    
    db = ip.db
    
    syscmds = db.get("syscmdlist",[] )
    if not syscmds:
        print textwrap.dedent("""
        System command list not initialized, probably the first run...
        running %rehashx to refresh the command list. Run %rehashx
        again to refresh command list (after installing new software etc.)
        """)
        ip.magic('rehashx')
        syscmds = db.get("syscmdlist")
    for cmd in syscmds:
        #print "al",cmd
        noext, ext = os.path.splitext(cmd)
        ip.IP.alias_table[noext] = (0,cmd)
    extend_shell_behavior(ip)

def extend_shell_behavior(ip):

    # Instead of making signature a global variable tie it to IPSHELL.
    # In future if it is required to distinguish between different
    # shells we can assign a signature per shell basis
    ip.IP.__sig__ = 0xa005
    # mark the IPSHELL with this signature
    ip.IP.user_ns['__builtins__'].__dict__['__sig__'] = ip.IP.__sig__

    from IPython.Itpl import ItplNS
    from IPython.genutils import shell
    # utility to expand user variables via Itpl
    # xxx do something sensible with depth?
    ip.IP.var_expand = lambda cmd, lvars=None, depth=2: \
        str(ItplNS(cmd.replace('#','\#'), ip.IP.user_ns, get_locals()))

    def get_locals():
        """ Substituting a variable through Itpl deep inside the IPSHELL stack
            requires the knowledge of all the variables in scope upto the last
            IPSHELL frame. This routine simply merges all the local variables
            on the IPSHELL stack without worrying about their scope rules
        """
        import sys
        # note lambda expression constitues a function call
        # hence fno should be incremented by one
        getsig = lambda fno: sys._getframe(fno+1).f_globals \
                             ['__builtins__'].__dict__['__sig__']
        getlvars = lambda fno: sys._getframe(fno+1).f_locals
        # trackback until we enter the IPSHELL
        frame_no = 1
        sig = ip.IP.__sig__
        fsig = ~sig
        while fsig != sig :
            try:
                fsig = getsig(frame_no)
            except (AttributeError, KeyError):
                frame_no += 1
            except ValueError:
                # stack is depleted
                # call did not originate from IPSHELL
                return {}
        first_frame = frame_no
        # walk further back until we exit from IPSHELL or deplete stack
        try:
            while(sig == getsig(frame_no+1)):
                frame_no += 1
        except (AttributeError, KeyError, ValueError):
            pass
        # merge the locals from top down hence overriding
        # any re-definitions of variables, functions etc.
        lvars = {}
        for fno in range(frame_no, first_frame-1, -1):
            lvars.update(getlvars(fno))
        #print '\n'*5, first_frame, frame_no, '\n', lvars, '\n'*5 #dbg
        return lvars

    def _runlines(lines):
        """Run a string of one or more lines of source.

        This method is capable of running a string containing multiple source
        lines, as if they had been entered at the IPython prompt.  Since it
        exposes IPython's processing machinery, the given strings can contain
        magic calls (%magic), special shell access (!cmd), etc."""

        # We must start with a clean buffer, in case this is run from an
        # interactive IPython session (via a magic, for example).
        ip.IP.resetbuffer()
        lines = lines.split('\n')
        more = 0
        command = ''
        for line in lines:
            # skip blank lines so we don't mess up the prompt counter, but do
            # NOT skip even a blank line if we are in a code block (more is
            # true)
            # if command is not empty trim the line
            if command != '' :
                line = line.strip()
            # add the broken line to the command
            if line and line[-1] == '\\' :
                command += line[0:-1] + ' '
                more = True
                continue
            else :
                # add the last (current) line to the command
                command += line
                if command or more:
                    more = ip.IP.push(ip.IP.prefilter(command,more))
                    command = ''
                    # IPython's runsource returns None if there was an error
                    # compiling the code.  This allows us to stop processing right
                    # away, so the user gets the error message at the right place.
                    if more is None:
                        break
        # final newline in case the input didn't have it, so that the code
        # actually does get executed
        if more:
            ip.IP.push('\n')

    ip.IP.runlines = _runlines

main()

# encoding: utf-8
"""
Enable pyglet to be used interacive by setting PyOS_InputHook.

Authors
-------

* Nicolas P. Rougier
* Fernando Perez
"""

#-----------------------------------------------------------------------------
#  Copyright (C) 2008-2011  The IPython Development Team
#
#  Distributed under the terms of the BSD License.  The full license is in
#  the file COPYING, distributed as part of this software.
#-----------------------------------------------------------------------------

# This has been modified from the Pyglet and GLUT event hooks to work with
# glfw.

#-----------------------------------------------------------------------------
# Imports
#-----------------------------------------------------------------------------

import os
import sys
import time

#-----------------------------------------------------------------------------
# Platform-dependent imports and functions
#-----------------------------------------------------------------------------

if os.name == 'posix':
    import select

    def stdin_ready():
        infds, outfds, erfds = select.select([sys.stdin],[],[],0)
        if infds:
            return True
        else:
            return False

elif sys.platform == 'win32':
    import msvcrt

    def stdin_ready():
        return msvcrt.kbhit()

#-----------------------------------------------------------------------------
# Code
#-----------------------------------------------------------------------------

def create_inputhook_glfw(mgr, render_loop):
    """Run the GLFW event loop by processing pending events only.

    This keeps processing pending events until stdin is ready.  After
    processing all pending events, a call to time.sleep is inserted.  This is
    needed, otherwise, CPU usage is at 100%.  This sleep time should be tuned
    though for best performance.
    """
    def inputhook_glfw():
        # We need to protect against a user pressing Control-C when IPython is
        # idle and this is running. We trap KeyboardInterrupt and pass.
        import cyglfw3 as glfw
        try:
            t = glfw.GetTime()
            while not stdin_ready():
                render_loop.next()

                used_time = glfw.GetTime() - t
                if used_time > 10.0:
                    # print 'Sleep for 1 s'  # dbg
                    time.sleep(1.0)
                elif used_time > 0.1:
                    # Few GUI events coming in, so we can sleep longer
                    # print 'Sleep for 0.05 s'  # dbg
                    time.sleep(0.05)
                else:
                    # Many GUI events coming in, so sleep only very little
                    time.sleep(0.001)
        except KeyboardInterrupt:
            pass
        return 0
    return inputhook_glfw

from IPython.lib.inputhook import inputhook_manager, InputHookBase

@inputhook_manager.register('glfw')
class GLFWInputHook(InputHookBase):
    def enable(self, app=None):
        """Enable event loop integration with GLFW.

        Parameters
        ----------
        app : ignored
           Ignored, it's only a placeholder to keep the call signature of all
           gui activation methods consistent, which simplifies the logic of
           supporting magics.

        Notes
        -----
        This methods sets the ``PyOS_InputHook`` for GLFW, which allows
        GLFW to integrate with terminal based applications like
        IPython.

        """
        inputhook_glfw = create_inputhook_glfw(self.manager, app)
        self.manager.set_inputhook(inputhook_glfw)
        return

# encoding: utf-8
"""
GLFW-based input hook, inspired by the IPython pyglet inputhook.

Original Authors
----------------

* Nicolas P. Rougier
* Fernando Perez

"""

# ----------------------------------------------------------------------------
# Copyright (c) 2017, yt Development Team.
#
# Distributed under the terms of the Modified BSD License.
#
# The full license is in the file COPYING.txt, distributed with this software.
# ----------------------------------------------------------------------------

# This is a part of the experimental Interactive Data Visualization 

import cyglfw3 as glfw
from IPython.terminal.pt_inputhooks import register

class InputHookGLFW:
    def __init__(self, rendering_context, long_used_time = 10.0,
                 used_time = 0.1, sleep_time = 1.0, short_sleep = 0.05,
                 long_sleep = 0.001):
        self.rendering_context = rendering_context
        self.used_time = used_time
        self.sleep_time = sleep_time
        self.short_sleep = short_sleep
        self.long_sleep = long_sleep

    def __call__(self, context):
        # We need to protect against a user pressing Control-C when IPython is
        # idle and this is running. We trap KeyboardInterrupt and pass.
        rc = self.rendering_context
        try:
            t = glfw.GetTime()
            while not context.input_is_ready():
                try:
                    next(self.rendering_context)
                except StopIteration:
                    break

                used_time = glfw.GetTime() - t
                if 1:
                    pass
                elif used_time > self.long_used_time:
                    # print 'Sleep for 1 s'  # dbg
                    time.sleep(self.sleep_time)
                    # Update our sleep times, in case they changed
                elif used_time > self.used_time:
                    # Few GUI events coming in, so we can sleep longer
                    # print 'Sleep for 0.05 s'  # dbg
                    time.sleep(self.long_sleep)
                else:
                    # Many GUI events coming in, so sleep only very little
                    time.sleep(self.short_sleep)
        except KeyboardInterrupt:
            pass
        return 0


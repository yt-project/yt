#!/usr/bin/python
# -*- coding: utf-8 -*-
#
# progressbar  - Text progress bar library for Python.
# Copyright (c) 2005 Nilton Volpato
#
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA

"""Main ProgressBar class."""

from __future__ import division

import math
import os
import signal
import sys
import time
import uuid

try:
    from fcntl import ioctl
    from array import array
    import termios
except ImportError:
    pass

import widgets

# Test to see if we are in an IPython session.
ipython = None
for key in ['KernelApp','IPKernelApp']:
  try:
    ipython = get_ipython().config[key]['parent_appname']
  except (NameError, KeyError):
    pass

ipython_notebook_css = """
td.pb_widget {
    width: auto;
}
td.pb_widget_fill {
    width: 100%;
}
table.pb {
    font-family: monospace;
    border: 0;
    margin: 0;
}
table.pb tr { border: 0; }
table.pb td {
    white-space: nowrap;
    border: 0;
}
div.pb {
    border: 1px solid #ddd;
    border-radius: 3px;
}
div.pb_bar {
    height: 1.5em;
}
""".replace('\n', ' ')

class UnknownLength: pass


class ProgressBar(object):
    """The ProgressBar class which updates and prints the bar.

    A common way of using it is like:
    >>> pbar = ProgressBar().start()
    >>> for i in range(100):
    ...    # do something
    ...    pbar.update(i+1)
    ...
    >>> pbar.finish()

    You can also use a ProgressBar as an iterator:
    >>> progress = ProgressBar()
    >>> for i in progress(some_iterable):
    ...    # do something
    ...

    Since the progress bar is incredibly customizable you can specify
    different widgets of any type in any order. You can even write your own
    widgets! However, since there are already a good number of widgets you
    should probably play around with them before moving on to create your own
    widgets.

    The term_width parameter represents the current terminal width. If the
    parameter is set to an integer then the progress bar will use that,
    otherwise it will attempt to determine the terminal width falling back to
    80 columns if the width cannot be determined.

    When implementing a widget's update method you are passed a reference to
    the current progress bar. As a result, you have access to the
    ProgressBar's methods and attributes. Although there is nothing preventing
    you from changing the ProgressBar you should treat it as read only.

    Useful methods and attributes include (Public API):
     - currval: current progress (0 <= currval <= maxval)
     - maxval: maximum (and final) value
     - finished: True if the bar has finished (reached 100%)
     - start_time: the time when start() method of ProgressBar was called
     - seconds_elapsed: seconds elapsed since start_time and last call to
                        update
     - percentage(): progress in percent [0..100]
    """

    __slots__ = ('currval', 'fd', 'finished', 'last_update_time',
                 'left_justify', 'maxval', 'next_update', 'num_intervals',
                 'poll', 'seconds_elapsed', 'signal_set', 'start_time',
                 'term_width', 'update_interval', 'widgets', '_time_sensitive',
                 '__iterable', 'attr', 'html_written', 'uuid')

    _DEFAULT_MAXVAL = 100
    _DEFAULT_TERMSIZE = 80
    _DEFAULT_WIDGETS = [widgets.Percentage, widgets.Bar]

    def __init__(self, maxval=None, widgets=None, term_width=None, poll=1,
                 left_justify=True, fd=sys.stdout, attr={}):
        """Initializes a progress bar with sane defaults."""

        # Don't share a reference with any other progress bars
        if widgets is None:
            widgets = [widget() for widget in self._DEFAULT_WIDGETS]

        self.maxval = maxval
        self.widgets = widgets
        self.fd = fd
        self.left_justify = left_justify

        self.signal_set = False
        if term_width is not None:
            self.term_width = term_width
        else:
            try:
                self._handle_resize()
                signal.signal(signal.SIGWINCH, self._handle_resize)
                self.signal_set = True
            except (SystemExit, KeyboardInterrupt): raise
            except:
                self.term_width = self._env_size()

        self.__iterable = None
        self._update_widgets()
        self.currval = 0
        self.finished = False
        self.last_update_time = None
        self.poll = poll
        self.seconds_elapsed = 0
        self.start_time = None
        self.update_interval = 1
        self.attr = attr

        # Set flag so we only write out the HTML once,
        # then update with javascript
        self.html_written = False

        self.uuid = str(uuid.uuid4())

        # Install our CSS if we are in an IPython notebook
        if ipython == 'ipython-notebook':
            from IPython.display import Javascript, display
            display(Javascript('//%s\n$("head").append("<style>%s</style>")' %
                               (self.uuid,ipython_notebook_css)))
            
            # Also add a function that removes progressbar output from the cells
            js = '''
                  // %s -- used to remove this code blob in the end
                  IPython.OutputArea.prototype.cleanProgressBar = function(uuids) {
                      // filter by uuid-strings 
                      var myfilter = function(output) { 
                          var nuids = uuids.length;
                          for (var i=0; i<nuids; i++) {
                              if (output.hasOwnProperty('html')) {
                                  if (output.html.indexOf(uuids[i]) != -1) {
                                      return false;
                                  }
                              }
                              if (output.hasOwnProperty('javascript')) {
                                  if (output.javascript.indexOf(uuids[i]) != -1) {
                                      return false;
                                  }
                              }
                          }
                          // keep all others
                          return true;
                      };

                      // Filter the ouputs
                      this.outputs = this.outputs.filter(myfilter);
                };
                ''' % self.uuid
            display(Javascript(js))

    def __call__(self, iterable):
        """Use a ProgressBar to iterate through an iterable."""

        try:
            self.maxval = len(iterable)
        except:
            if self.maxval is None:
                self.maxval = UnknownLength

        self.__iterable = iter(iterable)
        return self


    def __iter__(self):
        return self


    def __next__(self):
        try:
            value = next(self.__iterable)
            if self.start_time is None: self.start()
            else: self.update(self.currval + 1)
            return value
        except StopIteration:
            self.finish()
            raise


    # Create an alias so that Python 2.x won't complain about not being
    # an iterator.
    next = __next__


    def _env_size(self):
        """Tries to find the term_width from the environment."""

        return int(os.environ.get('COLUMNS', self._DEFAULT_TERMSIZE)) - 1


    def _handle_resize(self, signum=None, frame=None):
        """Tries to catch resize signals sent from the terminal."""

        h, w = array('h', ioctl(self.fd, termios.TIOCGWINSZ, '\0' * 8))[:2]
        self.term_width = w


    def percentage(self):
        """Returns the progress as a percentage."""
        return self.currval * 100.0 / self.maxval

    percent = property(percentage)


    def _format_widgets(self):
        result = []
        expanding = []
        width = self.term_width

        for index, widget in enumerate(self.widgets):
            if isinstance(widget, widgets.WidgetHFill):
                result.append(widget)
                expanding.insert(0, index)
            else:
                widget = widgets.format_updatable(widget, self)
                result.append(widget)
                width -= len(widget)

        count = len(expanding)
        while count:
            portion = max(int(math.ceil(width * 1. / count)), 0)
            index = expanding.pop()
            count -= 1

            widget = result[index].update(self, portion)
            width -= len(widget)
            result[index] = widget

        return result


    def _format_line(self):
        """Joins the widgets and justifies the line."""

        widgets = ''.join(self._format_widgets())

        if self.left_justify: return widgets.ljust(self.term_width)
        else: return widgets.rjust(self.term_width)


    def _format_html(self):
      html = '<div class="pb" id="%s"><table class="pb ui-widget"><tr>\n' % self.uuid
      for widget in self.widgets:
        if isinstance(widget, widgets.WidgetHFill):
          td_class = 'pb_widget_fill'
        else:
          td_class = 'pb_widget'

        html += ('<td class="%s">' % td_class) + \
                widgets.format_updatable_html(widget, self) + \
                '</td>\n'
      html += '</tr></table><div>'
      return html


    def _need_update(self):
        """Returns whether the ProgressBar should redraw the line."""
        if self.currval >= self.next_update or self.finished: return True

        delta = time.time() - self.last_update_time
        return self._time_sensitive and delta > self.poll


    def _update_widgets(self):
        """Checks all widgets for the time sensitive bit."""

        self._time_sensitive = any(getattr(w, 'TIME_SENSITIVE', False)
                                    for w in self.widgets)


    def update(self, value=None, attr={}):
        """Updates the ProgressBar to a new value."""

        if value is not None and value is not UnknownLength:
            if (self.maxval is not UnknownLength
                and not 0 <= value <= self.maxval):

                raise ValueError('Value out of range')

            self.currval = value

        self.attr.update(attr)

        if not self._need_update(): return
        if self.start_time is None:
            raise RuntimeError('You must call "start" before calling "update"')

        now = time.time()
        self.seconds_elapsed = now - self.start_time
        self.next_update = self.currval + self.update_interval

        if ipython == 'ipython-notebook':
            if not self.html_written:
                # We have yet to display the HTML, do that first
                from IPython.display import HTML, display
                display(HTML(self._format_html()))
                self.html_written = True
            else:
                # The HTML has been written once, now update with JS
                from IPython.display import Javascript, display
                for widget in self.widgets:
                    js = widgets.updatable_js(widget, self)
                    if js:
                        display(Javascript(js))
        else:
            self.fd.write('\r' + self._format_line())
            self.fd.flush()

        self.last_update_time = now


    def start(self):
        """Starts measuring time, and prints the bar at 0%.

        It returns self so you can use it like this:
        >>> pbar = ProgressBar().start()
        >>> for i in range(100):
        ...    # do something
        ...    pbar.update(i+1)
        ...
        >>> pbar.finish()
        """

        if self.maxval is None:
            self.maxval = self._DEFAULT_MAXVAL

        self.num_intervals = max(100, self.term_width)
        self.next_update = 0

        if self.maxval is not UnknownLength:
            if self.maxval < 0: raise ValueError('Value out of range')
            self.update_interval = self.maxval / self.num_intervals


        self.start_time = self.last_update_time = time.time()
        self.html_written = False
        self.finished = False
        self.update(0)

        return self


    def finish(self):
        """Puts the ProgressBar bar in the finished state."""

        self.finished = True
        self.update(self.maxval)
        self.start_time = None

        # Clean up notebook stuff, quite differently from regular
        if not ipython == 'ipython-notebook':
            self.fd.write('\n')
        else:
            from IPython.display import Javascript, display
            # First delete the node that held the progress bar from the page
            js = """var element = document.getElementById('%s');
                    element.parentNode.removeChild(element);""" % self.uuid
            display(Javascript(js))

            # Then also remove its trace from the cell output (so it doesn't get
            # stored with the notebook). This needs to be done for all widgets as
            # well as for progressBar
            uuids = [str(self.uuid)]
            uuids += [w.uuid for w in self.widgets if isinstance(w, widgets.Widget)]
            display(Javascript('this.cleanProgressBar(%s)' % uuids))

        if self.signal_set:
            signal.signal(signal.SIGWINCH, signal.SIG_DFL)

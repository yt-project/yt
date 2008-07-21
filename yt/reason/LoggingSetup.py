"""
Logging facilities for Reason

Author: Matthew Turk <matthewturk@gmail.com>
Affiliation: KIPAC/SLAC/Stanford
Homepage: http://yt.enzotools.org/
License:
  Copyright (C) 2008 Matthew Turk.  All Rights Reserved.

  This file is part of yt.

  yt is free software; you can redistribute it and/or modify
  it under the terms of the GNU General Public License as published by
  the Free Software Foundation; either version 3 of the License, or
  (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU General Public License for more details.

  You should have received a copy of the GNU General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.
"""

from yt.reason import *

import logging, sys
import wx.lib.mixins.listctrl as listmix

class LoggingListBox(wx.ListCtrl, listmix.ListCtrlAutoWidthMixin):
    def __init__(self, parent, dataSource):
        wx.ListCtrl.__init__(self, parent,
                             style=wx.LC_REPORT|wx.LC_SINGLE_SEL|wx.LC_VIRTUAL)
        listmix.ListCtrlAutoWidthMixin.__init__(self)
        columns = dataSource.GetColumnHeaders()
        for col, text in enumerate(columns):
            self.InsertColumn(col, text)
        self.SetItemCount(dataSource.GetCount())
        self.dataSource = dataSource

    def UpdateCounts(self):
        c = self.dataSource.GetCount()
        self.SetItemCount(c)
        self.EnsureVisible(c-1)
        self.Refresh()

    def OnGetItemText(self, item, col):
        return self.dataSource.GetItem(item)[col]

    def OnGetItemAttr(self, item): return None
    def OnGetItemImage(self, item): return -1

class LoggingWindowBox(wx.Panel):
    def __init__(self, parent, id):
        wx.Panel.__init__(self, parent, id)
        self.__redirect_logging()
        self.__setup_controls()
        mylog.info("Welcome to Reason.")

    def __setup_controls(self):
        self.event_list = LoggingListBox(self, self.handler.output)
        self.handler.setTextBox(self.event_list)
        self.sizer = wx.BoxSizer(wx.VERTICAL)
        self.sizer.Add(self.event_list, 1, wx.EXPAND)
        self.SetSizer(self.sizer)
        self.Fit()

    def __redirect_logging(self):
        yt.logger.disable_stream_logging()
        self.handler = OutputWindowStream()
        self.handler.setFormatter(logging.Formatter(yt.logger.fstring))
        yt.logger.rootLogger.addHandler(self.handler)

class OutputWindowStream(logging.Handler):
    # Shamelessly borrowed from:
    #   http://lists.wxwidgets.org/pipermail/wxpython-users/2007-October/069288.html
    def __init__(self):
        logging.Handler.__init__(self)
        self.output = LoggingDataSource()
        self.buffer = []
        self.events = []

    def setTextBox(self, textbox):
        """
        textbox presents the instance of textctrl in the mainframe
        Method is called after frame is up and running
        """
        self.textbox = textbox

    def shouldFlush(self):
        if self.textbox:
            return True
        else:
            return False

    def emit(self, record):
        self.buffer.append(record)
        if self.shouldFlush():
            self.flush()

    def flush(self):
        """
        Thread safe flush method which schedule the
        textbox update after a wx main loop ends
        """
        for record in self.buffer:
            self.format(record)
            wx.CallAfter(self.output.AppendItem, record)
            wx.CallAfter(self.textbox.UpdateCounts)
            wx.CallAfter(wx.YieldIfNeeded)
        self.buffer = []

    def close(self):
        self.flush()
        self.textbox = None
        logging.Handler.close(self)

class LoggingDataSource(object):
    cols = ("Component","Level","Time","Message")
    attrs = ("name","levelname","asctime","message")
    def __init__(self):
        self.records = []

    def GetColumnHeaders(self):
        return self.cols

    def GetCount(self):
        return len(self.records)

    def GetItem(self, index):
        return [getattr(self.records[index],att) for att in self.attrs]

    def UpdateCache(self, start, end):
        pass

    def AppendItem(self, rec):
        self.records.append(rec)

"""Streaming HTTP uploads module.

This module extends the standard httplib and urllib2 objects so that
iterable objects can be used in the body of HTTP requests.

In most cases all one should have to do is call :func:`register_openers()`
to register the new streaming http handlers which will take priority over
the default handlers, and then you can use iterable objects in the body
of HTTP requests.

**N.B.** You must specify a Content-Length header if using an iterable object
since there is no way to determine in advance the total size that will be
yielded, and there is no way to reset an iterator.

Example usage:

>>> from StringIO import StringIO
>>> import urllib2, poster.streaminghttp

>>> opener = poster.streaminghttp.register_openers()

>>> s = "Test file data"
>>> f = StringIO(s)

>>> req = urllib2.Request("http://localhost:5000", f,
...                       {'Content-Length': str(len(s))})
"""


import http.client as http_client
import socket
import sys
import urllib


def request_has_data(req):
    if sys.version_info >= (3, 0, 0):
        return hasattr(req, "data")
    else:
        return req.has_data()


def request_get_data(req):
    if sys.version_info >= (3, 0, 0):
        return req.data
    else:
        return req.get_data()


__all__ = [
    "StreamingHTTPConnection",
    "StreamingHTTPHandler",
    "register_openers",
]

if hasattr(http_client, "HTTPS"):
    __all__.extend(["StreamingHTTPSHandler", "StreamingHTTPSConnection"])


class _StreamingHTTPMixin:
    """Mixin class for HTTP and HTTPS connections that implements a streaming
    send method."""

    def send(self, value):
        """Send ``value`` to the server.

        ``value`` can be a string object, a file-like object that supports
        a .read() method, or an iterable object that supports a .next()
        method.
        """
        # Based on python 2.6's httplib.HTTPConnection.send()
        if self.sock is None:
            if self.auto_open:
                self.connect()
            else:
                raise http_client.NotConnected()

        # send the data to the server. if we get a broken pipe, then close
        # the socket. we want to reconnect when somebody tries to send again.
        #
        # NOTE: we DO propagate the error, though, because we cannot simply
        #       ignore the error... the caller will know if they can retry.
        if self.debuglevel > 0:
            print("send:", repr(value))
        try:
            blocksize = 8192
            if hasattr(value, "read"):
                if hasattr(value, "seek"):
                    value.seek(0)
                if self.debuglevel > 0:
                    print("sendIng a read()able")
                data = value.read(blocksize)
                while data:
                    self.sock.sendall(data)
                    data = value.read(blocksize)
            elif hasattr(value, "next"):
                if hasattr(value, "reset"):
                    value.reset()
                if self.debuglevel > 0:
                    print("sendIng an iterable")
                for data in value:
                    self.sock.sendall(data)
            else:
                self.sock.sendall(value)
        except socket.error as v:
            if v[0] == 32:  # Broken pipe
                self.close()
            raise


class StreamingHTTPConnection(_StreamingHTTPMixin, http_client.HTTPConnection):
    """Subclass of `httplib.HTTPConnection` that overrides the `send()` method
    to support iterable body objects"""

    pass


class StreamingHTTPHandler(urllib.request.HTTPHandler):
    """Subclass of `urllib2.HTTPHandler` that uses
    StreamingHTTPConnection as its http connection class."""

    handler_order = urllib.request.HTTPHandler.handler_order - 1

    def http_open(self, req):
        """Open a StreamingHTTPConnection for the given request"""
        return self.do_open(StreamingHTTPConnection, req)

    def http_request(self, req):
        """Handle a HTTP request.  Make sure that Content-Length is specified
        if we're using an iterable value"""
        # Make sure that if we're using an iterable object as the request
        # body, that we've also specified Content-Length
        if request_has_data(req):
            data = request_get_data(req)
            if hasattr(data, "read") or hasattr(data, "next"):
                if not req.has_header("Content-length"):
                    raise ValueError("No Content-Length specified for iterable body")
        return urllib.request.HTTPHandler.do_request_(self, req)


if hasattr(http_client, "HTTPS"):

    class StreamingHTTPSConnection(_StreamingHTTPMixin, http_client.HTTPSConnection):
        """Subclass of `httplib.HTTSConnection` that overrides the `send()`
        method to support iterable body objects"""

    class StreamingHTTPSHandler(urllib.request.HTTPSHandler):
        """Subclass of `urllib2.HTTPSHandler` that uses
        StreamingHTTPSConnection as its http connection class."""

        handler_order = urllib.request.HTTPSHandler.handler_order - 1

        def https_open(self, req):
            return self.do_open(StreamingHTTPSConnection, req)

        def https_request(self, req):
            # Make sure that if we're using an iterable object as the request
            # body, that we've also specified Content-Length
            if request_has_data(req):
                data = request_get_data(req)
                if hasattr(data, "read") or hasattr(data, "next"):
                    if not req.has_header("Content-length"):
                        raise ValueError(
                            "No Content-Length specified for iterable body"
                        )
            return urllib.request.HTTPSHandler.do_request_(self, req)


def get_handlers():
    handlers = [StreamingHTTPHandler]
    if hasattr(http_client, "HTTPS"):
        handlers.append(StreamingHTTPSHandler)
    return handlers


def register_openers():
    """Register the streaming http handlers in the global urllib2 default
    opener object.

    Returns the created OpenerDirector object."""
    opener = urllib.request.build_opener(*get_handlers())

    urllib.request.install_opener(opener)

    return opener

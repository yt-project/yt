.. _debug-drive:

Debugging yt
============

There are several different convenience functions that allow you to control yt
in perhaps unexpected and unorthodox manners.  These will allow you to conduct
in-depth debugging of processes that may be running in parallel on multiple
processors, as well as providing a mechanism of signalling to yt that you need
more information about a running process.  Additionally, yt has a built-in
mechanism for optional reporting of errors to a central server.  All of these
allow for more rapid development and debugging of any problems you might
encounter.

Additionally, yt is able to leverage existing developments in the IPython
community for parallel, interactive analysis.  This allows you to initialize
multiple yt processes through ``mpirun`` and interact with all of them from a
single, unified interactive prompt.  This enables and facilitates parallel
analysis without sacrificing interactivity and flexibility.

.. _pastebin:

Pastebin
--------

A pastebin is a website where you can easily copy source code and error
messages to share with yt developers or your collaborators. At
http://paste.yt-project.org/ a pastebin is available for placing scripts.  With
yt the script ``yt_lodgeit.py`` is distributed and wrapped with
the ``pastebin`` and ``pastebin_grab`` commands, which allow for commandline
uploading and downloading of pasted snippets.  To upload a script you
would supply it to the command:

.. code-block:: bash

   $ yt pastebin some_script.py

The URL will be returned.  If you'd like it to be marked 'private' and not show
up in the list of pasted snippets, supply the argument ``--private``.  All
snippets are given either numbers or hashes.  To download a pasted snippet, you
would use the ``pastebin_grab`` option:

.. code-block:: bash

   $ yt pastebin_grab 1768

The snippet will be output to the window, so output redirection can be used to
store it in a file.

Use the Python Debugger
-----------------------

yt is almost entirely composed of python code, so it makes sense to use
the `python debugger`_ as your first stop in trying to debug it.

.. _python debugger: https://docs.python.org/3/library/pdb.html

Signaling yt to Do Something
----------------------------

During startup, yt inserts handlers for two operating system-level signals.
These provide two diagnostic methods for interacting with a running process.
Signalling the python process that is running your script with these signals
will induce the requested behavior.

   SIGUSR1
     This will cause the python code to print a stack trace, showing exactly
     where in the function stack it is currently executing.
   SIGUSR1
     This will cause the python code to insert an IPython session wherever it
     currently is, with all local variables in the local namespace.  It should
     allow you to change the state variables.

If your yt-running process has PID 5829, you can signal it to print a
traceback with:

.. code-block:: bash

   $ kill -SIGUSR1 5829

Note, however, that if the code is currently inside a C function, the signal
will not be handled, and the stacktrace will not be printed, until it returns
from that function.

.. _remote-debugging:

Remote and Disconnected Debugging
---------------------------------

If you are running a parallel job that fails, often it can be difficult to do a
post-mortem analysis to determine what went wrong.  To facilitate this, yt
has implemented an `XML-RPC <https://en.wikipedia.org/wiki/XML-RPC>`_ interface
to the Python debugger (``pdb``) event loop.

Running with the ``--rpdb`` command will cause any uncaught exception during
execution to spawn this interface, which will sit and wait for commands,
exposing the full Python debugger.  Additionally, a frontend to this is
provided through the yt command.  So if you run the command:

.. code-block:: bash

   $ mpirun -np 4 python some_script.py --parallel --rpdb

and it reaches an error or an exception, it will launch the debugger.
Additionally, instructions will be printed for connecting to the debugger.
Each of the four processes will be accessible via:

.. code-block:: bash

   $ yt rpdb 0

where ``0`` here indicates the process 0.

For security reasons, this will only work on local processes; to connect on a
cluster, you will have to execute the command ``yt rpdb`` on the node on which
that process was launched.

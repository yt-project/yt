Making Movies
=============

The process of making movies is two-step; the first is to make the frames,
which one can do in yt, but the second requires an external software package
such as QTSuperImager or mencoder, where the frames are concatenated into a
single movie file.

Making Zoomin Movies
--------------------

There is a command-line script that comes with yt, ``yt_zoomin.py``, that will
create a set of frames given command line arguments.  However, if you need more
control over the process, the following snippet of code should be a good
starting place.

.. literalinclude:: ../../../examples/cookbook_make_zoomin.py
   :language: python
   :linenos:

Making Timeseries Movies
------------------------

There is a command-line script that comes with yt, ``yt_timeseries.py``, that will
create a set of frames given command line arguments.  However, the idiom for a
very simple set of timeseries frames is given here.

.. literalinclude:: ../../../examples/cookbook_make_timeseries.py
   :language: python
   :linenos:


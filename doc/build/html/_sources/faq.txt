FAQ
===

Why Python?
-----------

Well, the easiest answer is that I knew it, and it met the requirements.  The
more interesting answer, though, is a combination of things.  Python right now
has a lot of momentum behind it, particularly in the scientific community.
It's easy to compile, portable across many architectures, relatively simple to
write C-based exntesions for, and it's well-suited to rapid application
development.  With access to an interpreter, new avenues of data-exploration
are opened, and this can lead to much more rapid and interesting analysis.

Where can I learn more about Python?
------------------------------------

There are several good, free books about Python available on the web.  The best
place to start is with the `official tutorial <http://docs.python.org/tut/>`_,
but there's also `Dive Into Python <http://www.diveintopython.org/>`_, an
entire `collection of videos <http://showmedo.com/videos/python?topic=beginner_programming>`_ on 
`ShowMeDo.com <http://showmedo.com/>`_, and a much more specific guide to using
`NumPy <http://www.scipy.org/Tentative_NumPy_Tutorial>`_, which is the backend
on which all the math done in yt is based.

As far as books go, the only book I've found to be absolutely indispensible is
the `Beazley Book <http://www.amazon.com/exec/obidos/ASIN/0735710910>`_.

Who works on yt?
----------------

Matthew Turk is the lead developer, but Jeff Oishi, Britton Smith, Dave Collins
and Stephen Skory have all made substantive contributions.

What's up with the names?
-------------------------

In the book `Snow Crash <http://en.wikipedia.org/wiki/Snow_Crash>`_, yt is
Uncle Enzo's messenger.  Lagos is the keeper of the data, Raven is a master
slicer, and so on.

Are there any restrictions on my use of yt?
-------------------------------------------

yt has been released under Version 3 of the 
`GNU General Public License <http://www.gnu.org/licenses/gpl.html>`_.

If you found it useful, and have extended it in some meaningful way, of course
I'd love to see your contributions so they can be shared with the community.
Additionally, if you use yt in a paper, I'd love it if you'd drop me a line to
let me know.

How do I know what the units returned are?
------------------------------------------

This is a very important question.  The derived fields -- and the native data
types -- are returned as CGS units, to the best knowledge of the code; but if
you see something that looks way off, you should investigate.  To see,
specifically, what yt is returning for a given field, you can do: ::

   print lagos.fieldInfo[some_field].units

and it will show you the units that have been assigned to it.

If you are defining your own derived field, you should assume that the units
given to the function you define are already in CGS.

That being said, if for some reason yt is unable to determine the correct units
for your simulation, it will notify you.  It knows how to parse output from
all of the versions of Enzo I have used or encountered, and the newest public
release is a target platform.

What are all these .yt files?
-----------------------------

By default, yt attempts to serialize a couple pieces of data that help speed it
up in future invocations.  Specifically, the entire contents of the hierarchy,
the parent-child relationships between the grids, and any projections of the
entire volume that are made.

The numbers that make up the filenames are taken from one of two places.  If
your Enzo outputs the "CurrentTimeIdentifier" parameter, that is the number that it
uses.  (This is a string of digits representing the number of seconds since the
epoch.)  If that is unavailable, the creation time of the file as reported by
the file system is used.

This can cause problems.  Ticket #91 on the bug tracker records some of these.
If you are creating multiple parameter files per second, this can lead to
incorrect hierarchies, and thus incorrect behavior, inside your analysis.  In
such cases, it is recommended that the parameter 'serialize' in the section
'lagos' of your configuration file is set to 'False'.

How can I help?
---------------

If you find a bug, report it.  If you do something cool, write it up.  If you
find a place to improve the code, send in a patch. 

Something has gone wrong.  What do I do?
----------------------------------------

Well, first off, double check that you're giving the code what it needs, and
not asking it for something it can't provide.  Use the ``help()`` command on
an object or a method to get more information.

If you can't figure out what's up, please go ahead and copy the resultant
traceback information (the error message it prints out) along with any log
files, and either send an email to the ``yt-users`` mailing list (subscribe
first!) or attach them to a ticket at `<http://yt.enzotools.org/>`_.  

.. _axis-specification:

How do I specify an axis?
-------------------------

For now, axes are specified by integers -- 0,1,2 for x,y,z.  In the next
version, this will probably change, and allow for string-identification as
well.

Where can I go for support?
---------------------------

I've set up a ``yt-users`` mailing list.  There's more information about it
at the `yt homepage <http://yt.enzotools.org>`_.

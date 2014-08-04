.. _getting-involved:

Getting Involved
================

There are *lots* of ways to get involved with yt, as a community and as a
technical system -- not all of them just contributing code, but also
participating in the community, helping us with designing the websites, adding
documentation, and sharing your scripts with others.

Coding is only one way to be involved!

Communication Channels
----------------------

There are four main communication channels for yt:

 * We have an IRC channel, on ``irc.freenode.net`` in ``#yt``.
   You can connect through our web
   gateway without any special client, at http://yt-project.org/irc.html .
   *IRC is the first stop for conversation!*
 * `yt-users <http://lists.spacepope.org/listinfo.cgi/yt-users-spacepope.org>`_
   is a relatively high-traffic mailing list where people are encouraged to ask
   questions about the code, figure things out and so on.
 * `yt-dev <http://lists.spacepope.org/listinfo.cgi/yt-dev-spacepope.org>`_ is
   a much lower-traffic mailing list designed to focus on discussions of
   improvements to the code, ideas about planning, development issues, and so
   on.
 * `yt-svn <http://lists.spacepope.org/listinfo.cgi/yt-svn-spacepope.org>`_ is
   the (now-inaccurately titled) mailing list where all pushes to the primary
   repository are sent.

The easiest way to get involved with yt is to read the mailing lists, hang out
in IRC, and participate.  If someone asks a question you know the answer to (or
have your own question about!) write back and answer it.

If you have an idea about something, suggest it!  We not only welcome
participation, we encourage it.

.. _share-your-scripts:

Share Your Scripts
------------------

.. warning:: The yt Hub is currently offline due to some hosting problems.  We
             hope to have it back up online soon.

The next easiest way to get involved with yt is to participate in the `yt Hub
<http://hub.yt-project.org/>`_.  This is a place where scripts, paper
repositories, documents and so on can be submitted to share with the broader
community.

If you have a repository on `BitBucket <https://bitbucket.org/>`_ then you can
simply submit it through the ytHub submit link.   Otherwise, we provide the
``yt hubsubmit`` command, which will guide you through the process of creating
a mercurial repository, uploading it to BitBucket, and then submitting it
directly to the Hub.

This is one of the best ways to get involved in the community!  We would love
to have more examples that show complex or advanced behavior -- and if you have
used such scripts to write a paper, that too would be an amazing contribution.

Documentation 
-------------

The yt documentation -- which you are reading right now -- is constantly being
updated, and it is a task we would very much appreciate assistance with.
Whether that is adding a section, updating an outdated section, contributing
typo or grammatical fixes, adding a FAQ, or increasing coverage of
functionality, it would be very helpful if you wanted to help out.

The easiest way to help out is to fork the main yt repository (where the
documentation lives in the ``doc`` directory in the root of the yt mercurial
repository) and then make your changes in your own fork.  When you are done,
issue a pull request through the website for your new fork, and we can comment
back and forth and eventually accept your changes.

Gallery Images and Videos
-------------------------

If you have an image or video you'd like to display in the image or video
galleries, getting it included it easy!  You can either fork the `yt homepage
repository <http://bitbucket.org/yt_analysis/website>`_ and add it there, or
email it to us and we'll add it to the `Gallery
<http://yt-project.org/gallery.html>`_.

We're eager to show off the images and movies you make with yt, so please feel 
free to drop `us <http://lists.spacepope.org/listinfo.cgi/yt-dev-spacepope.org>`_ 
a line and let us know if you've got something great!

Technical Contributions
-----------------------

Contributing code is another excellent way to participate -- whether it's
bug fixes, new features, analysis modules, or a new code frontend.  See 
:ref:`creating_frontend` for more details.

The process is pretty simple: fork on BitBucket, make changes, issue a pull
request.  We can then go back and forth with comments in the pull request, but
usually we end up accepting.

For more information, see :ref:`contributing-code`, where we spell out how to
get up and running with a development environment, how to commit, and how to
use BitBucket.

Online Presence
---------------

Some of these fall under the other items, but if you'd like to help out with
the website or any of the other ways yt is presented online, please feel free!
Almost everything is kept in hg repositories on BitBucket, and it is very easy
to fork and contribute back changes.

Please feel free to dig in and contribute changes.

Word of Mouth
-------------

If you're using yt and it has increased your productivity, please feel
encouraged to share that information.  Cite our `paper
<http://adsabs.harvard.edu/abs/2011ApJS..192....9T>`_, tell your colleagues,
and just spread word of mouth.  By telling people about your successes, you'll
help bring more eyes and hands to the table -- in this manner, by increasing
participation, collaboration, and simply spreading the limits of what the code
is asked to do, we hope to help scale the utility and capability of yt with the
community size.

Feel free to `blog <http://blog.yt-project.org/>`_ about, `tweet
<http://twitter.com/yt_astro>`_ about and talk about what you are up to!

Long-Term Projects
------------------

There are some wild-eyed, out-there ideas that have been bandied about for the
future directions of yt -- some of them even written into the mission
statement.  The ultimate goal is to move past simple analysis and visualization
of data and begin to approach it from the other side, of generating data,
running solvers.  We also hope to increase its ability to act as an in situ
analysis code, by presenting a unified protocol.  Other projects include
interfacing with ParaView and VisIt, creating a web GUI for running
simulations, creating a run-tracker that follows simulations in progress, a
federated database for simulation outputs, and so on and so forth.

yt is an ambitious project.  Let's be ambitious together.

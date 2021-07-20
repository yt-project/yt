.. _asking-for-help:

What to do if you run into problems
===================================

If you run into problems with yt, there are a number of steps to follow
to come to a solution.  The first handful of options are things you can do
on your own, but if those don't yield results, we have provided a number of
ways to connect with our community of users and developers to solve the
problem together.

To summarize, here are the steps in order:

.. contents::
   :depth: 1
   :local:
   :backlinks: none

.. _dont-panic:

Don't panic and don't give up
-----------------------------

This may seem silly, but it's effective.  While yt is a robust code with
lots of functionality, like all actively-developed codes sometimes there
are bugs.  Chances are good that your problems have a quick fix, either
because someone encountered it before and fixed it, the documentation is
out of date, or some other simple solution.  Don't give up!  We want
to help you succeed!

.. _update-the-code:

Update to the latest version
----------------------------

Sometimes the pace of development is pretty fast on yt, particularly in the
development branch, so a fix to your problem may have already been developed by
the time you encounter it.  Many users' problems can simply be corrected by
updating to the latest version of the code and/or its dependencies. If you have
installed the latest stable release of yt then you should update yt using the
package manager appropriate for your python installation. See :ref:updating.

.. _search-the-documentation:

Search the documentation, FAQ, and mailing lists
------------------------------------------------

The documentation has a lot of the answers to everyday problems.  This doesn't
mean you have to read all of the docs top-to-bottom, but you should at least
run a search to see if relevant topics have been answered in the docs.  Click
on the search field to the right of this window and enter your text.  Another
good place to look for answers in the documentation is our :ref:`faq` page.

OK, so there was no obvious solution to your problem in the documentation.
It is possible that someone else experienced the problem before you did, and
wrote to the mailing list about it.  You can easily check the mailing list
archive with the other search field to the right of this window (or you can
use the search field below).

.. raw:: html

   <form action="http://www.google.com/cse" id="cse-search-box">
     <div>
       <input type="hidden" name="cx" value="010428198273461986377:xyfd9ztykqm" />
       <input type="hidden" name="ie" value="UTF-8" />
       <input type="text" name="q" size="31" />
       <input type="submit" name="sa" value="Search" />
     </div>
   </form>
   <script type="text/javascript" src="http://www.google.com/cse/brand?form=cse-search-box&lang=en"></script>

.. _look-at-the-source:

Look at the source code
-----------------------

We've done our best to make the source clean, and it is easily searchable from
your computer.

If you have not done so already (see :ref:`install-from-source`), clone a copy
of the yt git repository and make it the 'active' installation by doing


Once inside the yt git repository, you can then search for the class,
function, or keyword which is giving you problems with ``grep -r *``, which will
recursively search throughout the code base.  (For a much faster and cleaner
experience, we recommend ``grin`` instead of ``grep -r *``.  To install ``grin``
with python, just type ``python -m pip install grin``.)

So let's say that ``SlicePlot`` is giving you problems still, and you want to
look at the source to figure out what is going on.

.. code-block:: bash

  $ cd $YT_GIT/yt
  $ grep -r SlicePlot *         (or $ grin SlicePlot)

This will print a number of locations in the yt source tree where ``SlicePlot``
is mentioned.  You can now follow-up on this and open up the files that have
references to ``SlicePlot`` (particularly the one that defines SlicePlot) and
inspect their contents for problems or clarification.

.. _isolate_and_document:

Isolate and document your problem
---------------------------------

As you gear up to take your question to the rest of the community, try to distill
your problem down to the fewest number of steps needed to produce it in a
script.  This can help you (and us) to identify the basic problem.  Follow
these steps:

* Identify what it is that went wrong, and how you knew it went wrong.
* Put your script, errors, inputs and outputs online:

  * ``$ yt pastebin script.py`` - pastes script.py online
  * ``$ yt upload_image image.png`` - pastes image online
  * ``$ yt upload my_input.tar`` - pastes my_input.tar online

* Identify which version of the code you’re using.

  * ``$ yt version`` - provides version information, including changeset hash

It may be that through the mere process of doing this, you end up solving
the problem!

.. _irc:

Go on Slack to ask a question
-----------------------------

If you want a fast, interactive experience, you could try jumping into our Slack
to get your questions answered in a chatroom style environment.

To join our slack channel you will need to request an invite by going to
https://yt-project.org/development.html, click the "Join as @ Slack!" button, and
fill out the form. You will get an invite as soon as an administrator approves
your request.

.. _mailing-list:

Ask the mailing list
--------------------

If you still haven't yet found a solution, feel free to
write to the mailing list regarding your problems.  There are two mailing lists,
`yt-users <https://mail.python.org/archives/list/yt-users@python.org/>`_ and
`yt-dev <https://mail.python.org/archives/list/yt-dev@python.org/>`_.  The
first should be used for asking for help, suggesting features and so on, and
the latter has more chatter about the way the code is developed and discussions
of changes and feature improvements.

If you email ``yt-users`` asking for help, remember to include the information
about your problem you identified in :ref:`this step <isolate_and_document>`.

When you email the list, providing this information can help the developers
understand what you did, how it went wrong, and any potential fixes or similar
problems they have seen in the past.  Without this context, it can be very
difficult to help out!

.. _reporting-a-bug:

Submit a bug report
-------------------

If you have gone through all of the above steps, and you're still encountering
problems, then you have found a bug.  To submit a bug report, you can either
directly create one through the GitHub `web interface
<https://github.com/yt-project/yt/issues/new>`_.  Alternatively, email the
``yt-users`` mailing list and we will construct a new ticket in your stead.
Remember to include the information about your problem you identified in
:ref:`this step <isolate_and_document>`.

Special Issues
--------------

Installation Issues
^^^^^^^^^^^^^^^^^^^

If you are having installation issues and nothing from the
:ref:`installation instructions <installing-yt>` seems to work, you should
*definitely* email the ``yt-users`` email list.  You should provide information
about the host, the version of the code you are using, and the output of
``yt_install.log`` from your installation.  We are very interested in making
sure that yt installs everywhere!

Customization and Scripting Issues
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have customized yt in some way, or created your own plugins file (as
described in :ref:`plugin-file`) then it may be necessary to supply users
willing to help you (or the mailing list) with both your patches to the
source, the plugin file, and perhaps even the datafile on which you're running.

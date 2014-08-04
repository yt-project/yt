What is the yt Hub?
===================

.. warning:: The yt Hub is currently offline due to some hosting problems.  We
             hope to have it back up online soon.

The yt data hub is a mechanism by which images, data objects and projects can be
shared with other people.  For instance, one can upload projections and browse
them with an interface similar to Google Maps or upload notebooks and view them
from any web browser.

.. note:: All items posted on the hub are public!

Over time, more widgets will be added, and more datatypes will be able to be
uploaded.  If you are interested in adding more ways of sharing data, please
email the developers' list.  We would like to add support for 3D widgets such
as isocontours as well as interactive binning and rebinning of data from yt
data objects, to be displayed as phase plots and profiles.

Registering a User
------------------

Because of problems with spammers, registering a user can only be done from the
yt command line.  Once you have registered a user, you can log on to the
website and obtain an API key.

To register a user:

.. code-block:: bash

   $ yt hub_register

This will walk you through the process of registering.  You will need to supply
a name, a username, a password and an email address.  Once you have gotten that
out of the way, you can go to http://hub.yt-project.org/login and log in with
your new password.  You can then receive your API key by clicking on your
username in the upper left.

After you have gotten your API key, place it in in your ``~/.yt/config`` file:

.. code-block:: none

   [yt]
   hub_api_key = 3fd8de56c2884c13a2de4dd51a80974b

Replace ``3fd8de56c2884c13a2de4dd51a80974b`` with your API key.  At this point,
you're ready to go!

What Can Be Uploaded
--------------------

Currently, the yt hub can accept these types of data:

 * Projects and script repositories: these will be displayed with an optional
   image, a description, and a link to the source repository.
 * Projections and Slices: these will be displayed in a maps-like interface,
   for interactive panning and zooming
 * IPython notebooks: these are stored on the hub and are made available for
   download and via the IPython `nbviewer <http://nbviewer.ipython.org/>`_
   service.

How to Upload Data
------------------

Uploading data takes place inside scripts.  For the most part, it is relatively
simple to do: you construct the object you would like to share, and then you
upload it.

Uploading Projects
~~~~~~~~~~~~~~~~~~

For information on how to share a project or a set of scripts, see
:ref:`share-your-scripts`.

Uploading Projections and Slices
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Projections and slices both have a ``hub_upload`` method.  Here is an example
of uploading a projection:

.. code-block:: python

   from yt.mods import *
   ds = load("IsolatedGalaxy/galaxy0030/galaxy0030")
   proj = ds.proj(0, "density", weight="density")
   proj.hub_upload()

Here is an example of uploading a slice:

.. code-block:: python

   from yt.mods import *
   ds = load("JHK-DD0030/galaxy0030")
   sl = ds.slice(0, 0.5, fields=["density"])
   sl.hub_upload()

Uploading Notebooks
~~~~~~~~~~~~~~~~~~~

Notebooks can be uploaded from the bash command line:

.. code-block:: bash

   yt upload_notebook notebook_file.ipynb

After the notebook is finished uploading, yt will print a link to the raw
notebook as well as an nbviewer link to the same notebook.  Your notebooks will
be stored under your hub profile.

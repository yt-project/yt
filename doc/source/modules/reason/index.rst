.. _reason:

yt.reason - wxPython GUI
========================

When the GUI is onscreen, you can open up a shell to not only interact with the
data already in existence but to add new data objects.   The one instance of
ReasonMainWindow is known as *mainwindow* in the namespace of the interpreter.

Additionally, within the namespace of the Reason interpreter, you have access
to *outputs*, which is a list of the outputs open in the main window, and
*data_objects*, which is every single derived data object generated in the GUI.

.. currentmodule:: yt.reason

.. autoclass:: yt.reason.ReasonMainWindow
   :members: add_static_output, add_data_object, add_cutting_wrapper
   

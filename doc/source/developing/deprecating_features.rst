How to deprecate a feature
--------------------------

Since the 4.0.0 release, deprecation happens on a per-release basis.
A functionality can be marked as deprecated using 
``~yt._maintenance.deprecation.issue_deprecation_warning``, which takes a warning
message and two version numbers, indicating the earliest release deprecating the feature
and the one in which it will be removed completely.

The message should indicate a viable alternative to replace the deprecated feature at
the user level.
``since`` and ``removal`` are required [#]_ keyword-only arguments so as to enforce
readability of the source code.

Here's an example call.

.. code-block::python

    def old_function(*args, **kwargs):
        from yt._maintenance.deprecation import issue_deprecation_warning
        issue_deprecation_warning(
            "`old_function` is deprecated, use `replacement_function` instead."
            since="4.0.0",
            removal="4.1.0"
        )
        ...

.. [#] ``since`` is not required yet as of yt 4.0.0 because existing warnings predate its introduction. 
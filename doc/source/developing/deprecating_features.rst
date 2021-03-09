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

If a whole function or class is marked as deprecated, it should be removed from
``doc/source/reference/api/api.rst``.


.. [#] ``since`` is not required yet as of yt 4.0.0 because existing warnings predate its introduction.

Deprecating Derived Fields
--------------------------

Occasionally, one may want to deprecate a derived field in yt, normally
because naming conventions for fields have changed, or simply because a
field has outlived its usefulness. There are two ways to mark fields as
deprecated in yt.

The first way is if you simply want to mark a specific derived field as
deprecated. In that case, you call
:meth:`~yt.fields.field_info_container.FieldInfoContainer.add_deprecated_field`:

.. code-block:: python

    def _cylindrical_radial_absolute(field, data):
        """This field is deprecated and will be removed in a future version"""
        return np.abs(data[ftype, f"{basename}_cylindrical_radius"])


    registry.add_deprecated_field(
        (ftype, f"cylindrical_radial_{basename}_absolute"),
        sampling_type="local",
        function=_cylindrical_radial_absolute,
        since="4.0.0",
        removal="4.1.0",
        units=field_units,
        validators=[ValidateParameter("normal")],
    )

Note that the signature for
:meth:`~yt.fields.field_info_container.FieldInfoContainer.add_deprecated_field`
is the same as :meth:`~yt.fields.field_info_container.FieldInfoContainer.add_field`,
with the exception of the ``since`` and ``removal`` arguments which indicate in
what version the field was deprecated and in what version it will be removed.
The effect is to add a warning to the logger when the field is first used:

.. code-block:: python

    import yt

    ds = yt.load("GasSloshing/sloshing_nomag2_hdf5_plt_cnt_0100")
    sp = ds.sphere("c", (100.0, "kpc"))
    print(sp["gas", "cylindrical_radial_velocity_absolute"])

.. code-block:: pycon

    yt : [WARNING  ] 2021-03-09 16:30:47,460 The Derived Field
    ('gas', 'cylindrical_radial_velocity_absolute') is deprecated
    as of yt v4.0.0 and will be removed in yt v4.1.0

The second way to deprecate a derived field is to take an existing field
definition and change its name. In order to mark the original name as deprecated,
use the :meth:`~yt.fields.field_info_container.FieldInfoContainer.alias` method
and pass the ``since`` and ``removal`` arguments (see above) as a tuple in the
``deprecate`` keyword argument:

.. code-block:: python

    registry.alias(
        (ftype, "kinetic_energy"),
        (ftype, "kinetic_energy_density"),
        deprecate=("4.0.0", "4.1.0"),
    )

Note that the old field name which is to be deprecated goes first, and the new,
replacement field name goes second. In this case, the log message reports to
the user what field they should use:

.. code-block:: python

    print(sp["gas", "kinetic_energy"])

.. code-block:: pycon

    yt : [WARNING  ] 2021-03-09 16:29:12,911 The Derived Field
    ('gas', 'kinetic_energy') is deprecated as of yt v4.0.0 and will be removed
    in yt v4.1.0 Use ('gas', 'kinetic_energy_density') instead.

In most cases, the ``since`` and ``removal`` arguments should have a delta of
one minor release, and that should be the minimum value. However, the developer
is free to use their judgment about whether or not the delta should be multiple
minor releases if the field has a long provenance.

import numpy as np

from yt.utilities.on_demand_imports import _nibabel as nib

from .cartesian_coordinates import CartesianCoordinateHandler


class CartesianSubjectCoordinateHandler(CartesianCoordinateHandler):
    name = "cartesian_subject"

    def setup_fields(self, registry):
        super().setup_fields(registry)

        def _affine_transformed_x(field, data):
            x = data["index", "x"].in_units("code_length")
            y = data["index", "y"].in_units("code_length")
            z = data["index", "z"].in_units("code_length")
            pos = np.stack([x, y, z], axis=-1)
            new_x = (nib.apply_affine(data.ds.parameters["affine"], pos) * x.uq)[..., 0]
            return new_x

        registry.add_field(
            ("index", "affine_transformed_x"),
            sampling_type="cell",
            function=_affine_transformed_x,
            display_field=False,
            units="code_length",
        )

        def _affine_transformed_y(field, data):
            x = data["index", "x"].in_units("code_length")
            y = data["index", "y"].in_units("code_length")
            z = data["index", "z"].in_units("code_length")
            pos = np.stack([x, y, z], axis=-1)
            new_y = (nib.apply_affine(data.ds.parameters["affine"], pos) * x.uq)[..., 1]
            return new_y

        registry.add_field(
            ("index", "affine_transformed_y"),
            sampling_type="cell",
            function=_affine_transformed_y,
            display_field=False,
            units="code_length",
        )

        def _affine_transformed_z(field, data):
            x = data["index", "x"].in_units("code_length")
            y = data["index", "y"].in_units("code_length")
            z = data["index", "z"].in_units("code_length")
            pos = np.stack([x, y, z], axis=-1)
            new_z = (nib.apply_affine(data.ds.parameters["affine"], pos) * x.uq)[..., 2]
            return new_z

        registry.add_field(
            ("index", "affine_transformed_z"),
            sampling_type="cell",
            function=_affine_transformed_z,
            display_field=False,
            units="code_length",
        )

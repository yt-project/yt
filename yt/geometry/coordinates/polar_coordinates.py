from .cylindrical_coordinates import CylindricalCoordinateHandler


class PolarCoordinateHandler(CylindricalCoordinateHandler):
    name = "polar"

    def __init__(self, ds, ordering=("r", "theta", "z")):
        super().__init__(ds, ordering)
        # No need to set labels here

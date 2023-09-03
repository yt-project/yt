# These are some of the data object registries that are used in different places in the
# code. Not all of the self-registering objects are included in these.

from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from typing import Union

    from yt.data_objects.analyzer_objects import AnalysisTask
    from yt.data_objects.data_containers import YTDataContainer
    from yt.data_objects.derived_quantities import DerivedQuantity
    from yt.data_objects.static_output import Dataset
    from yt.data_objects.time_series import DatasetSeries
    from yt.visualization.volume_rendering.old_camera import Camera

analysis_task_registry: dict[str, type["AnalysisTask"]] = {}
derived_quantity_registry: dict[str, type["DerivedQuantity"]] = {}
output_type_registry: dict[str, type["Dataset"]] = {}
simulation_time_series_registry: dict[str, type["DatasetSeries"]] = {}

# TODO: split into 2 registries to avoid a typing.Union
data_object_registry: dict[str, "Union[type[YTDataContainer], type[Camera]]"] = {}

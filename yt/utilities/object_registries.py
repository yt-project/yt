# These are some of the data object registries that are used in different places in the
# code. Not all of the self-registering objects are included in these.

# type hints are simplified as raw Type (instead of, e.g., Type[Dataset])
# to workaround circular imports

# subclasses of yt.data_objects.analyzer_objects.AnalysisTask
analysis_task_registry: dict[str, type] = {}

# subclasses of yt.data_objects.data_containers.YTDataContainer
data_object_registry: dict[str, type] = {}

# suclasses of yt.data_objects.derived_quantity.DerivedQuantity
derived_quantity_registry: dict[str, type] = {}

# suclasses of yt.data_objects.static_outputs.Dataset
output_type_registry: dict[str, type] = {}

# subclasses of yt.data_objects.time_series.DatasetSeries
simulation_time_series_registry: dict[str, type] = {}

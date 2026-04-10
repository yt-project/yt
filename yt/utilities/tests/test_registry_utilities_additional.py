import copy

from numpy.testing import assert_allclose

from yt.utilities import definitions
from yt.utilities.object_registries import (
    analysis_task_registry,
    data_object_registry,
    derived_quantity_registry,
    output_type_registry,
    simulation_time_series_registry,
)
from yt.utilities.operator_registry import OperatorRegistry
from yt.utilities.tree_container import TreeContainer


class _DummyOp:
    def __init__(self, value):
        self.value = value
        self.args = ()
        self.kwargs = {}


class _TreeNode(TreeContainer):
    def __init__(self, name, children=None):
        super().__init__()
        self.name = name
        self.children = children


def test_definitions_match_conversion_constants():
    assert definitions.MAXLEVEL == 48
    assert_allclose(definitions.mpc_conversion["km"], definitions.km_per_mpc)
    assert_allclose(definitions.mpc_conversion["cm"], definitions.cm_per_mpc)
    assert definitions.formatted_length_unit_names["au"] == "AU"
    assert definitions.formatted_length_unit_names["rsun"] == r"R_\odot"
    assert definitions.formatted_length_unit_names["code_length"] == r"code\ length"
    assert_allclose(definitions.sec_conversion["days"], definitions.sec_per_day)
    assert_allclose(
        definitions.sec_conversion["Gyr"] / definitions.sec_conversion["Myr"],
        1.0e3,
    )


def test_tree_container_iterates_depth_first():
    leaf_a = _TreeNode("leaf_a")
    leaf_b = _TreeNode("leaf_b")
    root = _TreeNode("root", children=[_TreeNode("branch", [leaf_a]), leaf_b])

    assert [node.name for node in root] == ["root", "branch", "leaf_a", "leaf_b"]
    assert [node.name for node in leaf_a] == ["leaf_a"]


def test_operator_registry_find_deepcopies_registered_operators():
    registry = OperatorRegistry()
    registry["sum"] = _DummyOp(value=4.0)

    found = registry.find("sum", 1.0, 2.0, scale=3.0)

    assert found is not registry["sum"]
    assert found.value == 4.0
    assert found.args == (1.0, 2.0)
    assert found.kwargs == {"scale": 3.0}

    direct = _DummyOp(value=7.0)
    assert registry.find(direct) is direct


def test_object_registries_accept_and_retain_class_entries():
    registries = (
        analysis_task_registry,
        derived_quantity_registry,
        output_type_registry,
        simulation_time_series_registry,
        data_object_registry,
    )
    originals = [copy.copy(registry) for registry in registries]

    class AnalysisTask:
        pass

    class DerivedQuantity:
        pass

    class Dataset:
        pass

    class DatasetSeries:
        pass

    class DataObject:
        pass

    try:
        analysis_task_registry["analysis_task"] = AnalysisTask
        derived_quantity_registry["derived_quantity"] = DerivedQuantity
        output_type_registry["dataset"] = Dataset
        simulation_time_series_registry["series"] = DatasetSeries
        data_object_registry["data_object"] = DataObject

        assert analysis_task_registry["analysis_task"] is AnalysisTask
        assert derived_quantity_registry["derived_quantity"] is DerivedQuantity
        assert output_type_registry["dataset"] is Dataset
        assert simulation_time_series_registry["series"] is DatasetSeries
        assert data_object_registry["data_object"] is DataObject
    finally:
        for registry, original in zip(registries, originals, strict=True):
            registry.clear()
            registry.update(original)

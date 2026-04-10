from yt.utilities.hierarchy_inspection import get_classes_with_missing_requirements
from yt.utilities.object_registries import output_type_registry


class MissingBeta:
    @staticmethod
    def _missing_load_requirements():
        return ["beta_backend"]


class MissingAlpha:
    @staticmethod
    def _missing_load_requirements():
        return ["alpha_backend", "alpha_table"]


class FullyAvailable:
    @staticmethod
    def _missing_load_requirements():
        return []


def test_get_classes_with_missing_requirements_filters_and_sorts():
    original_registry = output_type_registry.copy()
    try:
        output_type_registry.clear()
        output_type_registry.update(
            {
                "beta": MissingBeta,
                "available": FullyAvailable,
                "alpha": MissingAlpha,
            }
        )

        missing = get_classes_with_missing_requirements()

        assert list(missing) == [MissingAlpha, MissingBeta]
        assert missing[MissingAlpha] == ["alpha_backend", "alpha_table"]
        assert missing[MissingBeta] == ["beta_backend"]
    finally:
        output_type_registry.clear()
        output_type_registry.update(original_registry)

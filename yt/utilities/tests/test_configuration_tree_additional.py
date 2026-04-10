import pytest

from yt.utilities.configuration_tree import ConfigLeaf, ConfigNode


def test_config_node_round_trip_update_and_display_helpers():
    root = ConfigNode(None)
    manual = ConfigNode("manual")
    root.add("manual", manual)
    assert manual.parent is root
    root.remove_child("manual")

    root.add_child("empty")
    assert "empty" in root
    assert repr(root) == "<Node None>"

    metadata = {"source": "unit-test.toml"}
    root.update(
        {"yt": {"internals": {"parallel": True}, "answer": 42}},
        extra_data=metadata,
    )

    assert root.get("yt", "answer").value == 42
    assert root.get_leaf("yt", "answer") == 42
    assert (
        root.get_leaf("yt", "answer", callback=lambda leaf: leaf.extra_data["source"])
        == "unit-test.toml"
    )
    assert root.get_deepest_leaf("yt", "internals", "missing_branch", "answer") == 42
    assert root.serialize() == {
        "empty": {},
        "yt": {"internals": {"parallel": True}, "answer": 42},
    }
    assert root.as_dict() == {"yt": {"internals": {"parallel": True}, "answer": 42}}
    assert root._repr_json_() == root.as_dict()
    assert root.as_dict(
        callback=lambda leaf: {"value": leaf.value, "source": leaf.extra_data["source"]}
    ) == {
        "yt": {
            "internals": {"parallel": {"value": True, "source": "unit-test.toml"}},
            "answer": {"value": 42, "source": "unit-test.toml"},
        }
    }

    root.pop_leaf(["yt", "answer"])
    assert root.as_dict() == {"yt": {"internals": {"parallel": True}}}

    rebuilt = ConfigNode.from_dict(
        {"section": {"subsection": 3}, "name": "yt"},
        extra_data={"source": "from-dict"},
    )
    assert rebuilt.serialize() == {"section": {"subsection": 3}, "name": "yt"}
    assert rebuilt.get_leaf("section", "subsection") == 3
    assert rebuilt.get("name").extra_data == {"source": "from-dict"}

    ordered = ConfigNode(None)
    ordered_branch = ConfigNode("branch")
    ordered_branch.add(
        "value",
        ConfigLeaf("value", parent=ordered_branch, value=7, extra_data={"source": "ordered"}),
    )
    ordered.add("branch", ordered_branch)
    ordered.add(
        "tail",
        ConfigLeaf("tail", parent=ordered, value=9, extra_data={"source": "ordered"}),
    )
    assert ordered.as_dict() == {"branch": {"value": 7}, "tail": 9}


def test_config_node_error_paths_cover_leaf_node_mismatches():
    root = ConfigNode(None)
    root.add("branch", ConfigNode("branch"))
    root.add("leaf", ConfigLeaf("leaf", parent=root, value=2, extra_data={}))

    with pytest.raises(KeyError, match="Cannot get key missing"):
        root.get_child("missing")

    with pytest.raises(RuntimeError, match="Expected a ConfigLeaf"):
        root.upsert_from_list(["branch"], 1, extra_data={})

    with pytest.raises(RuntimeError, match="Expected a ConfigNode"):
        root.upsert_from_list(["leaf", "child"], 3, extra_data={})

    search_root = ConfigNode(None)
    search_root.add("yt", ConfigNode("yt"))
    search_root.get("yt").add("internals", ConfigNode("internals"))

    with pytest.raises(KeyError, match="Cannot any node that contains the leaf answer"):
        search_root.get_deepest_leaf("yt", "internals", "answer")

    search_root.get("yt", "internals").add(
        "answer",
        ConfigNode("answer", parent=search_root.get("yt", "internals")),
    )
    with pytest.raises(RuntimeError, match="Expected a ConfigLeaf"):
        search_root.get_deepest_leaf("yt", "internals", "answer")


def test_config_leaf_helpers_and_type_validation_messages():
    root = ConfigNode(None)
    section = ConfigNode("yt")
    root.add("yt", section)

    leaf = ConfigLeaf("answer", parent=section, value=42, extra_data={"source": "cfg.toml"})
    section.add("answer", leaf)
    assert [node.key for node in leaf.get_tree()] == [None, "yt", "answer"]
    assert leaf.serialize() == 42
    assert repr(leaf) == "<Leaf answer: 42>"

    leaf.value = 43
    assert leaf.value == 43

    with pytest.raises(TypeError, match="Error when setting yt.answer") as exc_info:
        leaf.value = "forty-three"
    message = str(exc_info.value)
    assert "expected type <class 'int'>" in message
    assert "This entry was last modified in cfg.toml." in message

    leaf_without_source = ConfigLeaf("parallel", parent=section, value=True, extra_data={})
    with pytest.raises(TypeError, match="Error when setting yt.parallel") as exc_info:
        leaf_without_source.value = 1
    assert "last modified" not in str(exc_info.value)

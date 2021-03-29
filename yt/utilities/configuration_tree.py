class ConfigNode:
    def __init__(self, key, parent=None):
        self.key = key
        self.children = {}
        self.parent = parent

    def add(self, key, child):
        self.children[key] = child
        child.parent = self

    def update(self, other, extra_data=None):
        def _recursive_upsert(other_dict, keys):
            for key, val in other_dict.items():
                new_keys = keys + [key]
                if isinstance(val, dict):
                    _recursive_upsert(val, new_keys)
                else:
                    self.upsert_from_list(new_keys, val, extra_data)

        _recursive_upsert(other, keys=[])

    def get_child(self, key, constructor=None):
        if key in self.children:
            child = self.children[key]
        elif constructor is not None:
            child = self.children[key] = constructor()
        else:
            raise KeyError(f"Cannot get key {key}")
        return child

    def add_child(self, key):
        self.get_child(key, lambda: ConfigNode(key, parent=self))

    def remove_child(self, key):
        self.children.pop(key)

    def upsert_from_list(self, keys, value, extra_data=None):
        key, *next_keys = keys
        if len(next_keys) == 0:  # reach the end of the upsert
            leaf = self.get_child(
                key,
                lambda: ConfigLeaf(
                    key, parent=self, value=value, extra_data=extra_data
                ),
            )
            leaf.value = value
            leaf.extra_data = extra_data
            if not isinstance(leaf, ConfigLeaf):
                raise RuntimeError(f"Expected a ConfigLeaf, got {leaf}!")
        else:
            next_node = self.get_child(key, lambda: ConfigNode(key, parent=self))
            if not isinstance(next_node, ConfigNode):
                raise RuntimeError(f"Expected a ConfigNode, got {next_node}!")
            next_node.upsert_from_list(next_keys, value, extra_data)

    def get_from_list(self, key_list):
        next, *key_list_remainder = key_list
        child = self.get_child(next)
        if len(key_list_remainder) == 0:
            return child
        else:
            return child.get_from_list(key_list_remainder)

    def get(self, *keys):
        return self.get_from_list(keys)

    def get_leaf(self, *keys, callback=lambda leaf: leaf.value):
        leaf = self.get_from_list(keys)
        return callback(leaf)

    def pop_leaf(self, keys):
        *node_keys, leaf_key = keys
        node = self.get_from_list(node_keys)
        node.children.pop(leaf_key)

    def get_deepest_leaf(self, *keys, callback=lambda leaf: leaf.value):
        root_key, *keys, leaf_key = keys

        root_node = self.get_child(root_key)
        node_list = [root_node]
        node = root_node

        # Traverse the tree down following the keys
        for k in keys:
            try:
                node = node.get_child(k)
                node_list.append(node)
            except KeyError:
                break

        # For each node, starting from the deepest, try to find the leaf
        for node in reversed(node_list):
            try:
                leaf = node.get_child(leaf_key)
                if not isinstance(leaf, ConfigLeaf):
                    raise RuntimeError(f"Expected a ConfigLeaf, got {leaf}!")
                return callback(leaf)
            except KeyError:
                continue
        raise KeyError(f"Cannot any node that contains the leaf {leaf_key}.")

    def serialize(self):
        retval = {}
        for key, child in self.children.items():
            retval[key] = child.serialize()
        return retval

    @staticmethod
    def from_dict(other, parent=None, **kwa):
        me = ConfigNode(None, parent=parent)
        for key, val in other.items():
            if isinstance(val, dict):
                me.add(key, ConfigNode.from_dict(val, parent=me, **kwa))
            else:
                me.add(key, ConfigLeaf(key, parent=me, value=val, **kwa))
        return me

    def _as_dict_with_count(self, callback):
        data = {}
        total_count = 0
        for key, child in self.children.items():
            if isinstance(child, ConfigLeaf):
                total_count += 1
                data[key] = callback(child)
            elif isinstance(child, ConfigNode):
                child_data, count = child._as_dict_with_count(callback)
                total_count += count
                if count > 0:
                    data[key] = child_data

        return data, total_count

    def as_dict(self, callback=lambda child: child.value):
        data, _ = self._as_dict_with_count(callback)
        return data

    def __repr__(self):
        return f"<Node {self.key}>"

    def __contains__(self, item):
        return item in self.children

    # Add support for IPython rich display
    # see https://ipython.readthedocs.io/en/stable/config/integrating.html
    def _repr_json_(self):
        return self.as_dict()


class ConfigLeaf:
    def __init__(self, key, parent: ConfigNode, value, extra_data=None):
        self.key = key  # the name of the config leaf
        self._value = value
        self.parent = parent
        self.extra_data = extra_data

    def serialize(self):
        return self.value

    def get_tree(self):
        node = self
        parents = []
        while node is not None:
            parents.append(node)
            node = node.parent

        return reversed(parents)

    @property
    def value(self):
        return self._value

    @value.setter
    def value(self, new_value):
        if type(self.value) == type(new_value):
            self._value = new_value
        else:
            tree = self.get_tree()
            tree_str = ".".join(node.key for node in tree if node.key)
            msg = f"Error when setting {tree_str}.\n"
            msg += (
                "Tried to assign a value of type "
                f"{type(new_value)}, expected type {type(self.value)}."
            )
            source = self.extra_data.get("source", None)
            if source:
                msg += f"\nThis entry was last modified in {source}."
            raise TypeError(msg)

    def __repr__(self):
        return f"<Leaf {self.key}: {self.value}>"

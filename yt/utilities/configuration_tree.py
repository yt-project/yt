class ConfigNode:
    def __init__(self, parent=None):
        self.children = {}
        self.parent = parent

    def add(self, key, child):
        self.children[key] = child
        child.parent = self

    def update(self, other, extraData=None):
        def _helper(d, keys):
            for key, val in d.items():
                new_keys = keys + [key]
                if isinstance(val, dict):
                    _helper(val, new_keys)
                else:
                    self.upsert_from_list(new_keys, val, extraData)

        _helper(other, [])

    def get_child(self, key, constructor=None):
        if key in self.children:
            child = self.children[key]
        elif constructor is not None:
            child = self.children[key] = constructor()
        else:
            raise KeyError(f"Cannot get key {key}")
        return child

    def add_child(self, key):
        self.get_child(key, lambda: ConfigNode(parent=self))

    def remove_child(self, key):
        self.children.pop(key)

    def upsert_from_list(self, keys, value, extraData=None):
        key, *next_keys = keys
        if len(next_keys) == 0:  # reach the end of the upsert
            leaf = self.get_child(key, lambda: ConfigLeaf(self, value, extraData))
            leaf.value = value
            leaf.extraData = extraData
            if not isinstance(leaf, ConfigLeaf):
                raise RuntimeError(f"Expected a ConfigLeaf, got {leaf}!")
        else:
            next_node = self.get_child(key, lambda: ConfigNode(parent=self))
            if not isinstance(next_node, ConfigNode):
                raise RuntimeError(f"Expected a ConfigNode, got {next_node}!")
            next_node.upsert_from_list(next_keys, value, extraData)

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

    def get_most_common(self, keys, leaf):
        root_key, *keys = keys

        root_node = self.get_child(root_key)
        node_list = [root_node]
        node = root_node
        for k in keys:
            try:
                node = node.get_child(k)
                node_list.append(node)
            except KeyError:
                break

        for node in reversed(node_list):
            try:
                return node.get_child(leaf)
            except KeyError:
                continue
        raise KeyError(f"Cannot any node that contains the leaf {leaf}.")

    def serialize(self):
        retval = {}
        for key, child in self.children.items():
            retval[key] = child.serialize()
        return retval

    def __repr__(self):
        return f"<Node: {list(self.children.keys())}>"

    @staticmethod
    def from_dict(other, parent=None, **kwa):
        me = ConfigNode(parent=parent)
        for key, val in other.items():
            if isinstance(val, dict):
                me.add(key, ConfigNode.from_dict(val, parent=me, **kwa))
            else:
                me.add(key, ConfigLeaf(me, val, **kwa))
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


class ConfigLeaf:
    def __init__(self, parent: ConfigNode, value, extraData=None):
        self.value = value
        self.extra = extraData
        self.parent = parent

    def serialize(self):
        return self.value

    def __repr__(self):
        return f"<Leaf: {self.value}>"

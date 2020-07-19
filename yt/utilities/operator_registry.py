import copy


class OperatorRegistry(dict):
    def find(self, op, *args, **kwargs):
        if isinstance(op, str):
            # Lookup, assuming string or hashable object
            op = copy.deepcopy(self[op])
            op.args = args
            op.kwargs = kwargs
        return op

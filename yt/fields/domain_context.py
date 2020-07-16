domain_context_registry = {}


class DomainContext:
    class __metaclass__(type):
        def __init__(cls, name, b, d):
            type.__init__(cls, name, b, d)
            domain_context_registry[name] = cls

    _known_fluid_fields = ()
    _known_particle_fields = ()

    def __init__(self, ds):
        self.ds = ds

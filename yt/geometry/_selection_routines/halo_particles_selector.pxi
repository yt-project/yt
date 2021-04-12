cdef class HaloParticlesSelector(SelectorObject):
    cdef public object base_source
    cdef SelectorObject base_selector
    cdef object pind
    cdef public np.int64_t halo_id
    def __init__(self, dobj):
        self.base_source = dobj.base_source
        self.base_selector = self.base_source.selector
        self.pind = dobj.particle_indices

    def _hash_vals(self):
        return ("halo_particles", self.halo_id)

halo_particles_selector = HaloParticlesSelector

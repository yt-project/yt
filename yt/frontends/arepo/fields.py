from yt.frontends.gadget.api import GadgetFieldInfo

class ArepoFieldInfo(GadgetFieldInfo):
    known_particle_fields = GadgetFieldInfo.known_particle_fields + \
                            (("smoothing_length", ("code_length", [], None)),)

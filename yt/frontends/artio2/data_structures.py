#from _artio_reader import artio_is_valid
from _artio_reader import artio_is_valid, artio_fileset

from yt.data_objects.static_output import \
    StaticOutput

class Artio2StaticOutput(StaticOutput):
    _handle = None

    def __init__(self, filename, data_style='artio2'):
        if self._handle is not None : return
        self.fileset_prefix = filename[:-4]
        self._fileset = artio_fileset(self.fileset_prefix)
        StaticOutput.__init__(self, filename, data_style)

    #def _parse_parameter_file(self) :
        #self.parameters = _artio_reader.artio_read_parameters( self._handle )

    @classmethod
    def _is_valid(self, *args, **kwargs) :
        # a valid artio header file starts with a prefix and ends with .art
        if not args[0].endswith(".art"): return False
        return artio_is_valid(args[0][:-4])

from yt.utilities.io_handler import \
    BaseIOHandler

class IOHandlerArtio2(BaseIOHandler) :
    _data_style = "artio2"

    def _read_fluid_selection(self, chunks, selector, fields, size):
        return

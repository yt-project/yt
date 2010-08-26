class TigerFieldContainer(CodeFieldInfoContainer):
    """
    This is a container for Tiger-specific fields.
    """
    _shared_state = {}
    _field_list = {}
TigerFieldInfo = TigerFieldContainer()
add_tiger_field = TigerFieldInfo.add_field


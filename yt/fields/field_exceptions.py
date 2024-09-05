class ValidationException(Exception):
    pass


class NeedsGridType(ValidationException):
    def __init__(self, ghost_zones=0, fields=None):
        self.ghost_zones = ghost_zones
        self.fields = fields

    def __str__(self):
        s = "s" if self.ghost_zones != 1 else ""
        return f"fields {self.fields} require {self.ghost_zones} ghost zone{s}."


class NeedsOriginalGrid(NeedsGridType):
    def __init__(self):
        self.ghost_zones = 0


class NeedsDataField(ValidationException):
    def __init__(self, missing_fields):
        self.missing_fields = missing_fields

    def __str__(self):
        return f"({self.missing_fields})"


class NeedsProperty(ValidationException):
    def __init__(self, missing_properties):
        self.missing_properties = missing_properties

    def __str__(self):
        return f"({self.missing_properties})"


class NeedsParameter(ValidationException):
    def __init__(self, missing_parameters):
        self.missing_parameters = missing_parameters

    def __str__(self):
        return f"({self.missing_parameters})"


class NeedsConfiguration(ValidationException):
    def __init__(self, parameter, value):
        self.parameter = parameter
        self.value = value

    def __str__(self):
        return f"(Needs {self.parameter} = {self.value})"


class FieldUnitsError(Exception):
    pass

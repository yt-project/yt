# We don't need to import 'exceptions'
import os.path
from typing import List, Tuple

from unyt.exceptions import UnitOperationError


class YTException(Exception):
    def __init__(self, message=None, ds=None):
        Exception.__init__(self, message)
        self.ds = ds


# Data access exceptions:


class YTUnidentifiedDataType(YTException):
    def __init__(self, filename, *args, **kwargs):
        self.filename = filename
        self.args = args
        self.kwargs = kwargs

    def __str__(self):
        msg = f"Could not determine input format from {self.filename!r}"
        if self.args:
            msg += ", " + (", ".join(f"{a!r}" for a in self.args))
        if self.kwargs:
            msg += ", " + (", ".join(f"{k}={v!r}" for k, v in self.kwargs.items()))
        return msg


class YTAmbiguousDataType(YTUnidentifiedDataType):
    def __init__(self, filename, candidates):
        self.filename = filename
        self.candidates = candidates

    def __str__(self):
        msg = f"Multiple data type candidates for {self.filename}\n"
        msg += "The following independent classes were detected as valid :\n"
        for c in self.candidates:
            msg += f"{c}\n"
        msg += (
            "This degeneracy can be lifted using the `hint` keyword argument in yt.load"
        )
        return msg


class YTSphereTooSmall(YTException):
    def __init__(self, ds, radius, smallest_cell):
        YTException.__init__(self, ds=ds)
        self.radius = radius
        self.smallest_cell = smallest_cell

    def __str__(self):
        return f"{self.radius:0.5e} < {self.smallest_cell:0.5e}"


class YTAxesNotOrthogonalError(YTException):
    def __init__(self, axes):
        self.axes = axes

    def __str__(self):
        return f"The supplied axes are not orthogonal.  {self.axes}"


class YTNoDataInObjectError(YTException):
    def __init__(self, obj):
        self.obj_type = getattr(obj, "_type_name", "")

    def __str__(self):
        s = "The object requested has no data included in it."
        if self.obj_type == "slice":
            s += "  It may lie on a grid face.  Try offsetting slightly."
        return s


class YTFieldNotFound(YTException):
    def __init__(self, field, ds):
        self.field = field
        self.ds = ds

    def _get_suggestions(self) -> List[Tuple[str, str]]:
        from yt.funcs import levenshtein_distance

        field = self.field
        ds = self.ds

        suggestions = {}
        if not isinstance(field, tuple):
            ftype, fname = None, field
        elif field[1] is None:
            ftype, fname = None, field[0]
        else:
            ftype, fname = field

        # Limit the suggestions to a distance of 3 (at most 3 edits)
        # This is very arbitrary, but is picked so that...
        # - small typos lead to meaningful suggestions (e.g. `densty` -> `density`)
        # - we don't suggest unrelated things (e.g. `pressure` -> `density` has a distance
        #   of 6, we definitely do not want it)
        # A threshold of 3 seems like a good middle point.
        max_distance = 3

        # Suggest (ftype, fname), with alternative ftype
        for ft, fn in ds.derived_field_list:
            if fn.lower() == fname.lower() and (
                ftype is None or ft.lower() != ftype.lower()
            ):
                suggestions[ft, fn] = 0

        if ftype is not None:
            # Suggest close matches using levenshtein distance
            fields_str = {_: str(_).lower() for _ in ds.derived_field_list}
            field_str = str(field).lower()

            for (ft, fn), fs in fields_str.items():
                distance = levenshtein_distance(field_str, fs, max_dist=max_distance)
                if distance < max_distance:
                    if (ft, fn) in suggestions:
                        continue
                    suggestions[ft, fn] = distance

        # Return suggestions sorted by increasing distance (first are most likely)
        return [
            (ft, fn)
            for (ft, fn), distance in sorted(suggestions.items(), key=lambda v: v[1])
        ]

    def __str__(self):
        msg = f"Could not find field {self.field} in {self.ds}."
        try:
            suggestions = self._get_suggestions()
        except AttributeError:
            # This may happen if passing a field that is e.g. an Ellipsis
            # e.g. when using ds.r[...]
            suggestions = []
        if suggestions:
            msg += "\nDid you mean:\n\t"
            msg += "\n\t".join(str(_) for _ in suggestions)
        return msg


class YTParticleTypeNotFound(YTException):
    def __init__(self, fname, ds):
        self.fname = fname
        self.ds = ds

    def __str__(self):
        return f"Could not find particle_type {self.fname!r} in {self.ds}."


class YTSceneFieldNotFound(YTException):
    pass


class YTCouldNotGenerateField(YTFieldNotFound):
    def __str__(self):
        return f"Could field '{self.fname}' in {self.ds} could not be generated."


class YTFieldTypeNotFound(YTException):
    def __init__(self, ftype, ds=None):
        self.ftype = ftype
        self.ds = ds

    def __str__(self):
        if self.ds is not None and self.ftype in self.ds.particle_types:
            return (
                f"Could not find field type {self.ftype!r}. "
                "This field type is a known particle type for this dataset. "
                "Try adding this field with sampling_type='particle'."
            )
        else:
            return f"Could not find field type {self.ftype!r}."


class YTSimulationNotIdentified(YTException):
    def __init__(self, sim_type):
        YTException.__init__(self)
        self.sim_type = sim_type

    def __str__(self):
        return f"Simulation time-series type {self.sim_type!r} not defined."


class YTCannotParseFieldDisplayName(YTException):
    def __init__(self, field_name, display_name, mathtext_error):
        self.field_name = field_name
        self.display_name = display_name
        self.mathtext_error = mathtext_error

    def __str__(self):
        return (
            f"The display name {self.display_name!r} of the derived field {self.field_name!r} "
            f"contains the following LaTeX parser errors:\n{self.mathtext_error}"
        )


class YTCannotParseUnitDisplayName(YTException):
    def __init__(self, field_name, unit_name, mathtext_error):
        self.field_name = field_name
        self.unit_name = unit_name
        self.mathtext_error = mathtext_error

    def __str__(self):
        return (
            f"The unit display name {self.unit_name!r} of the derived field {self.field_name!r} "
            f"contains the following LaTeX parser errors:\n{self.mathtext_error}"
        )


class InvalidSimulationTimeSeries(YTException):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message


class MissingParameter(YTException):
    def __init__(self, ds, parameter):
        YTException.__init__(self, ds=ds)
        self.parameter = parameter

    def __str__(self):
        return f"dataset {self.ds} is missing {self.parameter} parameter."


class NoStoppingCondition(YTException):
    def __init__(self, ds):
        YTException.__init__(self, ds=ds)

    def __str__(self):
        return f"Simulation {self.ds} has no stopping condition. StopTime or StopCycle should be set."


class YTNotInsideNotebook(YTException):
    def __str__(self):
        return "This function only works from within an IPython Notebook."


class YTGeometryNotSupported(YTException):
    def __init__(self, geom):
        self.geom = geom

    def __str__(self):
        return f"We don't currently support {self.geom!r} geometry"


class YTCoordinateNotImplemented(YTException):
    def __str__(self):
        return "This coordinate is not implemented for this geometry type."


# define for back compat reasons for code written before yt 4.0
YTUnitOperationError = UnitOperationError


class YTUnitNotRecognized(YTException):
    def __init__(self, unit):
        self.unit = unit

    def __str__(self):
        return f"This dataset doesn't recognize {self.unit!r}"


class YTFieldUnitError(YTException):
    def __init__(self, field_info, returned_units):
        self.msg = (
            f"The field function associated with the field {field_info.name!r} returned "
            f"data with units {returned_units!r} but was defined with units {field_info.units!r}."
        )

    def __str__(self):
        return self.msg


class YTFieldUnitParseError(YTException):
    def __init__(self, field_info):
        self.msg = (
            f"The field {field_info.name!r} has unparsable units {field_info.units!r}."
        )

    def __str__(self):
        return self.msg


class YTSpatialFieldUnitError(YTException):
    def __init__(self, field):
        self.msg = (
            f"Field {field!r} is a spatial field but has unknown units but "
            "spatial fields must have explicitly defined units. Add the "
            "field with explicit 'units' to clear this error."
        )

    def __str__(self):
        return self.msg


class YTHubRegisterError(YTException):
    def __str__(self):
        return (
            "You must create an API key before uploading. See "
            "https://data.yt-project.org/getting_started.html"
        )


class YTNoFilenamesMatchPattern(YTException):
    def __init__(self, pattern):
        self.pattern = pattern

    def __str__(self):
        return f"No filenames were found to match the pattern: {self.pattern!r}"


class YTNoOldAnswer(YTException):
    def __init__(self, path):
        self.path = path

    def __str__(self):
        return f"There is no old answer available.\n{self.path!r}"


class YTNoAnswerNameSpecified(YTException):
    def __init__(self, message=None):
        if message is None or message == "":
            message = (
                "Answer name not provided for the answer testing test."
                "\n  Please specify --answer-name=<answer_name> in"
                " command line mode or in AnswerTestingTest.answer_name"
                " variable."
            )
        self.message = message

    def __str__(self):
        return str(self.message)


class YTCloudError(YTException):
    def __init__(self, path):
        self.path = path

    def __str__(self):
        return (
            f"Failed to retrieve cloud data. Connection may be broken.\n {self.path!r}"
        )


class YTEllipsoidOrdering(YTException):
    def __init__(self, ds, A, B, C):
        YTException.__init__(self, ds=ds)
        self._A = A
        self._B = B
        self._C = C

    def __str__(self):
        return "Must have A>=B>=C"


class EnzoTestOutputFileNonExistent(YTException):
    def __init__(self, filename):
        self.filename = filename
        self.testname = os.path.basename(os.path.dirname(filename))

    def __str__(self):
        return (
            f"Enzo test output file (OutputLog) not generated for: {self.testname!r}.\n"
            "Test did not complete."
        )


class YTNoAPIKey(YTException):
    def __init__(self, service, config_name):
        self.service = service
        self.config_name = config_name

    def __str__(self):
        from yt.config import config_dir

        try:
            conf = os.path.join(config_dir(), "yt", "yt.toml")
        except Exception:
            # this is really not a good time to raise another exception
            conf = "yt's configuration file"
        return f"You need to set an API key for {self.service!r} in {conf} as {self.config_name!r}"


class YTTooManyVertices(YTException):
    def __init__(self, nv, fn):
        self.nv = nv
        self.fn = fn

    def __str__(self):
        s = f"There are too many vertices ({self.nv}) to upload to Sketchfab. "
        s += f"Your model has been saved as {self.fn} .  You should upload manually."
        return s


class YTInvalidWidthError(YTException):
    def __init__(self, width):
        self.error = f"width ({str(width)}) is invalid"

    def __str__(self):
        return str(self.error)


class YTFieldNotParseable(YTException):
    def __init__(self, field):
        self.field = field

    def __str__(self):
        return f"Cannot identify field {self.field}"


class YTDataSelectorNotImplemented(YTException):
    def __init__(self, class_name):
        self.class_name = class_name

    def __str__(self):
        return f"Data selector {self.class_name!r} not implemented."


class YTParticleDepositionNotImplemented(YTException):
    def __init__(self, class_name):
        self.class_name = class_name

    def __str__(self):
        return f"Particle deposition method {self.class_name!r} not implemented."


class YTDomainOverflow(YTException):
    def __init__(self, mi, ma, dle, dre):
        self.mi = mi
        self.ma = ma
        self.dle = dle
        self.dre = dre

    def __str__(self):
        return (
            f"Particle bounds {self.mi} and {self.ma} "
            f"exceed domain bounds {self.dle} and {self.dre}"
        )


class YTIntDomainOverflow(YTException):
    def __init__(self, dims, dd):
        self.dims = dims
        self.dd = dd

    def __str__(self):
        return f"Integer domain overflow: {self.dims} in {self.dd}"


class YTIllDefinedFilter(YTException):
    def __init__(self, filter, s1, s2):
        self.filter = filter
        self.s1 = s1
        self.s2 = s2

    def __str__(self):
        return (
            f"Filter {self.filter!r} ill-defined. "
            f"Applied to shape {self.s1} but is shape {self.s2}."
        )


class YTIllDefinedParticleFilter(YTException):
    def __init__(self, filter, missing):
        self.filter = filter
        self.missing = missing

    def __str__(self):
        msg = (
            '\nThe fields\n\t{},\nrequired by the "{}" particle filter, '
            "are not defined for this dataset."
        )
        f = self.filter
        return msg.format("\n".join(str(m) for m in self.missing), f.name)


class YTIllDefinedBounds(YTException):
    def __init__(self, lb, ub):
        self.lb = lb
        self.ub = ub

    def __str__(self):
        v = f"The bounds {self.lb:0.3e} and {self.ub:0.3e} are ill-defined. "
        v += "Typically this happens when a log binning is specified "
        v += "and zero or negative values are given for the bounds."
        return v


class YTObjectNotImplemented(YTException):
    def __init__(self, ds, obj_name):
        self.ds = ds
        self.obj_name = obj_name

    def __str__(self):
        return f"The object type {self.obj_name!r} is not implemented for the dataset {self.ds!s}"


class YTParticleOutputFormatNotImplemented(YTException):
    def __str__(self):
        return "The particle output format is not supported."


class YTFileNotParseable(YTException):
    def __init__(self, fname, line):
        self.fname = fname
        self.line = line

    def __str__(self):
        return f"Error while parsing file {self.fname!r} at line {self.line}"


class YTRockstarMultiMassNotSupported(YTException):
    def __init__(self, mi, ma, ptype):
        self.mi = mi
        self.ma = ma
        self.ptype = ptype

    def __str__(self):
        v = f"Particle type '{self.ptype}' has minimum mass {self.mi:0.3e} and maximum "
        v += f"mass {self.ma:0.3e}.  Multi-mass particles are not currently supported."
        return v


class YTTooParallel(YTException):
    def __str__(self):
        return "You've used too many processors for this dataset."


class YTElementTypeNotRecognized(YTException):
    def __init__(self, dim, num_nodes):
        self.dim = dim
        self.num_nodes = num_nodes

    def __str__(self):
        return f"Element type not recognized - dim = {self.dim}, num_nodes = {self.num_nodes}"


class YTDuplicateFieldInProfile(Exception):
    def __init__(self, field, new_spec, old_spec):
        self.field = field
        self.new_spec = new_spec
        self.old_spec = old_spec

    def __str__(self):
        r = f"""Field {self.field} already exists with field spec:
               {self.old_spec}
               But being asked to add it with:
               {self.new_spec}"""
        return r


class YTInvalidPositionArray(Exception):
    def __init__(self, shape, dimensions):
        self.shape = shape
        self.dimensions = dimensions

    def __str__(self):
        r = f"""Position arrays must be length and shape (N,3).
               But this one has {self.dimensions} and {self.shape}."""
        return r


class YTIllDefinedCutRegion(Exception):
    def __init__(self, conditions):
        self.conditions = conditions

    def __str__(self):
        r = (
            "Can't mix particle/discrete and fluid/mesh conditions or quantities. "
            "Conditions specified:\n"
        )
        r += "\n".join(c for c in self.conditions)
        return r


class YTMixedCutRegion(Exception):
    def __init__(self, conditions, field):
        self.conditions = conditions
        self.field = field

    def __str__(self):
        r = f"""Can't mix particle/discrete and fluid/mesh conditions or
               quantities.  Field: {self.field} and Conditions specified:
            """
        r += "\n".join(c for c in self.conditions)
        return r


class YTGDFAlreadyExists(Exception):
    def __init__(self, filename):
        self.filename = filename

    def __str__(self):
        return f"A file already exists at {self.filename} and overwrite=False."


class YTNonIndexedDataContainer(YTException):
    def __init__(self, cont):
        self.cont = cont

    def __str__(self):
        class_name = self.cont.__class__.__name__
        return (
            f"The data container type ({class_name}) is an unindexed type. "
            "Operations such as ires, icoords, fcoords and fwidth will not work on it.\n"
            "Did you just attempt to perform an off-axis operation ? "
            "Be sure to consult the latest documentation to see whether the operation "
            "you tried is actually supported for your data type."
        )


class YTGDFUnknownGeometry(Exception):
    def __init__(self, geometry):
        self.geometry = geometry

    def __str__(self):
        return (
            """Unknown geometry %i. Please refer to GDF standard
                  for more information"""
            % self.geometry
        )


class YTInvalidUnitEquivalence(Exception):
    def __init__(self, equiv, unit1, unit2):
        self.equiv = equiv
        self.unit1 = unit1
        self.unit2 = unit2

    def __str__(self):
        return f"The unit equivalence {self.equiv!r} does not exist for the units {self.unit1!r} and {self.unit2!r}."


class YTPlotCallbackError(Exception):
    def __init__(self, callback):
        self.callback = "annotate_" + callback

    def __str__(self):
        return f"{self.callback} callback failed"


class YTUnsupportedPlotCallback(YTPlotCallbackError):
    def __init__(self, callback: str, plot_type: str) -> None:
        super().__init__(callback)
        self.plot_type = plot_type

    def __str__(self):
        return f"The `{self.plot_type}` class currently doesn't support the `{self.callback}` method."


class YTPixelizeError(YTException):
    def __init__(self, message):
        self.message = message

    def __str__(self):
        return self.message


class YTDimensionalityError(YTException):
    def __init__(self, wrong, right):
        self.wrong = wrong
        self.right = right

    def __str__(self):
        return f"Dimensionality specified was {self.wrong} but we need {self.right}"


class YTInvalidShaderType(YTException):
    def __init__(self, source):
        self.source = source

    def __str__(self):
        return f"Can't identify shader_type for file {self.source!r}"


class YTInvalidFieldType(YTException):
    def __init__(self, fields):
        self.fields = fields

    def __str__(self):
        return (
            "\nSlicePlot, ProjectionPlot, and OffAxisProjectionPlot can "
            "only plot fields that\n"
            "are defined on a mesh or for SPH particles, but received the "
            "following N-body\n"
            "particle fields:\n\n"
            f"    {self.fields!r}\n\n"
            "Did you mean to use ParticlePlot or plot a deposited particle "
            "field instead?"
        )


class YTUnknownUniformKind(YTException):
    def __init__(self, kind):
        self.kind = kind

    def __str__(self):
        return f"Can't determine kind specification for {self.kind!r}"


class YTUnknownUniformSize(YTException):
    def __init__(self, size_spec):
        self.size_spec = size_spec

    def __str__(self):
        return f"Can't determine size specification for {self.size_spec!r}"


class YTDataTypeUnsupported(YTException):
    def __init__(self, this, supported):
        self.supported = supported
        self.this = this

    def __str__(self):
        v = f"This operation is not supported for data of geometry {self.this!r}; "
        v += f"It supports data of geometries {self.supported!r}"
        return v


class YTBoundsDefinitionError(YTException):
    def __init__(self, message, bounds):
        self.bounds = bounds
        self.message = message

    def __str__(self):
        v = f"This operation has encountered a bounds error: {self.message} "
        v += f"\nSpecified bounds are {self.bounds!r}."
        return v


def screen_one_element_list(lis):
    if len(lis) == 1:
        return lis[0]
    return lis


class YTIllDefinedProfile(YTException):
    def __init__(self, bin_fields, fields, weight_field, is_pfield):
        nbin = len(bin_fields)
        nfields = len(fields)
        self.bin_fields = screen_one_element_list(bin_fields)
        self.bin_fields_ptype = screen_one_element_list(is_pfield[:nbin])
        self.fields = screen_one_element_list(fields)
        self.fields_ptype = screen_one_element_list(is_pfield[nbin : nbin + nfields])
        self.weight_field = weight_field
        if self.weight_field is not None:
            self.weight_field_ptype = is_pfield[-1]

    def __str__(self):
        msg = (
            "\nCannot create a profile object that mixes particle and mesh "
            "fields.\n\n"
            "Received the following bin_fields:\n\n"
            "   %s, particle_type = %s\n\n"
            "Profile fields:\n\n"
            "   %s, particle_type = %s\n"
        )
        msg = msg % (
            self.bin_fields,
            self.bin_fields_ptype,
            self.fields,
            self.fields_ptype,
        )

        if self.weight_field is not None:
            weight_msg = "\nAnd weight field:\n\n   %s, particle_type = %s\n"
            weight_msg = weight_msg % (self.weight_field, self.weight_field_ptype)
        else:
            weight_msg = ""

        return msg + weight_msg


class YTProfileDataShape(YTException):
    def __init__(self, field1, shape1, field2, shape2):
        self.field1 = field1
        self.shape1 = shape1
        self.field2 = field2
        self.shape2 = shape2

    def __str__(self):
        return (
            "Profile fields must have same shape: {self.field1!r} has "
            f"shape {self.shape1} and {self.field2!r} has shape {self.shape2}."
        )


class YTBooleanObjectError(YTException):
    def __init__(self, bad_object):
        self.bad_object = bad_object

    def __str__(self):
        v = f"Supplied:\n{self.bad_object}\nto a boolean operation"
        v += " but it is not a YTSelectionContainer3D object."
        return v


class YTBooleanObjectsWrongDataset(YTException):
    def __init__(self):
        pass

    def __str__(self):
        return "Boolean data objects must share a common dataset object."


class YTIllDefinedAMR(YTException):
    def __init__(self, level, axis):
        self.level = level
        self.axis = axis

    def __str__(self):
        return (
            f"Grids on the level {self.level} are not properly aligned with cell edges "
            f"on the parent level ({self.axis!r} axis)"
        )


class YTIllDefinedParticleData(YTException):
    pass


class YTIllDefinedAMRData(YTException):
    pass


class YTInconsistentGridFieldShape(YTException):
    def __init__(self, shapes):
        self.shapes = shapes

    def __str__(self):
        msg = "Not all grid-based fields have the same shape!\n"
        for name, shape in self.shapes:
            msg += f"    Field {name!r} has shape {shape}.\n"
        return msg


class YTInconsistentParticleFieldShape(YTException):
    def __init__(self, ptype, shapes):
        self.ptype = ptype
        self.shapes = shapes

    def __str__(self):
        msg = "Not all fields with field type {self.ptype!r} have the same shape!\n"
        for name, shape in self.shapes:
            field = (self.ptype, name)
            msg += f"    Field {field} has shape {shape}.\n"
        return msg


class YTInconsistentGridFieldShapeGridDims(YTException):
    def __init__(self, shapes, grid_dims):
        self.shapes = shapes
        self.grid_dims = grid_dims

    def __str__(self):
        msg = "Not all grid-based fields match the grid dimensions! "
        msg += f"Grid dims are {self.grid_dims}, "
        msg += "and the following fields have shapes that do not match them:\n"
        for name, shape in self.shapes:
            if shape != self.grid_dims:
                msg += f"    Field {name} has shape {shape}.\n"
        return msg


class YTCommandRequiresModule(YTException):
    def __init__(self, module: str):
        self.module = module

    def __str__(self):
        msg = f"This command requires {self.module!r} to be installed.\n\n"
        msg += f"Please install {self.module!r} with the package manager "
        msg += "appropriate for your python environment, e.g.:\n"
        msg += f"  conda install {self.module}\n"
        msg += "or:\n"
        msg += f" python -m pip install {self.module}\n"
        return msg


class YTModuleRemoved(Exception):
    def __init__(self, name, new_home=None, info=None):
        message = f"The {name} module has been removed from yt."
        if new_home is not None:
            message += f"\nIt has been moved to {new_home}."
        if info is not None:
            message += f"\nFor more information, see {info}."
        Exception.__init__(self, message)


class YTArrayTooLargeToDisplay(YTException):
    def __init__(self, size, max_size):
        self.size = size
        self.max_size = max_size

    def __str__(self):
        msg = f"The requested array is of size {self.size}.\n"
        msg += "We do not support displaying arrays larger\n"
        msg += f"than size {self.max_size}."
        return msg


class YTConfigurationError(YTException):
    pass


class GenerationInProgress(Exception):
    def __init__(self, fields):
        self.fields = fields
        super().__init__()


class MountError(Exception):
    def __init__(self, message):
        self.message = message

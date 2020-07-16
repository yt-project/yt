# We don't need to import 'exceptions'
import os.path

from unyt.exceptions import UnitOperationError


class YTException(Exception):
    def __init__(self, message=None, ds=None):
        Exception.__init__(self, message)
        self.ds = ds


# Data access exceptions:


class YTOutputNotIdentified(YTException):
    def __init__(self, args, kwargs):
        self.args = args
        self.kwargs = kwargs

    def __str__(self):
        return "Supplied %s %s, but could not load!" % (self.args, self.kwargs)


class YTSphereTooSmall(YTException):
    def __init__(self, ds, radius, smallest_cell):
        YTException.__init__(self, ds=ds)
        self.radius = radius
        self.smallest_cell = smallest_cell

    def __str__(self):
        return "%0.5e < %0.5e" % (self.radius, self.smallest_cell)


class YTAxesNotOrthogonalError(YTException):
    def __init__(self, axes):
        self.axes = axes

    def __str__(self):
        return "The supplied axes are not orthogonal.  %s" % (self.axes)


class YTNoDataInObjectError(YTException):
    def __init__(self, obj):
        self.obj_type = getattr(obj, "_type_name", "")

    def __str__(self):
        s = "The object requested has no data included in it."
        if self.obj_type == "slice":
            s += "  It may lie on a grid face.  Try offsetting slightly."
        return s


class YTFieldNotFound(YTException):
    def __init__(self, fname, ds):
        self.fname = fname
        self.ds = ds

    def __str__(self):
        return "Could not find field '%s' in %s." % (self.fname, self.ds)


class YTParticleTypeNotFound(YTException):
    def __init__(self, fname, ds):
        self.fname = fname
        self.ds = ds

    def __str__(self):
        return "Could not find particle_type '%s' in %s." % (self.fname, self.ds)


class YTSceneFieldNotFound(YTException):
    pass


class YTCouldNotGenerateField(YTFieldNotFound):
    def __str__(self):
        return "Could field '%s' in %s could not be generated." % (self.fname, self.ds)


class YTFieldTypeNotFound(YTException):
    def __init__(self, ftype, ds=None):
        self.ftype = ftype
        self.ds = ds

    def __str__(self):
        if self.ds is not None and self.ftype in self.ds.particle_types:
            return (
                "Could not find field type '%s'.  "
                + "This field type is a known particle type for this dataset.  "
                + "Try adding this field with particle_type=True."
            ) % self.ftype
        else:
            return "Could not find field type '%s'." % (self.ftype)


class YTSimulationNotIdentified(YTException):
    def __init__(self, sim_type):
        YTException.__init__(self)
        self.sim_type = sim_type

    def __str__(self):
        return "Simulation time-series type %s not defined." % self.sim_type


class YTCannotParseFieldDisplayName(YTException):
    def __init__(self, field_name, display_name, mathtext_error):
        self.field_name = field_name
        self.display_name = display_name
        self.mathtext_error = mathtext_error

    def __str__(self):
        return (
            'The display name "%s" '
            "of the derived field %s "
            "contains the following LaTeX parser errors:\n"
        ) % (self.display_name, self.field_name) + self.mathtext_error


class YTCannotParseUnitDisplayName(YTException):
    def __init__(self, field_name, unit_name, mathtext_error):
        self.field_name = field_name
        self.unit_name = unit_name
        self.mathtext_error = mathtext_error

    def __str__(self):
        return (
            'The unit display name "%s" '
            "of the derived field %s "
            "contains the following LaTeX parser errors:\n"
        ) % (self.unit_name, self.field_name) + self.mathtext_error


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
        return "dataset %s is missing %s parameter." % (self.ds, self.parameter)


class NoStoppingCondition(YTException):
    def __init__(self, ds):
        YTException.__init__(self, ds=ds)

    def __str__(self):
        return (
            "Simulation %s has no stopping condition.  StopTime or StopCycle should be set."
            % self.ds
        )


class YTNotInsideNotebook(YTException):
    def __str__(self):
        return "This function only works from within an IPython Notebook."


class YTGeometryNotSupported(YTException):
    def __init__(self, geom):
        self.geom = geom

    def __str__(self):
        return "We don't currently support %s geometry" % self.geom


class YTCoordinateNotImplemented(YTException):
    def __str__(self):
        return "This coordinate is not implemented for this geometry type."


# define for back compat reasons for code written before yt 4.0
YTUnitOperationError = UnitOperationError


class YTUnitNotRecognized(YTException):
    def __init__(self, unit):
        self.unit = unit

    def __str__(self):
        return "This dataset doesn't recognize %s" % self.unit


class YTFieldUnitError(YTException):
    def __init__(self, field_info, returned_units):
        self.msg = (
            "The field function associated with the field '%s' returned "
            "data with units '%s' but was defined with units '%s'."
        )
        self.msg = self.msg % (field_info.name, returned_units, field_info.units)

    def __str__(self):
        return self.msg


class YTFieldUnitParseError(YTException):
    def __init__(self, field_info):
        self.msg = "The field '%s' has unparseable units '%s'."
        self.msg = self.msg % (field_info.name, field_info.units)

    def __str__(self):
        return self.msg


class YTSpatialFieldUnitError(YTException):
    def __init__(self, field):
        msg = (
            "Field '%s' is a spatial field but has unknown units but "
            "spatial fields must have explicitly defined units. Add the "
            "field with explicit 'units' to clear this error."
        )
        self.msg = msg % (field,)

    def __str__(self):
        return self.msg


class YTHubRegisterError(YTException):
    def __str__(self):
        return (
            "You must create an API key before uploading.  See "
            + "https://data.yt-project.org/getting_started.html"
        )


class YTNoFilenamesMatchPattern(YTException):
    def __init__(self, pattern):
        self.pattern = pattern

    def __str__(self):
        return "No filenames were found to match the pattern: " + "'%s'" % (
            self.pattern
        )


class YTNoOldAnswer(YTException):
    def __init__(self, path):
        self.path = path

    def __str__(self):
        return "There is no old answer available.\n" + str(self.path)


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
        return "Failed to retrieve cloud data. Connection may be broken.\n" + str(
            self.path
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
            "Enzo test output file (OutputLog) not generated for: "
            + "'%s'" % (self.testname)
            + ".\nTest did not complete."
        )


class YTNoAPIKey(YTException):
    def __init__(self, service, config_name):
        self.service = service
        self.config_name = config_name

    def __str__(self):
        return "You need to set an API key for %s in ~/.config/yt/ytrc as %s" % (
            self.service,
            self.config_name,
        )


class YTTooManyVertices(YTException):
    def __init__(self, nv, fn):
        self.nv = nv
        self.fn = fn

    def __str__(self):
        s = "There are too many vertices (%s) to upload to Sketchfab. " % (self.nv)
        s += "Your model has been saved as %s .  You should upload manually." % (
            self.fn
        )
        return s


class YTInvalidWidthError(YTException):
    def __init__(self, width):
        self.error = "width (%s) is invalid" % str(width)

    def __str__(self):
        return str(self.error)


class YTFieldNotParseable(YTException):
    def __init__(self, field):
        self.field = field

    def __str__(self):
        return "Cannot identify field %s" % (self.field,)


class YTDataSelectorNotImplemented(YTException):
    def __init__(self, class_name):
        self.class_name = class_name

    def __str__(self):
        return "Data selector '%s' not implemented." % (self.class_name)


class YTParticleDepositionNotImplemented(YTException):
    def __init__(self, class_name):
        self.class_name = class_name

    def __str__(self):
        return "Particle deposition method '%s' not implemented." % (self.class_name)


class YTDomainOverflow(YTException):
    def __init__(self, mi, ma, dle, dre):
        self.mi = mi
        self.ma = ma
        self.dle = dle
        self.dre = dre

    def __str__(self):
        return "Particle bounds %s and %s exceed domain bounds %s and %s" % (
            self.mi,
            self.ma,
            self.dle,
            self.dre,
        )


class YTIntDomainOverflow(YTException):
    def __init__(self, dims, dd):
        self.dims = dims
        self.dd = dd

    def __str__(self):
        return "Integer domain overflow: %s in %s" % (self.dims, self.dd)


class YTIllDefinedFilter(YTException):
    def __init__(self, filter, s1, s2):
        self.filter = filter
        self.s1 = s1
        self.s2 = s2

    def __str__(self):
        return "Filter '%s' ill-defined.  Applied to shape %s but is shape %s." % (
            self.filter,
            self.s1,
            self.s2,
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
        return msg.format("\n".join([str(m) for m in self.missing]), f.name)


class YTIllDefinedBounds(YTException):
    def __init__(self, lb, ub):
        self.lb = lb
        self.ub = ub

    def __str__(self):
        v = "The bounds %0.3e and %0.3e are ill-defined. " % (self.lb, self.ub)
        v += "Typically this happens when a log binning is specified "
        v += "and zero or negative values are given for the bounds."
        return v


class YTObjectNotImplemented(YTException):
    def __init__(self, ds, obj_name):
        self.ds = ds
        self.obj_name = obj_name

    def __str__(self):
        v = r"The object type '%s' is not implemented for the dataset "
        v += r"'%s'."
        return v % (self.obj_name, self.ds)


class YTParticleOutputFormatNotImplemented(YTException):
    def __str__(self):
        return "The particle output format is not supported."


class YTFileNotParseable(YTException):
    def __init__(self, fname, line):
        self.fname = fname
        self.line = line

    def __str__(self):
        v = r"Error while parsing file %s at line %s"
        return v % (self.fname, self.line)


class YTRockstarMultiMassNotSupported(YTException):
    def __init__(self, mi, ma, ptype):
        self.mi = mi
        self.ma = ma
        self.ptype = ptype

    def __str__(self):
        v = "Particle type '%s' has minimum mass %0.3e and maximum " % (
            self.ptype,
            self.mi,
        )
        v += "mass %0.3e.  Multi-mass particles are not currently supported." % (
            self.ma
        )
        return v


class YTTooParallel(YTException):
    def __str__(self):
        return "You've used too many processors for this dataset."


class YTElementTypeNotRecognized(YTException):
    def __init__(self, dim, num_nodes):
        self.dim = dim
        self.num_nodes = num_nodes

    def __str__(self):
        return "Element type not recognized - dim = %s, num_nodes = %s" % (
            self.dim,
            self.num_nodes,
        )


class YTDuplicateFieldInProfile(Exception):
    def __init__(self, field, new_spec, old_spec):
        self.field = field
        self.new_spec = new_spec
        self.old_spec = old_spec

    def __str__(self):
        r = """Field %s already exists with field spec:
               %s
               But being asked to add it with:
               %s""" % (
            self.field,
            self.old_spec,
            self.new_spec,
        )
        return r


class YTInvalidPositionArray(Exception):
    def __init__(self, shape, dimensions):
        self.shape = shape
        self.dimensions = dimensions

    def __str__(self):
        r = """Position arrays must be length and shape (N,3).
               But this one has %s and %s.""" % (
            self.dimensions,
            self.shape,
        )
        return r


class YTIllDefinedCutRegion(Exception):
    def __init__(self, conditions):
        self.conditions = conditions

    def __str__(self):
        r = """Can't mix particle/discrete and fluid/mesh conditions or
               quantities.  Conditions specified:
            """
        r += "\n".join([c for c in self.conditions])
        return r


class YTMixedCutRegion(Exception):
    def __init__(self, conditions, field):
        self.conditions = conditions
        self.field = field

    def __str__(self):
        r = """Can't mix particle/discrete and fluid/mesh conditions or
               quantities.  Field: %s and Conditions specified:
            """ % (
            self.field,
        )
        r += "\n".join([c for c in self.conditions])
        return r


class YTGDFAlreadyExists(Exception):
    def __init__(self, filename):
        self.filename = filename

    def __str__(self):
        return "A file already exists at %s and overwrite=False." % self.filename


class YTNonIndexedDataContainer(YTException):
    def __init__(self, cont):
        self.cont = cont

    def __str__(self):
        return (
            "The data container (%s) is an unindexed type.  "
            "Operations such as ires, icoords, fcoords and fwidth "
            "will not work on it." % type(self.cont)
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
        return (
            "The unit equivalence '%s' does not exist for the units '%s' and '%s'."
            % (self.equiv, self.unit1, self.unit2)
        )


class YTPlotCallbackError(Exception):
    def __init__(self, callback, error):
        self.callback = "annotate_" + callback
        self.error = error

    def __str__(self):
        msg = "%s callback failed with the following error: %s"
        return msg % (self.callback, self.error)


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
        return "Dimensionality specified was %s but we need %s" % (
            self.wrong,
            self.right,
        )


class YTInvalidShaderType(YTException):
    def __init__(self, source):
        self.source = source

    def __str__(self):
        return "Can't identify shader_type for file '%s.'" % (self.source)


class YTInvalidFieldType(YTException):
    def __init__(self, fields):
        self.fields = fields

    def __str__(self):
        msg = (
            "\nSlicePlot, ProjectionPlot, and OffAxisProjectionPlot can "
            "only plot fields that\n"
            "are defined on a mesh or for SPH particles, but received the "
            "following N-body\n"
            "particle fields:\n\n"
            "    %s\n\n"
            "Did you mean to use ParticlePlot or plot a deposited particle "
            "field instead?" % self.fields
        )
        return msg


class YTUnknownUniformKind(YTException):
    def __init__(self, kind):
        self.kind = kind

    def __str__(self):
        return "Can't determine kind specification for %s" % (self.kind)


class YTUnknownUniformSize(YTException):
    def __init__(self, size_spec):
        self.size_spec = size_spec

    def __str__(self):
        return "Can't determine size specification for %s" % (self.size_spec)


class YTDataTypeUnsupported(YTException):
    def __init__(self, this, supported):
        self.supported = supported
        self.this = this

    def __str__(self):
        v = "This operation is not supported for data of geometry %s; " % self.this
        v += "It supports data of geometries %s" % (self.supported,)
        return v


class YTBoundsDefinitionError(YTException):
    def __init__(self, message, bounds):
        self.bounds = bounds
        self.message = message

    def __str__(self):
        v = "This operation has encountered a bounds error: "
        v += self.message
        v += " Specified bounds are '%s'." % (self.bounds,)
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
            "Profile fields must have same shape: %s has "
            + "shape %s and %s has shape %s."
        ) % (self.field1, self.shape1, self.field2, self.shape2)


class YTBooleanObjectError(YTException):
    def __init__(self, bad_object):
        self.bad_object = bad_object

    def __str__(self):
        v = "Supplied:\n%s\nto a boolean operation" % (self.bad_object)
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
        msg = (
            "Grids on the level {} are not properly aligned with cell edges "
            "on the parent level ({} axis)"
        ).format(self.level, self.axis)
        return msg


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
            msg += "    Field {} has shape {}.\n".format(name, shape)
        return msg


class YTInconsistentParticleFieldShape(YTException):
    def __init__(self, ptype, shapes):
        self.ptype = ptype
        self.shapes = shapes

    def __str__(self):
        msg = ("Not all fields with field type '{}' have the same shape!\n").format(
            self.ptype
        )
        for name, shape in self.shapes:
            field = (self.ptype, name)
            msg += "    Field {} has shape {}.\n".format(field, shape)
        return msg


class YTInconsistentGridFieldShapeGridDims(YTException):
    def __init__(self, shapes, grid_dims):
        self.shapes = shapes
        self.grid_dims = grid_dims

    def __str__(self):
        msg = "Not all grid-based fields match the grid dimensions! "
        msg += "Grid dims are {}, ".format(self.grid_dims)
        msg += "and the following fields have shapes that do not match them:\n"
        for name, shape in self.shapes:
            if shape != self.grid_dims:
                msg += "    Field {} has shape {}.\n".format(name, shape)
        return msg


class YTCommandRequiresModule(YTException):
    def __init__(self, module):
        self.module = module

    def __str__(self):
        msg = 'This command requires "%s" to be installed.\n\n' % self.module
        msg += 'Please install "%s" with the package manager ' % self.module
        msg += "appropriate for your python environment, e.g.:\n"
        msg += "  conda install %s\n" % self.module
        msg += "or:\n"
        msg += "  pip install %s\n" % self.module
        return msg


class YTModuleRemoved(Exception):
    def __init__(self, name, new_home=None, info=None):
        message = "The %s module has been removed from yt." % name
        if new_home is not None:
            message += "\nIt has been moved to %s." % new_home
        if info is not None:
            message += "\nFor more information, see %s." % info
        Exception.__init__(self, message)


class YTArrayTooLargeToDisplay(YTException):
    def __init__(self, size, max_size):
        self.size = size
        self.max_size = max_size

    def __str__(self):
        msg = "The requested array is of size %s.\n" % self.size
        msg += "We do not support displaying arrays larger\n"
        msg += "than size %s." % self.max_size
        return msg

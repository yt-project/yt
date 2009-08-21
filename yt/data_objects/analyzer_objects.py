class AnalysisTask(object):

    def __init__(self, *args, **kwargs):
        # This should only get called if the subclassed object
        # does not override
        if len(args) + len(kwargs) != len(self._params):
            raise RuntimeError
        self.__dict__.update(zip(self._params, args))
        self.__dict__.update(kwargs)

    def __repr__(self):
        # Stolen from AMRData.__repr__
        s = "%s: " % (self.__class__.__name__, self._analysis_type)
        s += ", ".join(["%s=%s" % (i, getattr(self,i))
                       for i in self._params])
        return s

class MaximumValue(AnalysisTask):
    _params = ['field']

    def eval(self, data_object):
        v = data_object.quantities["MaxLocation"](
                self.field, lazy_reader=True)[0]
        return v

class ParameterValue(AnalysisTask):
    _params = ['parameter']

    def __init__(self, parameter, cast=None):
        self.parameter = parameter
        if cast is None:
            cast = lambda a: a
        self.cast = cast

    def eval(self, pf):
        return self.cast(pf.get_parameter(parameter))

class CurrentTimeYears(AnalysisTask):
    _params = []

    def eval(self, pf):
        return pf["Time"] * pf["years"]

class SliceDataset(AnalysisTask):
    _params = ['field', 'axis']

    def eval(self, pf):
        pass

class SlicePlotDataset(AnalysisTask):
    _params = ['field', 'axis', 'center']

    def __init__(self, *args, **kwargs):
        import yt.raven
        self.raven = yt.raven
        AnalysisTask.__init__(self, *args, **kwargs)

    def eval(self, pf):
        pc = self.raven.PlotCollection(pf, center = self.center)
        pc.add_slice(self.field, self.axis)
        return pc.save()[0]

import numpy as np
import os

from functools import wraps

from yt.config import ytcfg
from yt.convenience import load
from yt.testing import assert_allclose
from yt.utilities.answer_testing.framework import \
    AnswerTestingTest, \
    can_run_ds, \
    FieldValuesTest, \
    GridValuesTest, \
    ProjectionValuesTest, \
    temp_cwd


class AssertWrapper(object):
    """
    Used to wrap a numpy testing assertion, in order to provide a useful name
    for a given assertion test.
    """
    def __init__(self, description, *args):
        # The key here is to add a description attribute, which nose will pick
        # up.
        self.args = args
        self.description = description

    def __call__(self):
        self.args[0](*self.args[1:])


def requires_outputlog(path = ".", prefix = ""):
    def ffalse(func):
        return lambda: None
    def ftrue(func):
        @wraps(func)
        def fyielder(*args, **kwargs):
            with temp_cwd(path):
                for t in func(*args, **kwargs):
                    if isinstance(t, AnswerTestingTest):
                        t.prefix = prefix
                    yield t
        return fyielder
    if os.path.exists("OutputLog"):
        return ftrue
    with temp_cwd(path):
        if os.path.exists("OutputLog"):
            return ftrue
    return ffalse
     
def standard_small_simulation(ds_fn, fields):
    if not can_run_ds(ds_fn): return
    dso = [None]
    tolerance = ytcfg.getint("yt", "answer_testing_tolerance")
    bitwise = ytcfg.getboolean("yt", "answer_testing_bitwise")
    for field in fields:
        if bitwise:
            yield GridValuesTest(ds_fn, field)
        if 'particle' in field: continue
        for dobj_name in dso:
            for axis in [0, 1, 2]:
                for weight_field in [None, "Density"]:
                    yield ProjectionValuesTest(
                        ds_fn, axis, field, weight_field,
                        dobj_name, decimals=tolerance)
            yield FieldValuesTest(
                    ds_fn, field, dobj_name, decimals=tolerance)
                    
class ShockTubeTest(object):
    def __init__(self, data_file, solution_file, fields, 
                 left_edges, right_edges, rtol, atol):
        self.solution_file = solution_file
        self.data_file = data_file
        self.fields = fields
        self.left_edges = left_edges
        self.right_edges = right_edges
        self.rtol = rtol
        self.atol = atol

    def __call__(self):
        # Read in the ds
        ds = load(self.data_file)
        ds.setup_deprecated_fields()
        exact = self.get_analytical_solution()

        ad = ds.all_data()
        position = ad['x']
        for k in self.fields:
            field = ad[k].d
            for xmin, xmax in zip(self.left_edges, self.right_edges):
                mask = (position >= xmin)*(position <= xmax)
                exact_field = np.interp(position[mask], exact['pos'], exact[k])
                myname = f"ShockTubeTest_{k}"
                # yield test vs analytical solution 
                yield AssertWrapper(myname, assert_allclose, field[mask], 
                                    exact_field, self.rtol, self.atol)

    def get_analytical_solution(self):
        # Reads in from file 
        pos, dens, vel, pres, inte = \
                np.loadtxt(self.solution_file, unpack=True)
        exact = {}
        exact['pos'] = pos
        exact['Density'] = dens
        exact['x-velocity'] = vel
        exact['Pressure'] = pres
        exact['ThermalEnergy'] = inte
        return exact

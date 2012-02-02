import matplotlib
matplotlib.use('Agg')
from yt.config import ytcfg
ytcfg["yt", "loglevel"] = "50"
ytcfg["yt", "serialize"] = "False"

from yt.utilities.answer_testing.api import \
    RegressionTestRunner, clear_registry, create_test, \
    TestFieldStatistics, TestAllProjections, registry_entries, \
    Xunit
from yt.utilities.command_line import get_yt_version

from yt.mods import *
import fnmatch
import imp
import optparse
import itertools
import time

#
# We assume all tests are to be run, unless explicitly given the name of a
# single test or something that can be run through fnmatch.
#
# Keep in mind that we use a different nomenclature here than is used in the
# Enzo testing system.  Our 'tests' are actually tests that are small and that
# run relatively quickly on a single dataset; in Enzo's system, a 'test'
# encompasses both the creation and the examination of data.  Here we assume
# the data is kept constant.
#

cwd = os.path.dirname(globals().get("__file__", os.getcwd()))


def load_tests(iname, idir):
    f, filename, desc = imp.find_module(iname, [idir])
    tmod = imp.load_module(iname, f, filename, desc)
    return tmod


def find_and_initialize_tests():
    mapping = {}
    for f in glob.glob(os.path.join(cwd, "*.py")):
        clear_registry()
        iname = os.path.basename(f[:-3])
        try:
            load_tests(iname, cwd)
            mapping[iname] = registry_entries()
            #print "Associating %s with" % (iname)
            #print "\n    ".join(registry_entries())
        except ImportError:
            pass
    return mapping

if __name__ == "__main__":
    clear_registry()
    mapping = find_and_initialize_tests()
    test_storage_directory = ytcfg.get("yt", "test_storage_dir")
    try:
        my_hash = get_yt_version()
    except:
        my_hash = "UNKNOWN%s" % (time.time())
    parser = optparse.OptionParser()
    parser.add_option("-f", "--parameter-file", dest="parameter_file",
        default=os.path.join(cwd, "DD0010/moving7_0010"),
        help="The parameter file value to feed to 'load' to test against")
    parser.add_option("-l", "--list", dest="list_tests", action="store_true",
        default=False, help="List all tests and then exit")
    parser.add_option("-t", "--tests", dest="test_pattern", default="*",
        help="The test name pattern to match.  Can include wildcards.")
    parser.add_option("-o", "--output", dest="storage_dir",
        default=test_storage_directory,
        help="Base directory for storing test output.")
    parser.add_option("-c", "--compare", dest="compare_name",
        default=None,
        help="The name against which we will compare")
    parser.add_option("-n", "--name", dest="this_name",
        default=my_hash,
        help="The name we'll call this set of tests")
    opts, args = parser.parse_args()

    if opts.list_tests:
        tests_to_run = []
        for m, vals in mapping.items():
            new_tests = fnmatch.filter(vals, opts.test_pattern)
            if len(new_tests) == 0: continue
            load_tests(m, cwd)
            keys = set(registry_entries())
            tests_to_run += [t for t in new_tests if t in keys]
        tests = list(set(tests_to_run))
        print "\n    ".join(tests)
        sys.exit(0)

    # Load the test pf and make sure it's good.
    pf = load(opts.parameter_file)
    if pf is None:
        print "Couldn't load the specified parameter file."
        sys.exit(1)

    # Now we modify our compare name and self name to include the pf.
    compare_id = opts.compare_name
    watcher = None
    if compare_id is not None:
        compare_id += "_%s_%s" % (pf, pf._hash())
        watcher = Xunit()
    this_id = opts.this_name + "_%s_%s" % (pf, pf._hash())

    rtr = RegressionTestRunner(this_id, compare_id,
                               results_path=opts.storage_dir,
                               compare_results_path=opts.storage_dir,
                               io_log=[opts.parameter_file])

    rtr.watcher = watcher
    tests_to_run = []
    for m, vals in mapping.items():
        new_tests = fnmatch.filter(vals, opts.test_pattern)

        if len(new_tests) == 0: continue
        load_tests(m, cwd)
        keys = set(registry_entries())
        tests_to_run += [t for t in new_tests if t in keys]
    for test_name in sorted(tests_to_run):
        print "RUNNING TEST", test_name
        rtr.run_test(test_name)
    if watcher is not None:
        rtr.watcher.report()
    failures = 0
    passes = 1
    for test_name, result in sorted(rtr.passed_tests.items()):
        print "TEST %s: %s" % (test_name, result)
        if result: passes += 1
        else: failures += 1
    print "Number of passes  : %s" % passes
    print "Number of failures: %s" % failures

import os, shelve, cPickle, sys
import yt.cmdln as cmdln
from output_tests import test_registry, MultipleOutputTest, \
                         RegressionTestException

class RegressionTestStorage(object):
    def __init__(self, id, path = "."):
        self.id = id
        self._path = os.path.join(path, "results_%s" % self.id)
        if not os.path.isdir(self._path): os.mkdir(self._path)
        if os.path.isfile(self._path): raise RuntimeError

    def _fn(self, tn):
        return os.path.join(self._path, tn)

    def __setitem__(self, test_name, result):
        # We have to close our shelf manually,
        # as the destructor does not necessarily do this.
        # Context managers would be more appropriate.
        f = open(self._fn(test_name), "wb")
        cPickle.dump(result, f, protocol=-1)
        f.close()

    def __getitem__(self, test_name):
        f = open(self._fn(test_name), "rb")
        tr = cPickle.load(f)
        f.close()
        return tr

class RegressionTestRunner(object):
    def __init__(self, id, compare_id = None,
                 results_path = ".", io_log = "OutputLog"):
        # This test runner assumes it has been launched with the current
        # working directory that of the test case itself.
        self.io_log = io_log
        self.id = id
        if compare_id is not None:
            self.old_results = RegressionTestStorage(
                                    compare_id, path=results_path)
        else:
            self.old_results = None
        self.results = RegressionTestStorage(id, path=results_path)
        self.plot_list = {}
        self.passed_tests = {}

    def run_all_tests(self):
        plot_list = []
        for i,name in enumerate(sorted(test_registry)):
            self.run_test(name)
        return plot_list

    def run_test(self, name):
        # We'll also need to call the "compare" operation,
        # but for that we'll need a data store.
        test = test_registry[name]
        plot_list = []
        if test.output_type == 'single':
            mot = MultipleOutputTest(self.io_log)
            for i,fn in enumerate(mot):
                # This next line is to keep the shelve module
                # from happily gobbling the disk
                #if i > 5: break 
                test_instance = test(fn)
                test_instance.name = "%s_%s" % (
                    os.path.basename(fn), test_instance.name )
                self._run(test_instance)

        elif test.output_type == 'multiple':
            test_instance = test(self.io_log)
            self._run(test_instance)

    def _run(self, test):
        print self.id, "Running", test.name,
        test.setup()
        test.run()
        self.plot_list[test.name] = test.plot()
        self.results[test.name] = test.result
        success = self._compare(test)
        if success == True: print "SUCCEEDED"
        else: print "FAILED"
        self.passed_tests[test.name] = success

    def _compare(self, test):
        if self.old_results is None:
            return True
        old_result = self.old_results[test.name]
        try:
            test.compare(old_result)
        except RegressionTestException, exc:
            return str(exc)
        return True

    def run_tests_from_file(self, filename):
        for line in open(filename):
            test_name = line.strip()
            if test_name not in test_registry:
                if test_name[0] != "#":
                    print "Test '%s' not recognized, skipping" % (test_name)
                continue
            print "Running '%s'" % (test_name)
            self.run_test(line.strip())

def run():
    # This should be made to work with the optparse library
    if sys.argv[-1] == "-f":
        first_runner = RegressionTestRunner("first")
        first_runner.run_all_tests()
    else:
        second_runner = RegressionTestRunner("second", "first")
        second_runner.run_all_tests()

class EnzoTestRunnerCommands(cmdln.Cmdln):
    name = "enzo_tests"

    def do_store(self, subcmd, opts, name):
        """
        Run and store a new dataset.

        ${cmd_option_list}
        """
        test_runner = RegressionTestRunner(name)
        test_runner.run_all_tests()

    def do_compare(self, subcmd, opts, reference, comparison):
        """
        Compare a reference dataset against a new dataset.  The new dataset
        will be run regardless of whether it exists or not.

        ${cmd_option_list}
        """
        test_runner = RegressionTestRunner(comparison, reference)
        test_runner.run_all_tests()

if __name__ == "__main__":
    etrc = EnzoTestRunnerCommands()
    sys.exit(etrc.main())

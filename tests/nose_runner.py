import argparse
import glob
import multiprocessing
import os
import shutil
import sys

import nose
import numpy
import yaml

from yt.config import ytcfg
from yt.utilities.answer_testing.framework import AnswerTesting

numpy.set_printoptions(threshold=5, edgeitems=1, precision=4)


class NoseWorker(multiprocessing.Process):
    def __init__(self, task_queue, result_queue):
        multiprocessing.Process.__init__(self)
        self.task_queue = task_queue
        self.result_queue = result_queue

    def run(self):
        proc_name = self.name
        while True:
            next_task = self.task_queue.get()
            if next_task is None:
                print(f"{proc_name}: Exiting")
                self.task_queue.task_done()
                break
            print(f"{proc_name}: {next_task}")
            result = next_task()
            self.task_queue.task_done()
            self.result_queue.put(result)
            if next_task.exclusive:
                print(f"{proc_name}: Exiting (exclusive)")
                break
        return


class NoseTask:
    def __init__(self, job, store_missing=False):
        argv, exclusive = job
        self.argv = argv
        self.name = argv[0]
        self.exclusive = exclusive
        self.store_missing = store_missing

    def _store_files(self):
        # The store lives at <answers_dir>/<name>/<name>[.db|.dat...] (shelve, backend-dependent
        # suffix). Test for the file, not the directory: the answer-testing plugin creates the
        # directory before any test runs, so a failed store run leaves an empty one behind.
        return glob.glob(os.path.join(self.answers_dir, self.name, self.name + "*"))

    def _store_answers(self):
        """Run the store pass for a group that has no store. Returns one of:

        "stored"   a usable store now exists (only with store_missing)
        "no-data"  the pass passed but wrote nothing: a dataset in the group is missing (the
                   plugin marks the results tainted and skips the dump), so there is nothing
                   to compare against until the data is seeded
        "failed"   the pass failed, or it generated a store without store_missing

        Anything but "stored" leaves nothing behind: the compare pass opens the store
        read-write and would create an empty one, and every later run would then compare
        against it.
        """
        print(
            f"{self.name}: no stored answers under {self.answers_dir}; "
            "running the store pass (--answer-store)",
            flush=True,
        )
        ok = nose.run(
            argv=self.argv + ["--answer-store"],
            addplugins=[AnswerTesting()],
            exit=False,
        )
        if ok and self._store_files() and self.store_missing:
            return "stored"
        if ok and self._store_files():
            outcome = "failed"
            why = (
                "generated a store, but store_missing is off; discarding it. Dispatch with "
                "store_missing=true to accept this ref's answers as the baseline"
            )
        elif ok:
            outcome = "no-data"
            why = "no store written: a dataset in the group is missing; skipping the group"
        else:
            outcome = "failed"
            why = "tests failed during the store pass; discarding the partial store"
        print(f"{self.name}: {why}", flush=True)
        shutil.rmtree(os.path.join(self.answers_dir, self.name), ignore_errors=True)
        return outcome

    def __call__(self):
        test_dir = ytcfg.get("yt", "test_data_dir")
        self.answers_dir = os.path.join(test_dir, "answers")
        answer_test = "--with-answer-testing" in self.argv
        stored = False
        if answer_test and not self._store_files():
            outcome = self._store_answers()
            if outcome != "stored":
                return (self.name, outcome, False)
            stored = True
        if os.path.isfile(f"{self.name}.xml"):
            os.remove(f"{self.name}.xml")
        if answer_test:
            print(f"{self.name}: comparing against stored answers", flush=True)
        success = nose.run(argv=self.argv, addplugins=[AnswerTesting()], exit=False)
        return (self.name, "passed" if success else "failed", stored)

    def __str__(self):
        return f"WILL DO self.name = {self.name}"


def parse_args(argv=None):
    parser = argparse.ArgumentParser(
        description="Run the nose answer-test groups in tests/tests.yaml in parallel."
    )
    parser.add_argument(
        "--groups",
        default=os.environ.get("YT_TEST_GROUPS", ""),
        help="comma-separated answer_tests keys to run (default: all; "
        "env fallback YT_TEST_GROUPS)",
    )
    parser.add_argument(
        "--unittests",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="run other_tests (the unittests job); default is yes unless --groups is given",
    )
    parser.add_argument(
        "--store-missing",
        action=argparse.BooleanOptionalAction,
        default=os.environ.get("YT_STORE_MISSING", "").lower() in ("1", "true", "yes"),
        help="keep the answers generated for groups that have no store, making them the "
        "baseline for later runs; default off (env fallback YT_STORE_MISSING), and such "
        "groups then fail",
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=int(os.environ.get("NUM_WORKERS", 6)),
        help="parallel nose processes (default: env NUM_WORKERS, else 6)",
    )
    args = parser.parse_args(argv)
    args.groups = [g.strip() for g in args.groups.split(",") if g.strip()]
    if args.unittests is None:
        args.unittests = not args.groups
    return args


def generate_tasks_input(groups=(), unittests=True):
    pyver = f"py{sys.version_info.major}{sys.version_info.minor}"
    test_dir = ytcfg.get("yt", "test_data_dir")
    answers_dir = os.path.join(test_dir, "answers")
    with open("tests/tests.yaml") as f:
        tests = yaml.load(f, Loader=yaml.FullLoader)

    all_answer_tests = tests["answer_tests"]
    if groups:
        unknown = sorted(set(groups) - set(all_answer_tests))
        if unknown:
            sys.exit(f"unknown answer-test groups: {unknown}")
        answer_tests = {k: all_answer_tests[k] for k in groups}
    else:
        answer_tests = all_answer_tests
    other_tests = tests["other_tests"] if unittests else {}

    base_argv = ["-s", "--nologcapture", "--with-xunit"]

    base_answer_argv = [
        f"--local-dir={answers_dir}",
        "--with-answer-testing",
        "--answer-big-data",
        "--local",
    ]

    args = []

    for test in list(other_tests.keys()):
        args.append(([test] + base_argv + other_tests[test], True))
    for answer in list(answer_tests.keys()):
        if answer_tests[answer] is None:
            continue
        argv = [f"{pyver}_{answer}"]
        argv += base_argv + base_answer_argv
        argv.append(f"--answer-name={argv[0]}")
        argv += answer_tests[answer]
        args.append((argv, False))

    # unittests must exclude every answer test, not just the groups selected for this run
    exclude_answers = []
    for key in all_answer_tests:
        for t in all_answer_tests[key] or []:
            exclude_answers.append(t.replace(".py:", ".").replace("/", "."))
    exclude_answers = [f"--exclude-test={ex}" for ex in exclude_answers]

    args = [
        (
            (item + [f"--xunit-file={item[0]}.xml"], exclusive)
            if item[0] != "unittests"
            else (item + ["--xunit-file=unittests.xml"] + exclude_answers, exclusive)
        )
        for item, exclusive in args
    ]
    return args


if __name__ == "__main__":
    # multiprocessing.log_to_stderr(logging.DEBUG)
    opts = parse_args()
    tasks = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()

    jobs = generate_tasks_input(opts.groups, opts.unittests)
    num_jobs = len(jobs)
    # no point spawning more workers than there are jobs
    num_consumers = min(opts.num_workers, num_jobs)
    print(f"{num_jobs} nose run(s) across {num_consumers} worker(s)", flush=True)
    consumers = [NoseWorker(tasks, results) for i in range(num_consumers)]
    for w in consumers:
        w.start()

    for job in jobs:
        if job[1]:
            num_consumers -= 1  # take into account exclusive jobs
        tasks.put(NoseTask(job, opts.store_missing))

    for _i in range(num_consumers):
        tasks.put(None)

    tasks.join()

    failed, stored, skipped = [], [], []
    while num_jobs:
        name, status, did_store = results.get()
        num_jobs -= 1
        if status == "failed":
            failed.append(name)
        elif status == "no-data":
            skipped.append(name)
        if did_store:
            stored.append(name)

    # Newly generated stores are the baseline from here on, and a skipped group is one the run
    # says nothing about; make both visible to whoever dispatched the run, in the log and in
    # the GitHub job summary when there is one.
    lines = []
    if stored:
        lines.append(f"STORED new answers ({len(stored)}): {' '.join(sorted(stored))}")
    if skipped:
        lines.append(
            f"SKIPPED, datasets missing ({len(skipped)}): {' '.join(sorted(skipped))}"
        )
    for line in lines:
        print(line)
    summary = os.environ.get("GITHUB_STEP_SUMMARY")
    if summary and lines:
        with open(summary, "a") as f:
            f.writelines(f"**{line}**\n" for line in lines)
    if skipped and opts.groups:
        # A skip in a full run is a known data gap. A skip in a group somebody asked for by
        # name, typically a new group whose data is not seeded yet, is not a pass.
        print("requested groups were skipped for missing datasets; failing the run")
        failed.extend(skipped)
    if failed:
        print(f"FAILED ({len(failed)}): {' '.join(sorted(failed))}")
        sys.exit(1)
    print("all nose runs passed")

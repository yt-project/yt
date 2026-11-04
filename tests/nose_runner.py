import argparse
import multiprocessing
import os
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
    def __init__(self, job):
        argv, exclusive = job
        self.argv = argv
        self.name = argv[0]
        self.exclusive = exclusive

    def __call__(self):
        test_dir = ytcfg.get("yt", "test_data_dir")
        answers_dir = os.path.join(test_dir, "answers")
        if "--with-answer-testing" in self.argv and not os.path.isdir(
            os.path.join(answers_dir, self.name)
        ):
            nose.run(
                argv=self.argv + ["--answer-store"],
                addplugins=[AnswerTesting()],
                exit=False,
            )
        if os.path.isfile(f"{self.name}.xml"):
            os.remove(f"{self.name}.xml")
        success = nose.run(argv=self.argv, addplugins=[AnswerTesting()], exit=False)
        return (self.name, success)

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

    num_consumers = opts.num_workers
    consumers = [NoseWorker(tasks, results) for i in range(num_consumers)]
    for w in consumers:
        w.start()

    num_jobs = 0
    for job in generate_tasks_input(opts.groups, opts.unittests):
        if job[1]:
            num_consumers -= 1  # take into account exclusive jobs
        tasks.put(NoseTask(job))
        num_jobs += 1

    for _i in range(num_consumers):
        tasks.put(None)

    tasks.join()

    failed = []
    while num_jobs:
        name, success = results.get()
        num_jobs -= 1
        if not success:
            failed.append(name)

    if failed:
        print(f"FAILED ({len(failed)}): {' '.join(sorted(failed))}")
        sys.exit(1)
    print("all nose runs passed")

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
        nose.run(argv=self.argv, addplugins=[AnswerTesting()], exit=False)
        return ""

    def __str__(self):
        return f"WILL DO self.name = {self.name}"


def generate_tasks_input():
    pyver = f"py{sys.version_info.major}{sys.version_info.minor}"
    test_dir = ytcfg.get("yt", "test_data_dir")
    answers_dir = os.path.join(test_dir, "answers")
    tests = yaml.load(open("tests/tests.yaml"), Loader=yaml.FullLoader)

    base_argv = ["-s", "--nologcapture", "--with-xunit"]

    base_answer_argv = [
        f"--local-dir={answers_dir}",
        "--with-answer-testing",
        "--answer-big-data",
        "--local",
    ]

    args = []

    for test in list(tests["other_tests"].keys()):
        args.append(([test] + base_argv + tests["other_tests"][test], True))
    for answer in list(tests["answer_tests"].keys()):
        if tests["answer_tests"][answer] is None:
            continue
        argv = [f"{pyver}_{answer}"]
        argv += base_argv + base_answer_argv
        argv.append(f"--answer-name={argv[0]}")
        argv += tests["answer_tests"][answer]
        args.append((argv, False))

    exclude_answers = []
    answer_tests = tests["answer_tests"]
    for key in answer_tests:
        for t in answer_tests[key]:
            exclude_answers.append(t.replace(".py:", ".").replace("/", "."))
    exclude_answers = [f"--exclude-test={ex}" for ex in exclude_answers]

    args = [
        (item + [f"--xunit-file={item[0]}.xml"], exclusive)
        if item[0] != "unittests"
        else (item + ["--xunit-file=unittests.xml"] + exclude_answers, exclusive)
        for item, exclusive in args
    ]
    return args


if __name__ == "__main__":
    # multiprocessing.log_to_stderr(logging.DEBUG)
    tasks = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()

    num_consumers = int(os.environ.get("NUM_WORKERS", 6))
    consumers = [NoseWorker(tasks, results) for i in range(num_consumers)]
    for w in consumers:
        w.start()

    num_jobs = 0
    for job in generate_tasks_input():
        if job[1]:
            num_consumers -= 1  # take into account exclusive jobs
        tasks.put(NoseTask(job))
        num_jobs += 1

    for _i in range(num_consumers):
        tasks.put(None)

    tasks.join()

    while num_jobs:
        result = results.get()
        num_jobs -= 1

import sys
import os
import yaml
import multiprocessing
import nose
from yt.extern.six import StringIO
from yt.config import ytcfg
from yt.utilities.answer_testing.framework import AnswerTesting

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
                print("%s: Exiting" % proc_name)
                self.task_queue.task_done()
                break
            print('%s: %s' % (proc_name, next_task))
            result = next_task()
            self.task_queue.task_done()
            self.result_queue.put(result)
        return

class NoseTask(object):
    def __init__(self, argv):
        self.argv = argv
        self.name = argv[0]

    def __call__(self):
        old_stderr = sys.stderr
        sys.stderr = mystderr = StringIO()
        test_dir = ytcfg.get("yt", "test_data_dir")
        answers_dir = os.path.join(test_dir, "answers")
        if '--with-answer-testing' in self.argv and \
                not os.path.isdir(os.path.join(answers_dir, self.name)):
            nose.run(argv=self.argv + ['--answer-store'],
                     addplugins=[AnswerTesting()], exit=False)
        if os.path.isfile("{}.xml".format(self.name)):
            os.remove("{}.xml".format(self.name))
        nose.run(argv=self.argv, addplugins=[AnswerTesting()], exit=False)
        sys.stderr = old_stderr
        return mystderr.getvalue()

    def __str__(self):
        return 'WILL DO self.name = %s' % self.name


def generate_tasks_input():
    pyver = "py{}{}".format(sys.version_info.major, sys.version_info.minor)
    if sys.version_info < (3, 0, 0):
        DROP_TAG = "py3"
    else:
        DROP_TAG = "py2"

    test_dir = ytcfg.get("yt", "test_data_dir")
    answers_dir = os.path.join(test_dir, "answers")
    with open('tests/tests.yaml', 'r') as obj:
        lines = obj.read()
    data = '\n'.join([line for line in lines.split('\n')
                      if DROP_TAG not in line])
    tests = yaml.load(data)

    base_argv = ['--local-dir=%s' % answers_dir, '-v',
                 '--with-answer-testing', '--answer-big-data', '--local']
    args = []

    for test in list(tests["other_tests"].keys()):
        args.append([test] + tests["other_tests"][test])
    for answer in list(tests["answer_tests"].keys()):
        if tests["answer_tests"][answer] is None:
            continue
        argv = ["{}_{}".format(pyver, answer)]
        argv += base_argv
        argv.append('--answer-name=%s' % argv[0])
        argv += tests["answer_tests"][answer]
        args.append(argv)

    args = [item + ['-s', '--nologcapture', '--xunit-file=%s.xml' % item[0]]
            for item in args]
    return args

if __name__ == "__main__":
    # multiprocessing.log_to_stderr(logging.DEBUG)
    tasks = multiprocessing.JoinableQueue()
    results = multiprocessing.Queue()

    num_consumers = 6  # TODO 
    consumers = [NoseWorker(tasks, results) for i in range(num_consumers)]
    for w in consumers:
        w.start()

    num_jobs = 0
    for job in generate_tasks_input():
        tasks.put(NoseTask(job))
        num_jobs += 1

    for i in range(num_consumers):
        tasks.put(None)

    tasks.join()

    while num_jobs:
        result = results.get()
        print(result)
        num_jobs -= 1

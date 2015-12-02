import sys
import os
import yaml
import multiprocessing as mp
import nose
from yt.config import ytcfg
from yt.utilities.answer_testing.framework import AnswerTesting

test_dir = ytcfg.get("yt", "test_data_dir")
answers_dir = os.path.join(test_dir, "answers")

with open('tests/tests_%i.%i.yaml' % sys.version_info[:2], 'r') as obj:
    tests = yaml.load(obj)

base_argv = ['fake', '--local-dir=%s' % answers_dir, '-v', '-s', '--nologcapture',
             '--with-answer-testing', '--answer-big-data', '--local']

def run_job(answer):
    global base_argv, answers_dir, tests
    argv = []
    argv += base_argv
    argv.append('--xunit-file=%s.xml' % answer)
    argv.append('--answer-name=%s' % answer)
    argv += tests[answer]
    if not os.path.isdir(os.path.join(answers_dir, answer)):
        nose.run(argv=argv + ['--answer-store'],
                 addplugins=[AnswerTesting()], exit=False)
    nose.run(argv=argv, addplugins=[AnswerTesting()], exit=False)

pool = mp.Pool(processes=mp.cpu_count())
results = pool.map(run_job, list(tests.keys()))

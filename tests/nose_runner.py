import sys
import os
import yaml
import multiprocessing as mp
import nose
import glob
from contextlib import closing
from yt.config import ytcfg
from yt.utilities.answer_testing.framework import AnswerTesting


def run_job(argv):
    with closing(open(str(os.getpid()) + ".out", "w")) as fstderr:
        cur_stderr = sys.stderr
        sys.stderr = fstderr
        answer = argv[0]
        test_dir = ytcfg.get("yt", "test_data_dir")
        answers_dir = os.path.join(test_dir, "answers")
        if not os.path.isdir(os.path.join(answers_dir, answer)):
            nose.run(argv=argv + ['--answer-store'],
                     addplugins=[AnswerTesting()], exit=False)
        nose.run(argv=argv, addplugins=[AnswerTesting()], exit=False)
    sys.stderr = cur_stderr

if __name__ == "__main__":
    test_dir = ytcfg.get("yt", "test_data_dir")
    answers_dir = os.path.join(test_dir, "answers")
    with open('tests/tests_%i.%i.yaml' % sys.version_info[:2], 'r') as obj:
        tests = yaml.load(obj)

    base_argv = ['--local-dir=%s' % answers_dir, '-v', '-s', '--nologcapture',
                 '--with-answer-testing', '--answer-big-data', '--local']
    args = [['unittests', '-v', '-s', '--nologcapture']]
    for answer in list(tests.keys()):
        argv = [answer]
        argv += base_argv
        argv.append('--xunit-file=%s.xml' % answer)
        argv.append('--answer-name=%s' % answer)
        argv += tests[answer]
        args.append(argv)
    
    processes = [mp.Process(target=run_job, args=(args[i],))
                 for i in range(len(args))]
    for p in processes:
        p.start()
    for p in processes:
        p.join(timeout=7200)
        if p.is_alive():
            p.terminate()
            p.join(timeout=30)
    for fname in glob.glob("*.out"):
        with open(fname, 'r') as fin:
            print(fin.read())
        os.remove(fname)

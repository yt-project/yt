from yt.enki import *

d = {'wd':'/usr/work/mturk/t', 'exe':'./enzo_red_i9_r16', 'restart':True, 'nproc':16, 'logFile':'test.log', 'parameterFile':'RedshiftOutput0000'}

hostRed.spawnProcess(d)

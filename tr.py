from yt.testing import *
pf = fake_random_pf(64)
p1 = [0.1, 0.2, 0.3]
p2 = [0.8, 0.1, 0.4]

ray = pf.ray(p1, p2)
dd = ray["dts"].sum()
print dd

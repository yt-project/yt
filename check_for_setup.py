import os

for dirpath, dirnames, filenames in os.walk("yt"):
    if "setup.py" not in filenames: print dirpath

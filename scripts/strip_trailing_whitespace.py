import re
from pathlib import Path

yt_dir = Path(__file__).parent.parent
for file in yt_dir.glob("yt/**/*.py"):
    matches = None
    with open(file, mode="r") as fileobj:
        lines = fileobj.readlines()
    with open(file, mode="w") as fileobj:
        striped_lines = [L.rstrip() for L in lines]
        fileobj.write("\n".join(striped_lines))
        if striped_lines and striped_lines[-1]:
            fileobj.write("\n")

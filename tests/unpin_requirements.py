import configparser
import io
import re
from collections.abc import MutableMapping

PINNED_VERSION_REGEXP = re.compile(r",?(<|<=|==)([0-9a-z]+\.?)+")


def unpin(s: str) -> str:
    return re.sub(PINNED_VERSION_REGEXP, "", s)


def unpin_mapping(m: MutableMapping, key) -> None:
    reqs = m[key].split("\n")
    m[key] = "\n".join(unpin(_) for _ in reqs)


if __name__ == "__main__":
    cp = configparser.ConfigParser()
    cp.read("setup.cfg")
    unpin_mapping(cp["options"], "install_requires")
    for target in cp["options.extras_require"].keys():
        unpin_mapping(cp["options.extras_require"], target)

    output = io.StringIO()
    cp.write(output)
    s = output.getvalue().replace("\t", " " * 4)
    s = "\n".join(_.rstrip() for _ in s.split("\n"))
    with open("setup_out.cfg", "w") as fh:
        fh.write(s)

# /// script
# requires-python = ">=3.10"
# dependencies = [
#     "tomli ; python_full_version < '3.11'",
#     "tomli-w",
# ]
# ///
import re
import sys

import tomli_w

if sys.version_info >= (3, 11):
    import tomllib
else:
    import tomli as tomllib

PINNED_VERSION_REGEXP = re.compile(r",?(<|<=|==)([0-9a-z]+\.?)+")


def unpin_requirements(requirements: list[str]) -> list[str]:
    return [re.sub(PINNED_VERSION_REGEXP, "", _) for _ in requirements]


if __name__ == "__main__":
    with open("pyproject.toml", "rb") as fr:
        config = tomllib.load(fr)

    config["project"]["dependencies"] = unpin_requirements(
        config["project"]["dependencies"]
    )
    for key, reqs in config["project"]["optional-dependencies"].items():
        config["project"]["optional-dependencies"][key] = unpin_requirements(reqs)

    with open("pyproject_out.toml", "wb") as fw:
        tomli_w.dump(config, fw)

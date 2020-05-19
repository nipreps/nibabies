#!/usr/bin/env python
"""Read current version."""
import sys
import os.path as op


def main():
    sys.path.insert(0, op.abspath("."))
    from nibabies._version import get_versions

    print(get_versions()["version"])


if __name__ == "__main__":
    main()
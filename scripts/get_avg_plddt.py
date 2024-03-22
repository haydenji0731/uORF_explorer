#!/usr/bin/env python

import json
import sys


def main(json_fn):
    f = open(json_fn)
    data = json.load(f)
    scores = list()
    for i in data['plddt']:
        scores.append(float(i))
    avg = sum(scores) / len(scores)
    print(avg)


if __name__ == "__main__":
    main(sys.argv[1])

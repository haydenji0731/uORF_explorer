#!/usr/bin/env python

import sys


def main(fn):
    fh = open(fn, 'r')
    out_fh = open("./reformatted.tsv", 'w')
    for line in fh:
        fields = line.strip().split("\t")
        for i in range(len(fields)):
            field = fields[i]
            if i == 0:
                tid = field
            else:
                out_fh.write(tid + "\t" + field + "\n")
    out_fh.close()
    fh.close()


if __name__ == "__main__":
    main(sys.argv[1])

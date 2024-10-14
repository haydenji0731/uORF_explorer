#!/usr/bin/env python

import sys

def load_tid(fn):
    tids = set()
    with open(fn, 'r') as f:
        for ln in f:
            tid = ln.strip()
            assert tid not in tids
            tids.add(tid)
    print(f'{len(tids)} transcript IDs loaded.')
    return tids

def get_tid(s):
    tmp = s.split(';')
    tid = None
    for x in tmp:
        kv_pair = x.strip().split(" ")
        if len(kv_pair) < 2:
            continue
        if kv_pair[0].strip() == 'transcript_id':
            tid = kv_pair[1].strip().replace('"', '')
            break
    return tid

def filter_gtf(fn, tids, out_fn):
    out_fh = open(out_fn, 'w')
    with open(fn, 'r') as f:
        for ln in f:
            clean_ln = ln.strip()
            fields = clean_ln.split("\t")
            feature = fields[2]
            if feature == "transcript":
                skip = False
                tid = get_tid(fields[8])
                if not tid: return -1
                if tid not in tids:
                    skip = True
                if not skip:
                    out_fh.write(ln)
            else:
                if not skip:
                    out_fh.write(ln)
    out_fh.close()
                
def main(in_fn, tid_fn, out_fn):
    tids = load_tid(tid_fn)
    res = filter_gtf(in_fn, tids, out_fn)
    if res == -1:
        print(f'An error occurred while loading the input GTF. Check format.')
    
if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])
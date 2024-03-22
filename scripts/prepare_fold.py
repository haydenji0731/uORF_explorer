#!/usr/bin/env python

import os
import pyfastx
import sys
import ast
from subprocess import call


def load_island(fn):
    fh = open(fn, 'r')
    pid_tbl = dict()
    for line in fh:
        fields = line.strip().split("\t")
        pid = fields[0]
        is_tx = ast.literal_eval(fields[1])
        pid_tbl[pid] = is_tx
    return pid_tbl


def write2fa(name, seq, out_fh):
    out_fh.write(">" + name + "\n")
    for i in range(0, len(seq), 60):
        out_fh.write(seq[i:i + 60] + "\n")


def mod_ref_fa(ref_prot_fn, out_dir):
    ref_prot_fa = pyfastx.Fasta(ref_prot_fn)
    tmp_fn = os.path.join(out_dir, "tmp.fa")
    tmp_fh = open(tmp_fn, 'w')
    for obj in ref_prot_fa:
        tid = obj.name.split("|")[1]
        write2fa(tid, obj.seq, tmp_fh)
    tmp_fh.close()
    return tmp_fn


def main(island_fn, ref_prot_fn, tx_prot_fn, out_dir, n):
    if not os.path.exists(out_dir):
        print("out directory not present. creating it...")
        os.makedirs(out_dir)
    pid_tbl = load_island(island_fn)

    mod_fn = mod_ref_fa(ref_prot_fn, out_dir)

    ref_prot_fa = pyfastx.Fasta(mod_fn)
    tx_prot_fa = pyfastx.Fasta(tx_prot_fn)

    ctr = 0
    fnum = 0
    fold_fn = "for_fold." + str(fnum) + ".fa"
    fold_fn = os.path.join(out_dir, fold_fn)
    fold_fh = open(fold_fn, 'w')
    for pid in pid_tbl:
        is_tx = pid_tbl[pid]
        if is_tx:
            pseq = tx_prot_fa[pid].seq
        else:
            pseq = ref_prot_fa[pid].seq

        if ctr < n:
            write2fa(pid, pseq, fold_fh)
            ctr += 1
        else:
            fold_fh.close()
            fnum += 1
            fold_fn = "for_fold." + str(fnum) + ".fa"
            fold_fn = os.path.join(out_dir, fold_fn)
            fold_fh = open(fold_fn, 'w')
            write2fa(pid, pseq, fold_fh)
            ctr = 1
    fold_fh.close()

    rm_cmd = "rm " + mod_fn
    print(rm_cmd)
    call(rm_cmd, shell=True)
    rm_cmd = "rm " + mod_fn + ".fxi"
    print(rm_cmd)
    call(rm_cmd, shell=True)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4], int(sys.argv[5]))

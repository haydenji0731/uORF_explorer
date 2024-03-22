#!/usr/bin/env python

import argparse
import pyfastx
from extract_novel_seq_tmp import load_dup_map, collate_uniq_info, extract_novel_cds_seq
from tqdm import tqdm
import sys


def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-a1", "--annot-1", type=str, help="", required=True)
    parser.add_argument("-a2", "--annot-2", type=str, help="", required=True)
    parser.add_argument("-m", "--map", type=str, help="", required=True)
    parser.add_argument("-pt", "--prot", type=str, help="", required=True)
    parser.add_argument("-op", "--out-prefix", type=str, help="", required=True)
    args = parser.parse_args()
    return args


def load_annot(fn):
    fh = open(fn, 'r')
    junc_tbl = dict()
    elen_tbl = dict()
    t2r_map = dict()
    for ln in tqdm(fh):
        clean_ln = ln.strip()
        fields = clean_ln.split("\t")
        if fields[2] == "transcript":
            strand = fields[6]
            kv_pairs = fields[8].split(";")
            for p in kv_pairs:
                tmp = p.strip().split(" ")
                if len(tmp) < 2:
                    continue
                k = tmp[0]
                if k == "transcript_id":
                    tid = tmp[1].replace('"', '')
                elif k == "reference_id":
                    rid = tmp[1].replace('"', '')
                    t2r_map[tid] = rid
                elif k == "splice":
                    # e1's last bp (upstream)
                    up_pos = tmp[2].strip().replace('(', '').replace(',', '')
                    up_pos = int(up_pos)
                    # e2's first bp - 1; correct if SPLAM is ON.
                    # e1 - e2 in forward; e2 - e1 in reverse
                    down_pos = tmp[3].strip().replace(')', '')
                    down_pos = int(down_pos)
                    junc = (up_pos, down_pos)
                    junc_tbl[tid] = (up_pos, down_pos)
                    elen = 0
                    break
            skip = False
        # fine novel exon length by iterating through exons
        elif fields[2] == "exon":
            if skip:
                continue
            l_pos = int(fields[3])
            r_pos = int(fields[4])
            if strand == '-':
                if l_pos == junc[1] + 1:
                    elen += r_pos - l_pos + 1
                    elen_tbl[tid] = elen
                    skip = True
                else:
                    elen += r_pos - l_pos + 1
            elif strand == '+':
                if r_pos == junc[0]:
                    elen += r_pos - l_pos + 1
                    elen_tbl[tid] = elen
                    skip = True
                else:
                    elen += r_pos - l_pos + 1
            else:
                print("Invalid strand information. Exiting...")
                sys.exit()
    fh.close()
    assert len(junc_tbl) == len(elen_tbl)
    # each tid paired up with a single reference
    return t2r_map, junc_tbl, elen_tbl


def main():
    args = get_args()
    # load gtfs
    t2r_map_1, junc_tbl_1, elen_tbl_1 = load_annot(args.annot_1)
    t2r_map_2, junc_tbl_2, elen_tbl_2 = load_annot(args.annot_2)

    # load dup mappings
    u2d_map = load_dup_map(args.map)
    u_tbl = collate_uniq_info(t2r_map_1, t2r_map_2, junc_tbl_1, junc_tbl_2, elen_tbl_1, elen_tbl_2, u2d_map)

    fn = args.out_prefix + "_novel_exons.tsv"
    fh = open(fn, 'w')
    for u_tid in u_tbl:
        tmp = u_tbl[u_tid]
        for rid, junc, elen, clen in tmp:
            fh.write(u_tid + "\t" + rid + "\t" + str(junc[0]) + "\t" + str(junc[1]) +
                     "\t" + str(elen) + "\t" + str(clen) + "\n")
    fh.close()

    fn = args.out_prefix + "_novel_CDSs.tsv"
    fh = open(fn, 'w')
    fa = pyfastx.Fasta(args.prot)
    cseq_tbl = extract_novel_cds_seq(u_tbl, fa)
    for u_tid in cseq_tbl:
        u_cseq_tbl = cseq_tbl[u_tid]
        for cseq in u_cseq_tbl:
            rid_lst = u_cseq_tbl[cseq]
            fh.write(u_tid + "\t" + cseq + "\t" + str(len(cseq)))
            for rid in rid_lst:
                fh.write("\t" + rid)
            fh.write("\n")
    fh.close()


if __name__ == "__main__":
    main()

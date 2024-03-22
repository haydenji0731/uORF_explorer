#!/usr/bin/env python

import argparse
import sys
import pyfastx
from tqdm import tqdm


def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-t2r", "--tx-to-ref-tsv", type=str, help="", required=True)
    parser.add_argument("-ref", "--ref-fa", type=str, help="", required=True)
    parser.add_argument("-a", "--annot", type=str, help="", required=True)
    parser.add_argument("-p", "--prot", type=str, help="", required=True)
    parser.add_argument("-op", "--out-prefix", type=str, help="", required=True)
    args = parser.parse_args()
    return args


def load_t2r_map(fn):
    fh = open(fn, 'r')
    t2r_lst = list()
    ctr = 0
    for ln in tqdm(fh):
        ctr += 1
        clean_ln = ln.strip()
        fields = clean_ln.split("\t")
        tid = fields[0].strip()
        rid = fields[1].strip()
        t2r_lst.append((tid, rid))
    assert ctr == len(t2r_lst)
    return t2r_lst


def load_tx_gtf(fn):
    fh = open(fn, 'r')
    # (tx, junc)
    junc_tbl = dict()
    # (tx, novel_exon_len)
    novel_elen_tbl = dict()
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
                elif k == "sites":
                    # e1 last bp
                    up_pos = tmp[1].strip().replace('(', '').replace(',', '')
                    up_pos = int(up_pos)
                    # e2 first bp - 1
                    # e1 - e2 in forward; e2 - e1 in reverse
                    down_pos = tmp[2].strip().replace(')', '')
                    down_pos = int(down_pos)
                    # assume all duplicates use the same junction
                    junc = (up_pos, down_pos)
                    junc_tbl[tid] = (up_pos, down_pos)
                    novel_elen = 0
                skip = False
        # find novel exon length
        elif fields[2] == "exon":
            if skip:
                continue
            l_pos = int(fields[3])
            r_pos = int(fields[4])
            if strand == '-':
                # SPLAM
                if l_pos == junc[1]:
                # TieBrush
                #if l_pos == junc[1] + 1:
                    novel_elen += r_pos - l_pos + 1
                    novel_elen_tbl[tid] = novel_elen
                    skip = True
                else:
                    novel_elen += r_pos - l_pos + 1
            elif strand == '+':
                if r_pos == junc[0]:
                    novel_elen += r_pos - l_pos + 1
                    novel_elen_tbl[tid] = novel_elen
                    skip = True
                else:
                    novel_elen += r_pos - l_pos + 1
            else:
                print("Invalid strand information. Exiting...")
                sys.exit()
    assert len(junc_tbl) == len(novel_elen_tbl)
    return junc_tbl, novel_elen_tbl


def find_lcs(s1, s2):
    lcs = ""
    min_len = min(len(s1), len(s2))

    for i in range(1, min_len + 1):
        if s1[-i:] == s2[-i:]:
            lcs = s1[-i:]

    return lcs


def get_novel_cds(t2r_set, ref_p_tbl, tx_p_tbl):
    novel_cds_tbl = dict()
    ident_prot = set()
    for tid, rid in tqdm(t2r_set):
        tid = tid.lstrip('\ufeff')
        ref_seq = ref_p_tbl[rid].seq
        tx_seq = tx_p_tbl[tid].seq
        lcs = find_lcs(ref_seq, tx_seq)
        novel_cds = tx_seq[:len(tx_seq) - len(lcs)]
        if novel_cds == "":
            ident_prot.add((tid, rid))
        novel_cds_tbl[(tid, rid)] = novel_cds
    return novel_cds_tbl, ident_prot


def main():
    args = get_args()
    # returns pyfastx seq object
    print("loading reference protein sequences")
    ref_p_tbl = pyfastx.Fasta(args.ref_fa)
    print("loading uorft-encoded protein sequences")
    tx_p_tbl = pyfastx.Fasta(args.prot)
    # sys.exit()
    print("loading tx-to-ref mappings")
    t2r_set = load_t2r_map(args.tx_to_ref_tsv)
    print("loading uorft annotation gtf")
    junc_tbl, novel_elen_tbl = load_tx_gtf(args.annot)
    print("writing out junction and novel exon lengths")
    etc_fn = args.out_prefix + "_tid_junc_elen.tsv"
    etc_fh = open(etc_fn, 'w')
    for tid in tqdm(junc_tbl):
        junc = junc_tbl[tid]
        novel_elen = novel_elen_tbl[tid]
        etc_fh.write(tid + "\t" + str(junc[0]) + "\t" + str(junc[1]) + "\t" + str(novel_elen) + "\n")
    etc_fh.close()
    print("getting novel CDS sequences")
    novel_cds_tbl, ident_prot = get_novel_cds(t2r_set, ref_p_tbl, tx_p_tbl)
    main_fn = args.out_prefix + "_novel_cds.tsv"
    main_fh = open(main_fn, 'w')
    print("writing out novel CDS sequences")
    for tid, rid in novel_cds_tbl:
        novel_cds = novel_cds_tbl[(tid, rid)]
        if novel_cds == "":
            novel_cds = "None"
        main_fh.write(tid + "\t" + rid + "\t" + novel_cds + "\n")
    main_fh.close()
    ident_fn = args.out_prefix + "_ident_prots.tsv"
    ident_fh = open(ident_fn, 'w')
    for tid, rid in ident_prot:
        ident_fh.write(tid + "\t" + rid + "\n")
    ident_fh.close()


if __name__ == "__main__":
    main()


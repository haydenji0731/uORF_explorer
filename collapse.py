#!/usr/bin/env python

import sys
import pyfastx
from subprocess import call


def filter_and_renum_gtf(gtf_fn, out_fn, tx_lst):
    gtf_fh = open(gtf_fn, 'r')
    out_fh = open(out_fn, 'w')
    tx_index = dict()
    ctr = 0
    for ln in gtf_fh:
        fields = ln.split("\t")
        feature = fields[2]
        if feature == "transcript":
            keep_print = False
            infos = fields[8].split(";")
            infos = infos[:-1]
            for info in infos:
                info_clean = info.strip()
                kv_pair = info_clean.split(" ")
                k = kv_pair[0]
                v = kv_pair[1].replace('"', '')
                if k == "transcript_id":
                    tid = v
                    if tid in tx_lst:
                        ctr += 1
                        tx_index[tid] = ctr
                        out_ln = renum_tx(ln, ctr)
                        out_fh.write(out_ln)
                        keep_print = True
                    break
        elif keep_print:
            out_ln = renum_tx(ln, ctr)
            out_fh.write(out_ln)
    gtf_fh.close()
    out_fh.close()
    print("total # of transcripts after collapse: " + str(ctr))
    return tx_index


def renum_tx(ln, ctr):
    fields = ln.strip().split("\t")
    out_ln = ""
    for i in range(len(fields)):
        if i < 8:
            out_ln += fields[i]
            out_ln += "\t"
        else:
            infos = fields[i].split(";")
            infos = infos[:-1]
            for j in range(len(infos)):
                info = infos[j]
                info_clean = info.strip()
                kv_pair = info_clean.split(" ")
                k = kv_pair[0]
                if k == "transcript_id":
                    out_ln = out_ln + 'transcript_id "uorft_denovo_' + str(ctr) + '"; '
                elif j == len(infos) - 1:
                    out_ln = out_ln + info + ";\n"
                else:
                    out_ln = out_ln + info + "; "
    return out_ln


# def collapse(fn):
#     prev_seq = None
#     prev_name = None
#     uniq_lst = list()
#     dup_tbl = dict()
#     for name, seq in pyfastx.Fasta(fn, build_index=False):
#         if prev_seq is None:
#             prev_seq = seq
#             prev_name = name
#             uniq_lst.append(name)
#             dup_tbl[name] = list()
#         else:
#             if seq == prev_seq:
#                 dup_tbl[prev_name].append(name)
#             else:
#                 prev_seq = seq
#                 prev_name = name
#                 uniq_lst.append(name)
#                 dup_tbl[name] = list()
#     return uniq_lst, dup_tbl


def collapse(fn):
    seq_tbl = dict()
    for name, seq in pyfastx.Fasta(fn, build_index=False):
        if seq not in seq_tbl:
            seq_tbl[seq] = [name]
        else:
            seq_tbl[seq].append(name)
    dup_tbl = dict()
    uniq_lst = list()
    for seq in seq_tbl:
        name_lst = seq_tbl[seq]
        key = name_lst[0]
        uniq_lst.append(key)
        for i in range(1, len(name_lst)):
            name = name_lst[i]
            dup_tbl[key].append(name)
    return uniq_lst, dup_tbl


def build_ref_tbl(gtf_fn):
    gtf_fh = open(gtf_fn, 'r')
    ref_tbl = dict()
    for ln in gtf_fh:
        fields = ln.split("\t")
        feature = fields[2]
        if feature == "transcript":
            infos = fields[8].split(";")
            infos = infos[:-1]
            for info in infos:
                info_clean = info.strip()
                kv_pair = info_clean.split(" ")
                k = kv_pair[0]
                v = kv_pair[1].replace('"', '')
                if k == "transcript_id":
                    tid = v
                elif k == "reference_id":
                    rid = v
                    ref_tbl[tid] = rid
                    break
    return ref_tbl


def write_dup_map(dup_tbl, ref_tbl, u2d_fn, u2dr_fn, dr2u_fn, tx_index):
    u2d_fh = open(u2d_fn, 'w')
    u2dr_fh = open(u2dr_fn, 'w')
    dr2u_fh = open(dr2u_fn, 'w')
    for tid in dup_tbl:
        new_tid = "uorft_denovo_" + str(tx_index[tid])
        rid_set = set()
        rid = ref_tbl[tid]
        rid_set.add(rid)
        dup_lst = dup_tbl[tid]
        if len(dup_lst) == 0:
            u2d_fh.write(new_tid + "\n")
        else:
            u2d_fh.write(new_tid + "\t")
            for i in range(len(dup_lst)):
                did = dup_lst[i]
                rid = ref_tbl[did]
                rid_set.add(rid)
                if i == len(dup_lst) - 1:
                    u2d_fh.write(did + "\n")
                else:
                    u2d_fh.write(did + "\t")
        u2dr_fh.write(new_tid)
        for rid in rid_set:
            u2dr_fh.write("\t" + rid)
            dr2u_fh.write(rid + "\t" + new_tid + "\n")
        u2dr_fh.write("\n")


def write2fa(name, seq, out_fh):
    out_fh.write(">" + name + "\n")
    for i in range(len(seq)):
        out_fh.write(seq[i:i+70] + "\n")


def main_collapse(nt_fn, prot_fn, gtf_fn, genome_fn, out_prefix):
    ref_tbl = build_ref_tbl(gtf_fn)
    nt_uniq_lst, nt_dup_tbl = collapse(nt_fn)
    prot_uniq_lst, prot_dup_tbl = collapse(prot_fn)
    is_edge = False
    if len(nt_uniq_lst) == len(prot_uniq_lst):
        print("no edge case to handle. continue...")
    else:
        print("handling a couple of edge cases...")
        is_edge = True
        to_del_set = set()
        consensus_tbl = dict()
        for nt_tid in nt_uniq_lst:
            nt_dup_lst = nt_dup_tbl[nt_tid]
            prot_dup_lst = prot_dup_tbl[nt_tid]
            if len(nt_dup_lst) != len(prot_dup_lst):
                if len(nt_dup_lst) > len(prot_dup_lst):
                    print("this can't happen")
                    sys.exit()
                diff = set(prot_dup_lst) - set(nt_dup_lst)
                consensus_tbl[nt_tid] = diff
                to_del_set |= diff

    filt_gtf_fn = out_prefix + "_dup_removed.gtf"
    tx_index = filter_and_renum_gtf(gtf_fn, filt_gtf_fn, nt_uniq_lst)
    filt_nt_fn = out_prefix + "_dup_removed_txs.fa"
    filt_prot_fn = out_prefix + "_dup_removed_prots.fa"

    gf_cmd = "gffread -w " + filt_nt_fn + " -g " + genome_fn + " " + filt_gtf_fn
    print(gf_cmd)
    call(gf_cmd, shell=True)

    gf_cmd = "gffread -S -y " + filt_prot_fn + " -g " + genome_fn + " " + filt_gtf_fn
    print(gf_cmd)
    call(gf_cmd, shell=True)

    # remove duplicate protein sequences that are NOT duplicates on nucleotide levels
    if is_edge:
        print("removing anomalous duplicate protein sequences...")
        final_prot_fn = out_prefix + "_dup_removed_prots_final.fa"
        final_prot_fh = open(final_prot_fn, 'w')
        for name, seq in pyfastx.Fasta(filt_prot_fn, build_index=False):
            if name in to_del_set:
                continue
            elif name in consensus_tbl:
                name_rev = name
                dup_lst = consensus_tbl[name]
                for did in dup_lst:
                    name_rev = name_rev + "|" + did
                write2fa(name_rev, seq, final_prot_fh)
            else:
                write2fa(name, seq, final_prot_fh)

    if is_edge:
        u2d_nt_fn = out_prefix + "_uniq2dup_nt.tsv"
        u2dr_nt_fn = out_prefix + "_uniq2dref_nt.tsv"
        u2d_prot_fn = out_prefix + "_uniq2dup_prot.tsv"
        u2dr_prot_fn = out_prefix + "_uniq2dref_prot.tsv"
        dr2u_nt_fn = out_prefix + "_dref2uniq_nt.tsv"
        dr2u_prot_fn = out_prefix + "_dref2uniq_prot.tsv"
        write_dup_map(nt_dup_tbl, ref_tbl, u2d_nt_fn, u2dr_nt_fn, dr2u_nt_fn, tx_index)
        write_dup_map(prot_dup_tbl, ref_tbl, u2d_prot_fn, u2dr_prot_fn, dr2u_prot_fn, tx_index)
    else:
        u2d_fn = out_prefix + "_uniq2dup.tsv"
        u2dr_fn = out_prefix + "_uniq2dref.tsv"
        dr2u_fn = out_prefix + "_dref2uniq.tsv"
        write_dup_map(nt_dup_tbl, ref_tbl, u2d_fn, u2dr_fn, dr2u_fn, tx_index)


def main(nt_fn, prot_fn, gtf_fn, genome_fn):
    ref_tbl = build_ref_tbl(gtf_fn)
    nt_uniq_lst, nt_dup_tbl = collapse(nt_fn)
    prot_uniq_lst, prot_dup_tbl = collapse(prot_fn)
    is_edge = False
    if len(nt_uniq_lst) == len(prot_uniq_lst):
        print("no edge case to handle. continue...")
    else:
        print("handling a couple of edge cases...")
        is_edge = True
        to_del_set = set()
        consensus_tbl = dict()
        for nt_tid in nt_uniq_lst:
            nt_dup_lst = nt_dup_tbl[nt_tid]
            prot_dup_lst = prot_dup_tbl[nt_tid]
            if len(nt_dup_lst) != len(prot_dup_lst):
                if len(nt_dup_lst) > len(prot_dup_lst):
                    print("this can't happen")
                    sys.exit()
                diff = set(prot_dup_lst) - set(nt_dup_lst)
                consensus_tbl[nt_tid] = diff
                to_del_set |= diff

    filt_gtf_fn = "./dup_removed.gtf"
    tx_index = filter_and_renum_gtf(gtf_fn, filt_gtf_fn, nt_uniq_lst)
    filt_nt_fn = "./dup_removed_txs.fa"
    filt_prot_fn = "./dup_removed_prots.fa"

    gf_cmd = "gffread -w " + filt_nt_fn + " -g " + genome_fn + " " + filt_gtf_fn
    print(gf_cmd)
    call(gf_cmd, shell=True)

    gf_cmd = "gffread -S -y " + filt_prot_fn + " -g " + genome_fn + " " + filt_gtf_fn
    print(gf_cmd)
    call(gf_cmd, shell=True)

    # remove duplicate protein sequences that are NOT duplicates on nucleotide levels
    if is_edge:
        print("removing anomalous duplicate protein sequences...")
        final_prot_fn = "./dup_removed_prots_final.fa"
        final_prot_fh = open(final_prot_fn, 'w')
        for name, seq in pyfastx.Fasta(filt_prot_fn, build_index=False):
            if name in to_del_set:
                continue
            elif name in consensus_tbl:
                name_rev = name
                dup_lst = consensus_tbl[name]
                for did in dup_lst:
                    name_rev = name_rev + "|" + did
                write2fa(name_rev, seq, final_prot_fh)
            else:
                write2fa(name, seq, final_prot_fh)

    if is_edge:
        u2d_nt_fn = "./uniq2dup_nt.tsv"
        u2dr_nt_fn = "./uniq2dref_nt.tsv"
        u2d_prot_fn = "./uniq2dup_prot.tsv"
        u2dr_prot_fn = "./uniq2dref_prot.tsv"
        dr2u_nt_fn = "./dref2uniq_nt.tsv"
        dr2u_prot_fn = "./dref2uniq_prot.tsv"
        write_dup_map(nt_dup_tbl, ref_tbl, u2d_nt_fn, u2dr_nt_fn, dr2u_nt_fn, tx_index)
        write_dup_map(prot_dup_tbl, ref_tbl, u2d_prot_fn, u2dr_prot_fn, dr2u_prot_fn, tx_index)
    else:
        u2d_fn = "./uniq2dup.tsv"
        u2dr_fn = "./uniq2dref.tsv"
        dr2u_fn = "./dref2uniq.tsv"
        write_dup_map(nt_dup_tbl, ref_tbl, u2d_fn, u2dr_fn, dr2u_fn, tx_index)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4])

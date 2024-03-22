#!/usr/bin/env python

import argparse
from time import strftime
import gffutils
import os
import pyfastx
from subprocess import call
from tqdm import tqdm
import sys
from collapse import main_collapse
from commons import load_gene_syn, write_gtf_line

# global vars
orf_d = dict()
junc_d = dict()
gname2gid_d = dict()
gid2tx_d = dict()
ref_cds_d = dict()
orf_cds_d = dict()
tx_cnt = 0
stop_codons = ["TAG", "TAA", "TGA"]


def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-orf", '--orf-gtf', type=str, help="", required=True)
    parser.add_argument("-junc", '--junc-bed', type=str, help="", required=True)
    parser.add_argument("--ref-gtf", type=str, help="", required=False, default=None)
    parser.add_argument("--ref-db", type=str, help="", required=True)
    parser.add_argument("--tmp", type=str, help="", required=True)
    parser.add_argument("-thres", "--threshold", type=float, help="", required=False, default=0.9)
    parser.add_argument("-g", "--genome", type=str, help="", required=True)
    parser.add_argument("-op", "--out-prefix", type=str, help="", required=True)
    parser.add_argument("--remove", help="", required=False, default=False, action='store_true')
    parser.add_argument("--collapse", help="", required=False, default=False, action='store_true')
    parser.add_argument('-syn', type=str, help="", required=False, default=None)
    args = parser.parse_args()
    return args


def load_orfs(args):
    orf_fh = open(args.orf_gtf, 'r')
    global orf_d
    gnames = set()
    for line in orf_fh:
        if line[0] == "#":
            continue
        fields = line.split("\t")
        feature = fields[2]
        if feature == "transcript":
            ctg = fields[0]
            strand = fields[6]
            start = int(fields[3])
            end = int(fields[4])
            infos = fields[8].split(";")
            for info in infos:
                info = info.strip()
                key = info.split(" ")[0]
                val = info.split(" ")[1].replace('"', '')
                if key == "transcript_id":
                    oid = val
                elif key == "associated_gene":
                    gname = val
                elif key == "transcript_biotype":
                    btype = val
                    break
            # for later matching gnames to gids
            gnames.add(gname)
            orf_d[oid] = (ctg, start, end, strand, gname, btype)
    orf_fh.close()
    return gnames


def load_juncs(args):
    junc_fh = open(args.junc_bed, 'r')
    global junc_d
    first = True
    for line in junc_fh:
        if first:
            first = False
            continue
        fields = line.split("\t")
        ctg = fields[0]
        # junction between e1 and e2 (start, end) = last bp of e1 and first bp of e2 - 1
        start = int(fields[1])
        end = int(fields[2])
        jid = fields[3]
        cov = int(fields[4])
        strand = fields[5].strip()
        # don't load if strand in unknown
        if strand == '.':
            continue
        if ctg not in junc_d.keys():
            junc_d[ctg] = [(start, end, jid, cov, strand)]
        else:
            junc_d[ctg].append((start, end, jid, cov, strand))
    junc_fh.close()


def match_gname2gid(args, db, gnames, gsyn_tbl):
    global gname2gid_d
    ref_db = gffutils.FeatureDB(db)
    out_fn = os.path.join(args.tmp, "gene_name2gene_id.tsv")
    out_fh = open(out_fn, 'w')
    out_fh.write("gene_name\tgene_id\n")
    tgt_gene_tbl = dict()
    for gene in ref_db.features_of_type('gene'):
        gname = gene.attributes['gene_name'][0]
        gid = gene.attributes['gene_id'][0]
        tgt_gene_tbl[gname] = gid

    gname_lst = list(gnames)
    for i in tqdm(range(len(gname_lst))):
        gname = gname_lst[i]
        if gname in gsyn_tbl:
            src_gname = gsyn_tbl[gname]
        else:
            src_gname = gname
        if src_gname in tgt_gene_tbl:
            gid = tgt_gene_tbl[src_gname]
            gname2gid_d[gname] = gid
            out_fh.write(gname + "\t" + gid + "\n")
    out_fh.close()


def match_tx2gid(db):
    global gname2gid_d
    global gid2tx_d
    ref_db = gffutils.FeatureDB(db)
    for gname in gname2gid_d.keys():
        gid = gname2gid_d[gname]
        tx_lst = list(ref_db.children(gid, featuretype="transcript", order_by='start'))
        gid2tx_d[gid] = tx_lst


def find_cand_junc(start, end, ctg, strand):
    global junc_d
    juncs = junc_d[ctg]
    cand_juncs = list()
    for jst, jen, jid, jcov, js in juncs:
        # skip unknown strands
        if js != strand:
            continue

        # sanity check
        if jst >= jen:
            print("WARNING: junction start and end positions are not numerically sorted. Exiting...")
            sys.exit()

        if strand == '+':
            if start <= jst <= end <= jen:
                # (sd, sa) = (jst, jen)
                cand_juncs.append((jst, jen, jid, jcov, js))
        else:
            if jst <= (start - 1) <= jen < end:
                # (sa, sd) = (jst, jen)
                cand_juncs.append((jst, jen, jid, jcov, js))
    return cand_juncs


def build_txs(db, oid, ost, oen, ctg, strand, gname, junc, out_fh):
    global gname2gid_d
    global tx_cnt
    global ref_cds_d
    global orf_cds_d
    ref_db = gffutils.FeatureDB(db)

    if gname not in gname2gid_d:
        return False

    gid = gname2gid_d[gname]
    tx_lst = gid2tx_d[gid]
    (jst, jen, jid, jcov, js) = junc
    if strand == '+':
        for tx in tx_lst:
            ref_cds_len = 0
            orf_cds_len = 0
            tid = tx.attributes['transcript_id'][0]
            cds_lst = list(ref_db.children(tid, featuretype='CDS', order_by='start'))
            keep_print = False
            for i in range(len(cds_lst)):
                cds = cds_lst[i]
                ref_cds_len += (cds.end - cds.start + 1)
                if (cds.start - 1) <= jen < cds.end:
                    # sanity check
                    if keep_print:
                        print("Warning: something might have gone wrong here.")
                        sys.exit()
                    tx_cnt += 1
                    orf_tid = "uorft_" + str(tx_cnt)
                    end_wstop = cds_lst[len(cds_lst) - 1].end + 3
                    write_gtf_line(out_fh, ctg, "transcript", ost, end_wstop, strand, gid, oid, orf_tid, tid,
                                   True, (jst, jen))
                    write_gtf_line(out_fh, ctg, "exon", ost, jst, strand, gid, oid, orf_tid, tid, True, (jst, jen))
                    write_gtf_line(out_fh, ctg, "CDS", ost, jst, strand, gid, oid, orf_tid, tid, True, (jst, jen))
                    if jst < ost:
                        print("Error: Invalid feature coordinates")
                        sys.exit()
                    orf_cds_len += (jst - ost + 1)
                    # single exon Txs
                    if i == (len(cds_lst) - 1):
                        write_gtf_line(out_fh, ctg, "exon", jen + 1, cds.end + 3, strand, gid, oid, orf_tid, tid,
                                       True, (jst, jen))
                        write_gtf_line(out_fh, ctg, "CDS", jen + 1, cds.end + 3, strand, gid, oid, orf_tid, tid,
                                       True, (jst, jen))
                    else:
                        write_gtf_line(out_fh, ctg, "exon", jen + 1, cds.end, strand, gid, oid, orf_tid, tid,
                                       True, (jst, jen))
                        write_gtf_line(out_fh, ctg, "CDS", jen + 1, cds.end, strand, gid, oid, orf_tid, tid,
                                       True, (jst, jen))
                    if cds.end < (jen + 1):
                        print("Error: Invalid feature coordinates")
                        sys.exit()
                    orf_cds_len += (cds.end - (jen + 1) + 1)
                    keep_print = True
                elif keep_print:
                    if i == (len(cds_lst) - 1):
                        write_gtf_line(out_fh, ctg, "exon", cds.start, cds.end + 3, strand, gid, oid, orf_tid, tid)
                        write_gtf_line(out_fh, ctg, "CDS", cds.start, cds.end + 3, strand, gid, oid, orf_tid, tid)
                    else:
                        write_gtf_line(out_fh, ctg, "exon", cds.start, cds.end, strand, gid, oid, orf_tid, tid)
                        write_gtf_line(out_fh, ctg, "CDS", cds.start, cds.end, strand, gid, oid, orf_tid, tid)
                    orf_cds_len += (cds.end - cds.start + 1)
            if keep_print:
                if orf_cds_len % 3 == 0:
                    if tid not in ref_cds_d:
                        ref_cds_d[tid] = ref_cds_len
                    if oid not in orf_cds_d.keys():
                        orf_cds_d[oid] = [(orf_tid, orf_cds_len, tid)]
                    else:
                        orf_cds_d[oid].append((orf_tid, orf_cds_len, tid))
    else:
        for tx in tx_lst:
            tid = tx.attributes['transcript_id'][0]
            cds_lst = list(ref_db.children(tid, featuretype='CDS', order_by='start'))
            # last CDS becomes first
            cds_lst.reverse()
            keep_print = False

            ref_cds_len = 0
            orf_cds_len = 0

            for i in range(len(cds_lst)):
                cds = cds_lst[i]
                ref_cds_len += (cds.end - cds.start + 1)
                if cds.start <= jst <= cds.end:
                    if keep_print:
                        print("Warning: something might have gone wrong here.")
                        sys.exit()
                    tx_cnt += 1
                    orf_tid = "uorft_" + str(tx_cnt)
                    start_wstop = cds_lst[len(cds_lst) - 1].start - 3
                    write_gtf_line(out_fh, ctg, "transcript", start_wstop, oen, strand, gid, oid, orf_tid, tid,
                                   True, (jst, jen))
                    write_gtf_line(out_fh, ctg, "exon", jen + 1, oen, strand, gid, oid, orf_tid, tid, True, (jst, jen))
                    write_gtf_line(out_fh, ctg, "CDS", jen + 1, oen, strand, gid, oid, orf_tid, tid, True, (jst, jen))
                    if oen < jen:
                        print("Error: Invalid feature coordinates")
                        sys.exit()
                    orf_cds_len += (oen - (jen + 1) + 1)
                    if i == (len(cds_lst) - 1):
                        write_gtf_line(out_fh, ctg, "exon", cds.start - 3, jst, strand, gid, oid, orf_tid, tid,
                                       True, (jst, jen))
                        write_gtf_line(out_fh, ctg, "CDS", cds.start - 3, jst, strand, gid, oid, orf_tid, tid,
                                       True, (jst, jen))
                    else:
                        write_gtf_line(out_fh, ctg, "exon", cds.start, jst, strand, gid, oid, orf_tid, tid,
                                       True, (jst, jen))
                        write_gtf_line(out_fh, ctg, "CDS", cds.start, jst, strand, gid, oid, orf_tid, tid,
                                       True, (jst, jen))
                    if jst < cds.start:
                        print("Error: Invalid feature coordinates")
                        sys.exit()
                    orf_cds_len += (jst - cds.start + 1)
                    keep_print = True
                elif keep_print:
                    if i == (len(cds_lst) - 1):
                        write_gtf_line(out_fh, ctg, "exon", cds.start - 3, cds.end, strand, gid, oid, orf_tid, tid)
                        write_gtf_line(out_fh, ctg, "CDS", cds.start - 3, cds.end, strand, gid, oid, orf_tid, tid)
                    else:
                        write_gtf_line(out_fh, ctg, "exon", cds.start, cds.end, strand, gid, oid, orf_tid, tid)
                        write_gtf_line(out_fh, ctg, "CDS", cds.start, cds.end, strand, gid, oid, orf_tid, tid)
                    orf_cds_len += (cds.end - cds.start + 1)
            if keep_print:
                if orf_cds_len % 3 == 0:
                    if tid not in ref_cds_d:
                        ref_cds_d[tid] = ref_cds_len
                    if oid not in orf_cds_d.keys():
                        orf_cds_d[oid] = [(orf_tid, orf_cds_len, tid)]
                    else:
                        orf_cds_d[oid].append((orf_tid, orf_cds_len, tid))
    out_fh.flush()
    return True


def write2fa(name, seq, out_fh):
    out_fh.write(">" + name + "\n")
    for i in range(0, len(seq), 60):
        out_fh.write(seq[i:i + 60] + "\n")


def filter_gtf(gtf_fh, out_fh, tx_lst):
    for line in gtf_fh:
        fields = line.split("\t")
        feature = fields[2]
        if feature == "transcript":
            keep_print = False
            infos = fields[8].split(";")
            for info in infos:
                info_clean = info.strip()
                kv_pair = info_clean.split(" ")
                if len(kv_pair) < 2:
                    break
                k = kv_pair[0]
                v = kv_pair[1].replace('"', '')
                if k == "transcript_id":
                    tid = v
                    if tid in tx_lst:
                        out_fh.write(line)
                        keep_print = True
                    break
        elif keep_print:
            out_fh.write(line)


def main():
    args = get_args()

    if not os.path.exists(args.tmp):
        print("specified tmp dir not present. creating it...")
        os.makedirs(args.tmp)

    # load reference db
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Locating reference DB.")
    db = os.path.join(args.tmp, args.ref_db)
    if os.path.isfile(db):
        if args.ref_gtf:
            print(strftime("%Y-%m-%d %H:%M:%S | ") + "Reference DB already present in the tmp dir.\n"
                                                     "You should only provide a reference gtf if you need a new DB "
                                                     "built.\nIf so, please remove this existing DB and re-run.\n"
                                                     "If not, please remove the reference gtf argument.")
            sys.exit()
        print(strftime("%Y-%m-%d %H:%M:%S | ") + "Reference DB found in the tmp dir.")
    else:
        if args.ref_gtf is None:
            print("Reference DB not detected. Please provide a reference gtf to build one.")
        print(strftime("%Y-%m-%d %H:%M:%S | ") + "Reference DB not present. Creating one in the tmp dir.")
        gffutils.create_db(args.ref_gtf, db, id_spec="ID", merge_strategy="create_unique", keep_order=True,
                           disable_infer_transcripts=True, disable_infer_genes=True)

    if args.syn is not None:
        print("loading gene synonyms...")
        gsyn_tbl = load_gene_syn(args)

    # load orfs
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Loading ORFs.")
    gnames = load_orfs(args)
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Finished loading.")

    # matching gene_names to gene_ids
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Fetching gene_name2gene_id mappings.")
    gname2gid_fn = os.path.join(args.tmp, "gene_name2gene_id.tsv")
    if os.path.isfile(gname2gid_fn):
        print(strftime("%Y-%m-%d %H:%M:%S | ") + "Gene names have previously matched to gene ids. Loading...")
        gname2gid_fh = open(gname2gid_fn, 'r')
        first = True
        global gname2gid_d
        for line in gname2gid_fh:
            if first:
                first = False
                continue
            gname = line.split("\t")[0]
            gid = line.split("\t")[1].strip()
            gname2gid_d[gname] = gid
        gname2gid_fh.close()
        print(strftime("%Y-%m-%d %H:%M:%S | ") + "Finished loading.")
    else:
        print(strftime("%Y-%m-%d %H:%M:%S | ") + "Matching gene names to gene ids.")
        match_gname2gid(args, db, gnames, gsyn_tbl)
        print(strftime("%Y-%m-%d %H:%M:%S | ") + "Finished matching.")

    # obtaining txs belonging to each gene_ids
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Fetching gene_id2Tx mappings")
    match_tx2gid(db)

    # load junctions
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Loading junctions.")
    load_juncs(args)

    # iterate through oids to (1) find candidate junctions and (2) build novel Txs
    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Finding candidate splice junctions and building corresponding Txs.")
    global orf_d

    # counters
    tot = 0
    oids = list(orf_d.keys())

    # out files
    tmp_gtf_fn = os.path.join(args.tmp, "tmp.gtf")
    tmp_gtf_fh = open(tmp_gtf_fn, 'w')
    unmatched_fn = os.path.join(args.tmp, "unmatched.txt")
    unmatched_fh = open(unmatched_fn, 'w')
    unmatched_fh.write("orf_id\tgene_name\n")

    for i in tqdm(range(len(oids))):
        oid = oids[i]
        tot += 1
        ctg, start, end, strand, gname, btype = orf_d[oid]
        cand_juncs = find_cand_junc(start, end, ctg, strand)
        for junc in cand_juncs:
            matched = build_txs(db, oid, start, end, ctg, strand, gname, junc, tmp_gtf_fh)
            if not matched:
                unmatched_fh.write(oid + "\t" + gname + "\n")
                break

    tmp_gtf_fh.close()
    unmatched_fh.close()

    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Applying CDS len filter")
    print("Threshold: %f" % args.threshold)

    # apply cds len filter
    global ref_cds_d
    global orf_cds_d
    passed_lst = list()

    # counters
    spliced = 0
    passed = 0

    cds_fn = os.path.join(args.tmp, "cds_len.tsv")
    cds_fh = open(cds_fn, 'w')
    cds_fh.write("orf_id\ttx_id\tref_id\ttx_cds_len\tref_cds_len\tratio\tstatus\n")

    for oid in orf_cds_d.keys():
        spliced += 1
        orf_tx_lst = orf_cds_d[oid]
        qual = False
        for orf_tid, orf_cds_len, ref_tid in orf_tx_lst:
            ref_cds_len = ref_cds_d[ref_tid]
            ratio = orf_cds_len / ref_cds_len
            if ratio >= args.threshold:
                passed_lst.append(orf_tid)
                qual = True
                cds_fh.write(oid + "\t" + orf_tid + "\t" + ref_tid + "\t" + str(orf_cds_len) + "\t" +
                             str(ref_cds_len) + "\t" + str(ratio) + "\tTrue\n")
            else:
                cds_fh.write(oid + "\t" + orf_tid + "\t" + ref_tid + "\t" + str(orf_cds_len) + "\t" +
                             str(ref_cds_len) + "\t" + str(ratio) + "\tFalse\n")
        if qual:
            passed += 1
    cds_fh.close()
    print("%d / %d orfs had >= 1 valid splice junctions" % (spliced, tot))
    print("%d / %d orfs passed the cds len filter" % (passed, spliced))
    print("%d Txs built in total" % len(passed_lst))

    tmp_gtf_fh = open(tmp_gtf_fn, 'r')
    tmp_gtf_filt_fn = os.path.join(args.tmp, "tmp.filt.gtf")
    tmp_gtf_filt_fh = open(tmp_gtf_filt_fn, 'w')
    filter_gtf(tmp_gtf_fh, tmp_gtf_filt_fh, passed_lst)
    tmp_gtf_filt_fh.close()

    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Extracting the transcript / CDS sequences using gffread")
    nt_fa_fn = os.path.join(args.tmp, "txs.fa")
    nt_gf_cmd = "gffread -w " + nt_fa_fn + " -g " + args.genome + " " + tmp_gtf_filt_fn
    print(nt_gf_cmd)
    call(nt_gf_cmd, shell=True)
    prot_fa_fn = os.path.join(args.tmp, "prots.fa")
    prot_gf_cmd = "gffread -S -y " + prot_fa_fn + " -g " + args.genome + " " + tmp_gtf_filt_fn
    print(prot_gf_cmd)
    call(prot_gf_cmd, shell=True)

    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Checking for stop codons at the end of each Tx seqs")
    missing_stop = list()
    passed_lst_2 = list()
    global stop_codons
    for name, seq in pyfastx.Fasta(nt_fa_fn, build_index=False):
        last_codon = seq[-3:].upper()
        tid = name.split(" ")[0].replace('>', '')
        if last_codon not in stop_codons:
            missing_stop.append(tid)
        else:
            passed_lst_2.append(tid)
    # TODO: fix this message
    print("%d / %d Txs were missing stop codons at the end "
          "(i.e., incomplete ORFs annotated in the reference GTF" % (len(missing_stop), len(passed_lst)))
    print("Proceeding with %d Txs with complete ORFs" % (len(passed_lst_2)))

    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Removing Txs with premature stop codons.")
    passed_lst_3 = list()
    orf_completed = list()
    first = True
    for name, seq in pyfastx.Fasta(prot_fa_fn, build_index=False):
        tid = name.split(" ")[0].replace('>', '')
        if first:
            print(tid)
            first = False
        stop_cnt = seq.count("*")

        # # TODO: try out this logic instead
        # if tid in passed_lst_2:
        #     if stop_cnt == 0:
        #         passed_lst_3.append(tid)
        #         passed_lst_3.append(tid)
        # # tid must be in missing_stop
        # else:
        #     if tid not in missing_stop:
        #         print("ERROR: something must have gone in the 2nd filter. Exiting")
        #         sys.exit()
        #     if stop_cnt == 1:
        #         orf_completed.append(tid)

        if stop_cnt > 0:
            if tid in missing_stop:
                if stop_cnt == 1:
                    orf_completed.append(tid)
        else:
            if tid in passed_lst_2:
                passed_lst_3.append(tid)

    print("%d / %d Txs passed the premature stop codon filter" % (len(passed_lst_3), len(passed_lst_2)))
    print("%d / %d Txs w/ missing stop codons has gained exactly one." % (len(orf_completed), len(missing_stop)))

    stop_gained_fn = os.path.join(args.tmp, "stop_gained.txt")
    stop_gained_fh = open(stop_gained_fn, 'w')
    for tid in orf_completed:
        stop_gained_fh.write(tid + "\n")
    stop_gained_fh.close()

    if args.collapse:
        print(strftime("%Y-%m-%d %H:%M:%S | ") + "Further filtering gtf files and obtaining fasta sequences")
        pre_col_gtf_fn = os.path.join(args.tmp, "pre_collapse.gtf")
        pre_col_gtf_fh = open(pre_col_gtf_fn, 'w')
        tmp_gtf_filt_fh = open(tmp_gtf_filt_fn, 'r')
        filter_gtf(tmp_gtf_filt_fh, pre_col_gtf_fh, passed_lst_3)
        pre_col_gtf_fh.close()
        tmp_gtf_filt_fh.close()

        # outputting fasta files
        pre_col_nt_fa_fn = os.path.join(args.tmp, "pre_collapse_txs.fa")
        pre_col_prot_fa_fn = os.path.join(args.tmp, "pre_collapse_prots.fa")
        gf_cmd = "gffread -w " + pre_col_nt_fa_fn + " -g " + args.genome + " " + pre_col_gtf_fn
        print(gf_cmd)
        call(gf_cmd, shell=True)
        gf_cmd = "gffread -S -y " + pre_col_prot_fa_fn + " -g " + args.genome + " " + pre_col_gtf_fn
        print(gf_cmd)
        call(gf_cmd, shell=True)

        print(strftime("%Y-%m-%d %H:%M:%S | ") + "Collapsing duplicate nucleotide / protein sequences")
        # calling the collapse module
        main_collapse(pre_col_nt_fa_fn, pre_col_prot_fa_fn, pre_col_gtf_fn, args.genome, args.out_prefix)
    else:
        print(strftime("%Y-%m-%d %H:%M:%S | ") + "Outputting final gtf and fasta files")
        final_gtf_fn = args.out_prefix + "_filtered.gtf"
        final_gtf_fh = open(final_gtf_fn, 'w')
        tmp_gtf_filt_fh = open(tmp_gtf_filt_fn, 'r')
        filter_gtf(tmp_gtf_filt_fh, final_gtf_fh, passed_lst_3)
        final_gtf_fh.close()
        tmp_gtf_filt_fh.close()

        # outputting fasta files
        final_nt_fa_fn = args.out_prefix + "_filtered_txs.fa"
        final_prot_fa_fn = args.out_prefix + "_filtered_prots.fa"
        gf_cmd = "gffread -w " + final_nt_fa_fn + " -g " + args.genome + " " + final_gtf_fn
        print(gf_cmd)
        call(gf_cmd, shell=True)
        gf_cmd = "gffread -S -y " + final_prot_fa_fn + " -g " + args.genome + " " + final_gtf_fn
        print(gf_cmd)
        call(gf_cmd, shell=True)

    print(strftime("%Y-%m-%d %H:%M:%S | ") + "Finished analysis. Check output files.")


if __name__ == "__main__":
    main()

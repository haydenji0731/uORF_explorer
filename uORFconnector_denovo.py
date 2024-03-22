#!/usr/bin/env python

import argparse
import sys
import os
import gffutils
from commons import write_gtf_line, load_gene_syn
from tqdm import tqdm
import pyfastx
from subprocess import call
import pickle
from collapse import main_collapse


def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-orf", "--orf-gtf", type=str, help="", required=True)
    parser.add_argument("--ref-gtf", type=str, help="", required=False, default=None)
    parser.add_argument("--ref-db", type=str, help="", required=True)
    parser.add_argument("-syn", "--gene-syn", type=str, help="", required=False, default=None)
    parser.add_argument("-map", "--gene-map", type=str, help="", required=True)
    parser.add_argument("-junc", "--junc-bed", type=str, help="", required=False, default=None)
    parser.add_argument("-g", "--genome", type=str, help="", required=True)
    parser.add_argument("-op", "--out-prefix", type=str, help="", required=True)
    parser.add_argument("-c", "--cds-len", type=float, help="", required=False, default=0.9)
    parser.add_argument("-s", "--junc-score", type=float, help="", required=False, default=0.9)
    parser.add_argument("-tmp", "--tmp-dir", type=str, help="", required=False, default="./tmp")
    parser.add_argument("--cont", help="", required=False, default=False, action='store_true')
    parser.add_argument("--split", help="", required=False, default=False, action='store_true')
    parser.add_argument("-b", "--batch-size", type=int, help="", required=False, default=500000)
    parser.add_argument("--collapse", help="", required=False, default=False, action='store_true')
    args = parser.parse_args()
    return args


def load_orfs(args):
    fh = open(args.orf_gtf, 'r')
    orf_tbl = dict()
    coord_tbl = dict()
    gnames = set()
    for line in fh:
        if line[0] == "#":
            continue
        fields = line.split("\t")
        feature = fields[2]
        if feature == "transcript":
            ctg = fields[0]
            strand = fields[6]
            # start = int(fields[3])
            # end = int(fields[4])
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
            coord_tbl[oid] = []
            # orf_tbl[oid] = (ctg, start, end, strand, gname, btype)
            orf_tbl[oid] = (ctg, strand, gname, btype)

            # for later matching gnames to gids
            gnames.add(gname)
        elif feature == "CDS":
            st_pos = int(fields[3])
            en_pos = int(fields[4])
            coord_tbl[oid].append((st_pos, en_pos))

    fh.close()
    return orf_tbl, coord_tbl, gnames


def load_juncs(args):
    junc_fh = open(args.junc_bed, 'r')
    junc_tbl = dict()
    first = True
    for line in tqdm(junc_fh):
        if first:
            first = False
            continue
        fields = line.split("\t")
        ctg = fields[0]
        # junction between e1 and e2 (start, end) = last bp of e1 and first bp of e2 - 1
        start = int(fields[1])
        end = int(fields[2])
        strand = fields[5].strip()
        # don't load if strand in unknown
        if strand == '.':
            continue
        jc = (ctg, start, end, strand)
        if ctg not in junc_tbl:
            junc_tbl[ctg] = {jc}
        else:
            junc_tbl[ctg].add(jc)
    junc_fh.close()
    return junc_tbl


def load_gene_map(args):
    fh = open(args.gene_map, 'r')
    gname2gid_tbl = dict()
    first = True
    for line in fh:
        if first:
            first = False
            continue
        clean_ln = line.strip()
        gname = clean_ln.split("\t")[0]
        gid = clean_ln.split("\t")[1]
        gname2gid_tbl[gname] = gid
    fh.close()
    return gname2gid_tbl


def match_gname2gid(args, db, gnames, gsyn_tbl):
    gname2gid_tbl = dict()
    out_fh = open(args.gene_map, 'w')
    out_fh.write("gene_name\tgene_id\n")
    tgt_gene_tbl = dict()
    for gene in db.features_of_type('gene'):
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
            gname2gid_tbl[gname] = gid
            out_fh.write(gname + "\t" + gid + "\n")
    out_fh.close()
    return gname2gid_tbl


def match_gid2tx(db, gname2gid_tbl):
    gid2tx_tbl = dict()
    ctr = 0
    for gname in tqdm(gname2gid_tbl.keys()):
        gid = gname2gid_tbl[gname]
        tx_lst = list(db.children(gid, featuretype="transcript", order_by='start'))
        mane_tx = None
        for tx in tx_lst:
            if "tag" not in tx.attributes:
                continue
            tag_lst = tx.attributes["tag"]
            if "MANE_Select" in tag_lst:
                mane_tx = tx
        if mane_tx is None:
            ctr += 1
        gid2tx_tbl[gid] = mane_tx

    print("# of gene loci without a MANE transcript: " + str(ctr))
    return gid2tx_tbl


def search_sd(coords, ctg, strand, g_fa):
    spliced = False
    if len(coords) > 1:
        spliced = True

    if not spliced:
        start, end = coords[0]
        s = g_fa.fetch(ctg, (start, end), strand=strand)
        sd_lst = list()
        for i in range(len(s) - 1):
            ss = s[i:i + 2]
            if ss == "GT":
                if strand == "+":
                    sd_lst.append((start + i, start + i + 1))
                elif strand == "-":
                    sd_lst.append((end - i - 1, end - i))
                else:
                    print("ERROR: invalid strand info.")
                    sys.exit()
    else:
        sd_lst = list()
        for i in range(len(coords)):
            coord = coords[i]
            start, end = coord
            # record the annotated splice donor sites
            if strand == '+':
                if i != len(coords) - 1:
                    sd_lst.append((end + 1, end + 2))
            else:
                if i != 0:
                    sd_lst.append((start - 2, start - 1))

            s = g_fa.fetch(ctg, (start, end), strand=strand)
            sd_lst = list()
            for j in range(len(s) - 1):
                ss = s[j:j + 2]
                if ss == "GT":
                    if strand == "+":
                        sd_lst.append((start + j, start + j + 1))
                    elif strand == "-":
                        sd_lst.append((end - j - 1, end - j))
                    else:
                        print("ERROR: invalid strand info.")
                        sys.exit()
    return sd_lst


# functions for completely de novo (sa, sd) search + pairing
def search_sa(tx, g_fa, ctg, strand, db):
    tid = tx.attributes['transcript_id'][0]
    cds_lst = list(db.children(tid, featuretype='CDS', order_by='start'))
    sa_lst = list()
    up_pos = cds_lst[0].start
    down_pos = cds_lst[len(cds_lst) - 1].end
    seq = g_fa.fetch(ctg, (up_pos, down_pos), strand=strand)
    for i in range(len(seq) - 1):
        dint = seq[i:i + 2]
        if dint == "AG":
            if strand == '+':
                sa_lst.append((up_pos + i, up_pos + i + 1))
            elif strand == '-':
                sa_lst.append((down_pos - i - 1, down_pos - i))
            else:
                print("ERROR: invalid strand info.")
                sys.exit()
    return sa_lst


def pair_sa_sd(sd_lst, sa_lst, ctg, strand):
    jc_lst = list()
    if strand == '+':
        for sd in sd_lst:
            for sa in sa_lst:
                if sd[1] < sa[0]:
                    jc = (sd[0] - 2, sa[1] + 1)
                    jc_lst.append((ctg, strand, jc))
    elif strand == '-':
        for sd in sd_lst:
            for sa in sa_lst:
                if sd[0] > sa[1]:
                    jc = (sa[0] - 2, sd[1] + 1)
                    jc_lst.append((ctg, strand, jc))
    return jc_lst


# functions for completely de novo (sa, sd) search + pairing


def extract_jc(ctg, strand, coords, sd, tx, g_fa, db, thres):
    tid = tx.attributes['transcript_id'][0]
    cds_lst = list(db.children(tid, featuretype='CDS', order_by='start'))
    jc_lst = list()
    if tx.strand != strand:
        print("Transcript and u/uo-ORF located on different strands")
        return None

    # first, get the CDS length of the reference tx
    clen = 0
    for cds in cds_lst:
        clen += cds.end - cds.start + 1

    skip_len = 0

    # OLD
    # if strand == '+':
    #     # novel exon length
    #     novel_elen = sd[0] - o_st
    # elif strand == '-':
    #     novel_elen = o_en - sd[1]

    # compute novel exon length
    if len(coords) == 1:
        o_st, o_en = coords[0]
        if strand == '+':
            novel_elen = sd[0] - o_st
        else:
            novel_elen = o_en - sd[1]
    else:
        novel_elen = 0
        if strand == '-':
            coords.reverse()
            for coord in coords:
                st_pos, en_pos = coord
                if st_pos - 1 <= sd[1] <= en_pos:
                    tmp = en_pos - sd[1] + 1
                    assert tmp > 0
                    novel_elen += tmp
                    break
                else:
                    novel_elen += en_pos - st_pos + 1
        else:
            for coord in coords:
                st_pos, en_pos = coord
                if st_pos <= sd[0] <= en_pos + 1:
                    tmp = sd[0] - st_pos + 1
                    assert tmp > 0
                    novel_elen += tmp
                    break
                else:
                    novel_elen += en_pos - st_pos + 1

    if strand == '+':
        for i in range(len(cds_lst)):
            cds = cds_lst[i]
            valid_sa = False

            # check if this position passes the clen filter
            loss_perc = (skip_len - novel_elen) / clen
            if loss_perc > (1 - thres):
                break

            if i == 0:
                s = g_fa.fetch(ctg, (cds.start - 2, cds.start - 1), strand=strand)
                if s == "AG":
                    valid_sa = True
            else:
                valid_sa = True
            if valid_sa:
                up_e_pos = sd[0] - 1
                down_e_pos = cds.start
                if down_e_pos > up_e_pos:
                    # last bp of upstream exon - 1 --> SPLAM coord sys
                    # this way, building tx becomes easier
                    jc = (up_e_pos - 1, down_e_pos, i)
                    if jc in jc_lst:
                        print("this shouldn't happen; like never")
                        sys.exit()
                    jc_lst.append((ctg, strand, jc))

            # update skip length
            skip_len += cds.end - cds.start + 1
    elif strand == '-':
        cds_lst.reverse()
        for i in range(len(cds_lst)):
            cds = cds_lst[i]
            valid_sa = False

            # check if this position passes the clen filter
            loss_perc = (skip_len - novel_elen) / clen
            if loss_perc > (1 - thres):
                break

            if i == 0:
                s = g_fa.fetch(ctg, (cds.end + 1, cds.end + 2), strand=strand)
                if s == "AG":
                    valid_sa = True
            else:
                valid_sa = True
            if valid_sa:
                up_e_pos = cds.end
                down_e_pos = sd[1] + 1
                if down_e_pos > up_e_pos:
                    # this way, building tx becomes easier
                    jc = (up_e_pos - 1, down_e_pos, i)
                    if jc in jc_lst:
                        sys.exit()
                    jc_lst.append((ctg, strand, jc))

            # update skip length
            skip_len += cds.end - cds.start + 1
    else:
        print("ERROR: invalid strand info.")
        sys.exit()
    return jc_lst


def build_tx(ctg, tx, coords, jc, strand, db, ctr, out_fh, gid, oid):
    ref_clen = 0
    denovo_clen = 0
    up_pos = jc[0]
    down_pos = jc[1]
    cds_index = jc[2]
    tid = tx.attributes['transcript_id'][0]
    cds_lst = list(db.children(tid, featuretype='CDS', order_by='start'))
    orf_tid = "uorft_denovo_" + str(ctr)
    keep_print = False
    if strand == '+':
        st_pos, _ = coords[0]
        en_pos = cds_lst[len(cds_lst) - 1].end + 3
        write_gtf_line(out_fh, ctg, "transcript", st_pos, en_pos, strand, gid, oid, orf_tid, tid,
                       True, (up_pos + 1, down_pos - 1))
        for coord in coords:
            o_st, o_en = coord
            if o_st <= up_pos + 1 <= o_en:
                # record splice junction in TieBrush coordinate system
                write_gtf_line(out_fh, ctg, "exon", o_st, (up_pos + 1), strand, gid, oid, orf_tid, tid,
                               True, (up_pos + 1, down_pos - 1))
                write_gtf_line(out_fh, ctg, "CDS", o_st, (up_pos + 1), strand, gid, oid, orf_tid, tid,
                               True, (up_pos + 1, down_pos - 1))
                denovo_clen += ((up_pos + 1) - o_st + 1)
                break
            else:
                write_gtf_line(out_fh, ctg, "exon", o_st, o_en, strand, gid, oid, orf_tid, tid)
                write_gtf_line(out_fh, ctg, "CDS", o_st, o_en, strand, gid, oid, orf_tid, tid)
                denovo_clen += (o_en - o_st + 1)

        for i in range(len(cds_lst)):
            cds = cds_lst[i]
            ref_clen += (cds.end - cds.start + 1)
            if i == cds_index:
                if i == (len(cds_lst) - 1):
                    # fixed the splice junction annotations; has NO impact on the transcript sequences
                    write_gtf_line(out_fh, ctg, "exon", down_pos, en_pos, strand, gid, oid, orf_tid, tid,
                                   True, (up_pos + 1, down_pos - 1))
                    write_gtf_line(out_fh, ctg, "CDS", down_pos, en_pos, strand, gid, oid, orf_tid, tid,
                                   True, (up_pos + 1, down_pos - 1))
                else:
                    write_gtf_line(out_fh, ctg, "exon", down_pos, cds.end, strand, gid, oid, orf_tid, tid,
                                   True, (up_pos + 1, down_pos - 1))
                    write_gtf_line(out_fh, ctg, "CDS", down_pos, cds.end, strand, gid, oid, orf_tid, tid,
                                   True, (up_pos + 1, down_pos - 1))
                if cds.end < down_pos:
                    print("ERROR: Invalid feature coordinates")
                    sys.exit()
                # exclude the stop codon in cds len
                denovo_clen += (cds.end - down_pos + 1)
                keep_print = True
            elif keep_print:
                if i == (len(cds_lst) - 1):
                    write_gtf_line(out_fh, ctg, "exon", cds.start, cds.end + 3, strand, gid, oid, orf_tid, tid)
                    write_gtf_line(out_fh, ctg, "CDS", cds.start, cds.end + 3, strand, gid, oid, orf_tid, tid)
                else:
                    write_gtf_line(out_fh, ctg, "exon", cds.start, cds.end, strand, gid, oid, orf_tid, tid)
                    write_gtf_line(out_fh, ctg, "CDS", cds.start, cds.end, strand, gid, oid, orf_tid, tid)
                denovo_clen += (cds.end - cds.start + 1)
    else:
        cds_lst.reverse()
        coords.reverse()
        _, en_pos = coords[0]
        st_pos = cds_lst[len(cds_lst) - 1].start - 3
        write_gtf_line(out_fh, ctg, "transcript", st_pos, en_pos, strand, gid, oid, orf_tid, tid,
                       True, (up_pos + 1, down_pos - 1))
        for coord in coords:
            o_st, o_en = coord
            if o_st <= down_pos <= o_en:
                write_gtf_line(out_fh, ctg, "exon", down_pos, o_en, strand, gid, oid, orf_tid, tid,
                               True, (up_pos + 1, down_pos - 1))
                write_gtf_line(out_fh, ctg, "CDS", down_pos, o_en, strand, gid, oid, orf_tid, tid,
                               True, (up_pos + 1, down_pos - 1))
                denovo_clen += (o_en - down_pos + 1)
                break
            else:
                write_gtf_line(out_fh, ctg, "exon", o_st, o_en, strand, gid, oid, orf_tid, tid)
                write_gtf_line(out_fh, ctg, "CDS", o_st, o_en, strand, gid, oid, orf_tid, tid)
                denovo_clen += (o_en - o_st + 1)

        for i in range(len(cds_lst)):
            cds = cds_lst[i]
            ref_clen += (cds.end - cds.start + 1)
            if i == cds_index:
                if i == (len(cds_lst) - 1):
                    write_gtf_line(out_fh, ctg, "exon", st_pos, up_pos + 1, strand, gid, oid, orf_tid, tid,
                                   True, (up_pos + 1, down_pos - 1))
                    write_gtf_line(out_fh, ctg, "CDS", st_pos, up_pos + 1, strand, gid, oid, orf_tid, tid,
                                   True, (up_pos + 1, down_pos - 1))
                else:
                    write_gtf_line(out_fh, ctg, "exon", cds.start, up_pos + 1, strand, gid, oid, orf_tid, tid,
                                   True, (up_pos + 1, down_pos - 1))
                    write_gtf_line(out_fh, ctg, "CDS", cds.start, up_pos + 1, strand, gid, oid, orf_tid, tid,
                                   True, (up_pos + 1, down_pos - 1))
                if (up_pos + 1) < cds.start:
                    print("ERROR: Invalid feature coordinates")
                    sys.exit()
                denovo_clen += ((up_pos + 1) - cds.start + 1)
                keep_print = True
            elif keep_print:
                if i == (len(cds_lst) - 1):
                    write_gtf_line(out_fh, ctg, "exon", cds.start - 3, cds.end, strand, gid, oid, orf_tid, tid)
                    write_gtf_line(out_fh, ctg, "CDS", cds.start - 3, cds.end, strand, gid, oid, orf_tid, tid)
                else:
                    write_gtf_line(out_fh, ctg, "exon", cds.start, cds.end, strand, gid, oid, orf_tid, tid)
                    write_gtf_line(out_fh, ctg, "CDS", cds.start, cds.end, strand, gid, oid, orf_tid, tid)
                denovo_clen += (cds.end - cds.start + 1)
    # check for divisibility by 3
    if denovo_clen % 3 == 0:
        return ref_clen, denovo_clen
    return False


def filter_gtf(gtf_fh, out_fh, tx_set):
    for line in tqdm(gtf_fh):
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
                    if tid not in tx_set:
                        out_fh.write(line)
                        keep_print = True
                    break
        elif keep_print:
            out_fh.write(line)


def get_clen(db):
    clen_tbl = dict()
    tx_lst = list(db.features_of_type('transcript', order_by='start'))
    for tx in tqdm(tx_lst):
        tid = tx.attributes['transcript_id'][0]
        cds_lst = list(db.children(tid, featuretype='CDS', order_by='start'))
        clen = 0
        for cds in cds_lst:
            clen += (cds.end - cds.start + 1)
        clen_tbl[tid] = clen
    return clen_tbl


def main():
    args = get_args()

    if not os.path.exists(args.tmp_dir):
        os.makedirs(args.tmp_dir)

    if os.path.isfile(args.ref_db):
        if args.ref_gtf:
            print("WARNING: Reference db detected. If you want to rebuild it, please remove the --ref-db argument.")
    else:
        print("Creating reference db...")
        gffutils.create_db(args.ref_gtf, args.ref_db, id_spec="ID", merge_strategy="create_unique", keep_order=True,
                           disable_infer_transcripts=True, disable_infer_genes=True)
    db = gffutils.FeatureDB(args.ref_db)
    load_syn = False
    if args.gene_syn:
        print("loading gene synonyms...")
        gsyn_tbl = load_gene_syn(args.gene_syn)
        load_syn = True
    print("loading orfs...")
    orf_tbl, coord_tbl, gnames = load_orfs(args)
    if os.path.isfile(args.gene_map):
        print("loading gene name-to-gene id mappings...")
        gname2gid_tbl = load_gene_map(args)
    else:
        print("matching gene names to reference gene ids...")
        if load_syn:
            gname2gid_tbl = match_gname2gid(args, db, gnames, gsyn_tbl)
        else:
            gname2gid_tbl = match_gname2gid(args, db, gnames, dict())

    print("match gene ids to transcripts...")
    gid2tx_tbl = match_gid2tx(db, gname2gid_tbl)

    if args.junc_bed:
        print("Loading known junctions...")
        filt_junc_tbl = load_juncs(args)
    else:
        filt_junc_tbl = None
        print("No junctions.bed provided. Keeping all de novo splice sites...")

    print("loading the reference genome...")
    g_fa = pyfastx.Fasta(args.genome)
    oids = list(orf_tbl.keys())

    if not args.cont:
        print("compiling a list of de novo junctions")
        jc_tbl = dict()
        tot = 0

        dnv_jc_set = set()

        sd_tbl = dict()
        cnt_tbl = dict()
        for i in tqdm(range(len(oids))):
            tot += 1
            oid = oids[i]
            ctg, strand, gname, btype = orf_tbl[oid]
            coords = coord_tbl[oid]

            sd_lst = search_sd(coords, ctg, strand, g_fa)
            sd_tbl[oid] = len(sd_lst)
            gid = gname2gid_tbl[gname]
            tx = gid2tx_tbl[gid]
            cnt_tbl[oid] = 0

            # no MANE tx at this gene locus
            if tx is None:
                continue

            # assert len(sd_lst) == len(set(sd_lst))
            for sd in sd_lst:
                jc_lst = extract_jc(ctg, strand, coords, sd, tx, g_fa, db, args.cds_len)
                assert len(jc_lst) == len(set(jc_lst))
                cnt_tbl[oid] += len(jc_lst)
                dnv_jc_lst = list()
                if filt_junc_tbl:
                    for ctg, strand, jc in jc_lst:
                        # revert it back to TieBrush coordinate system
                        jc_obj = (ctg, jc[0] + 1, jc[1] - 1, strand)
                        if jc_obj not in filt_junc_tbl[ctg]:
                            # still in SPLAM coordinate system
                            dnv_jc_lst.append((ctg, strand, jc))
                else:
                    dnv_jc_lst = jc_lst

                if oid in jc_tbl:
                    jc_tbl[oid].append((tx, dnv_jc_lst))
                else:
                    jc_tbl[oid] = [(tx, dnv_jc_lst)]

                for ctg, strand, jc in dnv_jc_lst:
                    tmp = (ctg, strand, (jc[0], jc[1]))
                    dnv_jc_set.add(tmp)

        print("# of novel de novo junctions: " + str(len(dnv_jc_set)))

        # TODO: useless outfile; delete
        cnt_fn = args.out_prefix + ".denovo.junc.tsv"
        cnt_fh = open(cnt_fn, 'w')
        for oid in cnt_tbl:
            cnt_fh.write(oid + "\t" + str(cnt_tbl[oid]) + "\n")
        cnt_fh.close()

        if len(dnv_jc_set) > 0:
            novel_fn = args.out_prefix + ".denovo.junc.bed"
            novel_fh = open(novel_fn, 'w')
            ctr = 0
            for ctg, strand, jc in dnv_jc_set:
                ctr += 1
                pad = 8 - len(str(ctr))
                junc_name = "JUNC" + ('0' * pad) + str(ctr)
                novel_fh.write(ctg + "\t" + str(jc[0]) + "\t" + str(jc[1]) + "\t" + junc_name + "\t0\t" + strand + "\n")
            novel_fh.close()

        junc_index = args.out_prefix + ".junc.index"
        with open(junc_index, "wb") as f:
            pickle.dump(jc_tbl, f)
    else:
        junc_index = args.out_prefix + ".junc.index"
        if os.path.isfile(junc_index):
            with open(junc_index, "rb") as f:
                jc_tbl = pickle.load(f)
        else:
            print("ARGUMENT ERROR: no preliminary junction index provided.")
            sys.exit()
        print("loading preliminary junctions")
        junc_fn = args.out_prefix + ".denovo.junc.bed"
        junc_fh = open(junc_fn, 'r')
        name_tbl = dict()
        glob_jc_set = set()
        for line in tqdm(junc_fh):
            clean_ln = line.strip()
            fields = clean_ln.split("\t")
            ctg = fields[0]
            up_pos = int(fields[1])
            down_pos = int(fields[2])
            name = fields[3]
            strand = fields[4]
            name_tbl[name] = (ctg, strand, (up_pos, down_pos))
            tmp = (ctg, strand, (up_pos, down_pos))
            glob_jc_set.add(tmp)

    # evaluate junctions w/ SPLAM
    print("scoring with SPLAM...")
    splam_fn = args.tmp_dir + "/junction_score.bed"
    # TODO: fix this later...

    skip = False
    if args.cont:
        if os.path.isfile(splam_fn):
            skip = True
        else:
            print("SPLAM output bed not detected... running SPLAM to score denovo junctions")

    if not skip:
        if args.split:
            print("splitting input .bed to fit into memory")
            tot_jc = len(glob_jc_set)
            bsize = args.batch_size
            num_it = int(tot_jc / bsize)
            i = 0
            merged_fn = args.tmp_dir + "/merged.bed"
            while i < num_it:
                ext_cmd = "head -n " + str((i + 1) * bsize) + " " + junc_fn + " | tail -n " + str(bsize) + " > " + \
                          args.tmp_dir + "/tmp.bed"
                print(ext_cmd)
                call(ext_cmd, shell=True)
                splam_cmd = "splam score -G " + args.genome + " -m " + "./splam/model/splam_script.pt -o " + \
                            args.tmp_dir + " " + args.tmp_dir + "/tmp.bed"
                print(splam_cmd)
                call(splam_cmd, shell=True)
                pipe_cmd = "cat " + splam_fn + " >> " + merged_fn
                print(pipe_cmd)
                call(pipe_cmd, shell=True)
                i += 1
            remaining = tot_jc % bsize
            if remaining > 0:
                ext_cmd = "tail -n " + str(remaining) + " > " + args.tmp_dir + "/tmp.bed"
                print(ext_cmd)
                call(ext_cmd, shell=True)
                splam_cmd = "splam score -G " + args.genome + " -m " + "./splam/model/splam_script.pt -o " + \
                            args.tmp_dir + " " + args.tmp_dir + "/tmp.bed"
                print(splam_cmd)
                call(splam_cmd, shell=True)
                pipe_cmd = "cat " + splam_fn + " >> " + merged_fn
                print(pipe_cmd)
                call(pipe_cmd, shell=True)

            rm_cmd = "rm " + args.tmp_dir + "/tmp.bed"
            print(rm_cmd)
            call(rm_cmd, shell=True)
            mv_cmd = "mv " + merged_fn + " " + splam_fn
            print(mv_cmd)
            call(mv_cmd, shell=True)
        else:
            splam_cmd = "splam score -G " + args.genome + " -m " + "./splam/model/splam_script.pt -o " + \
                        args.tmp_dir + " " + junc_fn
            print(splam_cmd)
            call(splam_cmd, shell=True)

    print("filtering out de novo junctions below the score threshold...")
    print("splice score threshold: " + str(args.junc_score))

    junc_filt_fn = args.out_prefix + ".denovo.junc.kept.bed"
    skip = False
    if args.cont:
        if os.path.isfile(junc_filt_fn):
            print("Loading junctions with passing scores...")
            junc_filt_fh = open(junc_filt_fn, 'r')
            kept_junc_set = set()
            for line in junc_filt_fh:
                clean_ln = line.strip()
                fields = clean_ln.split("\t")
                ctg = fields[0]
                strand = fields[5]
                up_pos = int(fields[1])
                down_pos = int(fields[2])
                kept_junc_set.add((ctg, strand, (up_pos, down_pos)))
            junc_filt_fh.close()
            skip = True
        else:
            print("There must be a BED file containing junctions with passing scores. "
                  "Creating it based on previous SPLAM output.")

    # TODO: check if this is correct
    if not skip:
        fh = open(splam_fn, 'r')
        junc_filt_fh = open(junc_filt_fn, 'w')
        ctr = 0
        kept = 0
        kept_junc_set = set()
        for line in tqdm(fh):
            ctr += 1
            clean_ln = line.strip()
            fields = clean_ln.split("\t")

            ctg = fields[0]
            up_pos = int(fields[1])
            down_pos = int(fields[2])
            strand = fields[5]

            up_score = float(fields[6])
            down_score = float(fields[7])
            avg_score = up_score + down_score / 2
            if avg_score >= args.junc_score:
                kept += 1
                junc_filt_fh.write(line)
                kept_junc_set.add((ctg, strand, (up_pos, down_pos)))
        junc_filt_fh.close()
        fh.close()
        print("%d / %d junctions passed the score filter" % (kept, ctr))
        assert kept == len(kept_junc_set)

    print("Building transcripts with passing junctions")
    init_fn = args.out_prefix + ".init.gtf"
    init_fh = open(init_fn, 'w')
    ctr = 0
    to_del_lst = list()

    for oid in tqdm(jc_tbl):
        # ctg, start, end, strand, gname, btype = orf_tbl[oid]
        ctg, strand, gname, btype = orf_tbl[oid]
        coords = coord_tbl[oid]
        gid = gname2gid_tbl[gname]
        item = jc_tbl[oid]
        for tx, jc_lst in item:

            for ctg, strand, jc in jc_lst:
                # only keep junctions with qualifying scores
                # jc_tbl still contains junctions written in SPLAM coordinate system
                tmp = (ctg, strand, (jc[0], jc[1]))
                if tmp not in kept_junc_set:
                    continue
                ctr += 1
                tid = "uorft_denovo_" + str(ctr)
                res = build_tx(ctg, tx, coords, jc, strand, db, ctr, init_fh, gid, oid)
                if not res:
                    to_del_lst.append(tid)

    print("%d Txs were constructed" % ctr)
    print("%d / %d Txs passed the frameshift filter" % (ctr - len(to_del_lst), ctr))

    print("filtering out the initial gtf file...")
    init_fh = open(init_fn, 'r')
    filt_fn = args.out_prefix + ".filt.gtf"
    filt_fh = open(filt_fn, 'w')
    filter_gtf(init_fh, filt_fh, set(to_del_lst))
    filt_fh.close()
    init_fh.close()

    print("Extracting the transcript / protein sequences... Please make sure that gffread is in your PATH variable.")
    nt_fa_fn = args.out_prefix + ".filt.nt.fa"
    nt_gf_cmd = "gffread -w " + nt_fa_fn + " -g " + args.genome + " " + filt_fn
    print(nt_gf_cmd)
    call(nt_gf_cmd, shell=True)
    pt_fa_fn = args.out_prefix + ".filt.pt.fa"
    pt_gf_cmd = "gffread -S -y " + pt_fa_fn + " -g " + args.genome + " " + filt_fn
    print(pt_gf_cmd)
    call(pt_gf_cmd, shell=True)

    print("Checking for stop codons at the end...")
    stop_codons = ["TAG", "TAA", "TGA"]
    missing_stop_lst = list()
    tot = 0
    for name, seq in pyfastx.Fasta(nt_fa_fn, build_index=False):
        tot += 1
        last_codon = seq[-3:].upper()
        tid = name.split(" ")[0].replace('>', '')
        if last_codon not in stop_codons:
            missing_stop_lst.append(tid)
    print("Checking for premature stop codons...")
    ptc_lst = list()
    orf_completed = 0
    for name, seq in pyfastx.Fasta(pt_fa_fn, build_index=False):
        tid = name.split(" ")[0].replace('>', '')
        stop_cnt = seq.count("*")
        if stop_cnt > 0:
            ptc_lst.append(tid)
            if tid in missing_stop_lst and stop_cnt == 1:
                orf_completed += 1
    to_del_lst = list(set(missing_stop_lst) | set(ptc_lst))
    print("%d / %d Txs passed the premature stop codon + missing stop codon filter" % (tot - len(to_del_lst), tot))

    if args.collapse:
        print("Collapsing duplicate protein sequences...")
        pcol_gtf_fn = args.out_prefix + ".pcol.gtf"
        pcol_gtf_fh = open(pcol_gtf_fn, 'w')
        fh = open(filt_fn, 'r')
        filter_gtf(fh, pcol_gtf_fh, to_del_lst)
        pcol_gtf_fh.close()
        fh.close()

        # create fasta files
        pcol_nt_fa_fn = args.out_prefix + ".final.nt.fa"
        pcol_pt_fa_fn = args.out_prefix + ".final.pt.fa"
        gf_cmd = "gffread -w " + pcol_nt_fa_fn + " -g " + args.genome + " " + pcol_gtf_fn
        print(gf_cmd)
        call(gf_cmd, shell=True)
        gf_cmd = "gffread -S -y " + pcol_pt_fa_fn + " -g " + args.genome + " " + pcol_gtf_fn
        print(gf_cmd)
        call(gf_cmd, shell=True)
        main_collapse(pcol_nt_fa_fn, pcol_pt_fa_fn, pcol_gtf_fn, args.genome, args.out_prefix)
    else:
        print("Outputting final gtf and fasta files")
        final_gtf_fn = args.out_prefix + ".final.gtf"
        final_gtf_fh = open(final_gtf_fn, 'w')
        fh = open(filt_fn, 'r')
        filter_gtf(fh, final_gtf_fh, to_del_lst)
        final_gtf_fh.close()
        fh.close()

        # create fasta files
        final_nt_fa_fn = args.out_prefix + ".final.nt.fa"
        final_pt_fa_fn = args.out_prefix + ".final.pt.fa"
        gf_cmd = "gffread -w " + final_nt_fa_fn + " -g " + args.genome + " " + final_gtf_fn
        print(gf_cmd)
        call(gf_cmd, shell=True)
        gf_cmd = "gffread -S -y " + final_pt_fa_fn + " -g " + args.genome + " " + final_gtf_fn
        print(gf_cmd)
        call(gf_cmd, shell=True)


if __name__ == "__main__":
    main()

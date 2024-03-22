#!/usr/bin/env python

import argparse
import os
import gffutils
from subprocess import call
from tqdm import tqdm
import sys
import pyfastx
import re
import mappy as mp
import pickle
import orfipy_core
from commons import load_gene_syn, translate

USE_EQX = 0x4000000
synonymous_mut = dict()
no_gname = set()


def build_target_tbls(db_fn):
    db = gffutils.FeatureDB(db_fn)
    gene_lst = list(db.features_of_type('gene'))
    gene_tbl = dict()
    tx_tbl = dict()
    gene_set = set()
    global no_gname

    for i in tqdm(range(len(gene_lst))):
        gene = gene_lst[i]
        # TODO: fix this better
        try:
            gene_type = gene.attributes['gene_biotype'][0]
            # gene_type = gene.attributes['gene_biotype'][0] == "pseudogene"
        except:
            gene_type = gene.attributes['gene_type'][0]
            # gene_type = gene.attributes['gene_type'][0] == "pseudogene"
        if gene_type == "pseudogene":
            gid = gene.id
            gene_tbl[gid] = list()
            if 'gene_name' in gene.attributes:
                gname = gene.attributes['gene_name'][0]
                gene_tbl[gid].append(gname)
                gene_set.add(gname)
            elif 'gene' in gene.attributes:
                gname = tx.attributes['gene'][0]
                gene_tbl[gid].append(gname)
                gene_set.add(gname)
            continue
        tx_lst = list(db.children(gene, featuretype="transcript", order_by='start'))

        if len(tx_lst) == 0:
            continue

        tx = tx_lst[0]
        gid = tx.attributes['Parent'][0]
        gene_tbl[gid] = set()

        try:
            gid = gene.id
        except:
            print("investigate this edge case more closely!")

        # TODO: revisit
        if 'gene_synonym' in gene.attributes:
            gname = gene.attributes['gene_synonym'][0]
            gene_tbl[gid].add(gname)
            gene_set.add(gname)

        if 'source_gene_common_name' in tx.attributes:
            gname = tx.attributes['source_gene_common_name'][0]
            gene_tbl[gid].add(gname)
            gene_set.add(gname)
        elif 'gene_name' in tx.attributes:
            gname = tx.attributes['gene_name'][0]
            gene_tbl[gid].add(gname)
            gene_set.add(gname)
        elif 'gene' in tx.attributes:
            gname = tx.attributes['gene'][0]
            gene_tbl[gid].add(gname)
            gene_set.add(gname)

        if len(gene_tbl[gid]) == 0:
            print("WARNING: no gene name detected for gene: " + gid)
            print("adding the gene id as is...")
            gene_tbl[gid].add(gid)
            no_gname.add(gid)

        for tx in tx_lst:
            tid = tx.attributes['ID'][0]
            tx_tbl[tid] = (gid, tx)

    return gene_tbl, tx_tbl, gene_set


def get_cds(args):
    fa = pyfastx.Fasta(args.target_fasta)
    pat = "CDS"
    cds_tbl = dict()
    for seq in fa:
        res = re.search(pat, seq.description)
        if res is not None:
            cds_info = seq.description.split("=")[1]
            # convert to 0-based, half-open-half-closed coordinate system
            cds_st = int(cds_info.split("-")[0]) - 1
            cds_en = int(cds_info.split("-")[1])
            cds_tbl[seq.name] = (cds_st, cds_en)
    return cds_tbl


def load_orfs(args):
    fh = open(args.orf_gtf, 'r')
    orf_tbl = dict()
    for line in fh:
        fields = line.split("\t")
        feature = fields[2]
        if feature == "transcript":
            ctg = fields[0]
            strand = fields[6]
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
            orf_tbl[oid] = (ctg, strand, gname, btype)
    return orf_tbl


def map_to_gene_ome(args, orf_tbl):
    global synonymous_mut
    full_match_tbl = dict()
    partial_match_tbl = dict()
    if args.genome_index is None:
        index_fn = os.path.join(args.tmp_dir, "target_genome.mmi")
        print("mm2 genome index argument not provided. will be looking for it in tmp location: %s" % index_fn)
        if not os.path.exists(args.tmp_dir):
            print("tmp directory wasn't detected. creating it...")
            os.makedirs(args.tmp_dir)
    else:
        index_fn = args.genome_index

    if not os.path.isfile(index_fn):
        print("mm2 genome index not found. creating it...")
        a = mp.Aligner(fn_idx_in=args.genome, fn_idx_out=index_fn, preset="sr", n_threads=args.threads,
                       k=args.kmer_size)
    else:
        a = mp.Aligner(fn_idx_in=index_fn, preset="sr", n_threads=args.threads, k=args.kmer_size)

    for name, seq, qual in mp.fastx_read(args.orf_fasta):
        perfect = str(len(seq)) + "="
        orf_strand = orf_tbl[name][1]
        for hit in a.map(seq, cs=True, extra_flags=USE_EQX):
            if hit.strand == -1 and orf_strand != '-':
                continue
            elif hit.strand == 1 and orf_strand != '+':
                continue
            if hit.cigar_str == perfect:
                if name not in full_match_tbl:
                    full_match_tbl[name] = [hit]
                else:
                    full_match_tbl[name].append(hit)
            else:
                r_seq = a.seq(hit.ctg, start=hit.r_st, end=hit.r_en).upper()
                if hit.strand == -1:
                    r_seq = mp.revcomp(r_seq)
                q_aa = translate(seq.upper())
                # checking for syn mutations
                for start, stop, strand, desc in orfipy_core.orfs(r_seq, minlen=0, maxlen=1000, starts=['ATG']):
                    # since it's reverse complemented, all orfs must be in '+' direction
                    if strand == '-':
                        continue
                    r_orf_nt = r_seq[start:stop + 3]
                    # never enters this part of the code...
                    if len(r_orf_nt) != len(seq):
                        break
                    r_orf_aa = translate(r_orf_nt)
                    ctr = 0
                    for i in range(len(q_aa)):
                        if r_orf_aa[i] == q_aa[i]:
                            ctr += 1
                        else:
                            break
                    if ctr == len(q_aa):
                        if name in synonymous_mut:
                            synonymous_mut[name].append(hit)
                        else:
                            synonymous_mut[name] = [hit]
                        if name not in full_match_tbl:
                            full_match_tbl[name] = [hit]
                        else:
                            full_match_tbl[name].append(hit)
                # storing the partial match
                if name not in partial_match_tbl:
                    partial_match_tbl[name] = [hit]
                else:
                    partial_match_tbl[name].append(hit)
    return full_match_tbl, partial_match_tbl


def map_to_tx_ome(args):
    global synonymous_mut
    full_match_tbl = dict()
    partial_match_tbl = dict()
    if args.transcript_index is None:
        index_fn = os.path.join(args.tmp_dir, "target_transcriptome.mmi")
        print("mm2 transcriptome index argument not provided. will be looking for it in tmp location: %s" % index_fn)
        if not os.path.exists(args.tmp_dir):
            print("tmp directory wasn't detected. creating it...")
            os.makedirs(args.tmp_dir)
    else:
        index_fn = args.transcript_index

    if not os.path.isfile(index_fn):
        print("mm2 transcriptome index not found. creating it...")
        a = mp.Aligner(fn_idx_in=args.target_fasta, fn_idx_out=index_fn, preset="sr", n_threads=args.threads,
                       k=args.kmer_size)
    else:
        a = mp.Aligner(fn_idx_in=index_fn, preset="sr", n_threads=args.threads, k=args.kmer_size)

    for name, seq, qual in mp.fastx_read(args.orf_fasta):
        perfect = str(len(seq)) + "="
        for hit in a.map(seq, cs=True, extra_flags=USE_EQX):
            if hit.strand != 1:
                continue
            if hit.cigar_str == perfect:
                if name not in full_match_tbl:
                    full_match_tbl[name] = [hit]
                else:
                    full_match_tbl[name].append(hit)
            else:
                r_seq = a.seq(hit.ctg, start=hit.r_st, end=hit.r_en).upper()
                q_aa = translate(seq.upper())
                # checking for synonymous mutation, preserving underlying amino acid seq
                for start, stop, strand, desc in orfipy_core.orfs(r_seq, minlen=0, maxlen=1000, starts=['ATG']):
                    if strand == '-':
                        continue
                    r_orf_nt = r_seq[start:stop + 3]
                    if len(r_orf_nt) != len(seq):
                        break
                    r_orf_aa = translate(r_orf_nt)
                    ctr = 0
                    for i in range(len(q_aa)):
                        if r_orf_aa[i] == q_aa[i]:
                            ctr += 1
                        else:
                            break
                    if ctr == len(q_aa):
                        if name in synonymous_mut:
                            synonymous_mut[name].append(hit)
                        else:
                            synonymous_mut[name] = [hit]
                        if name not in full_match_tbl:
                            full_match_tbl[name] = [hit]
                        else:
                            full_match_tbl[name].append(hit)
                # storing the partial match
                if name not in partial_match_tbl:
                    partial_match_tbl[name] = [hit]
                else:
                    partial_match_tbl[name].append(hit)
    return full_match_tbl, partial_match_tbl


def check_missing_gene(orf_tbl, gene_set, gsyn_tbl, gene_tbl=None, db_fn=None):
    if db_fn is not None:
        if gene_tbl is None:
            print("ARGUMENT ERROR: please provide the gene table for check_missing_gene() when --recheck flag is ON")
            sys.exit()
        print("re-matching missing gene names by iterating through transcripts as well")
        db = gffutils.FeatureDB(db_fn)
        tx_lst = list(db.features_of_type('transcript'))
        # TODO: is this right? :0
        for i in tqdm(range(len(tx_lst))):
            tx = tx_lst[i]
            gid = tx.attributes["Parent"][0]
            try:
                gname = tx.attributes["source_gene_common_name"][0]
            except:
                continue
            if gid not in gene_tbl:
                print("no gene entry for this transcript's parent?")
                sys.exit()
            gene_tbl[gid].add(gname)
            gene_set.add(gname)

    missing_gene_set = set()
    for oid in orf_tbl:
        ctg, strand, gname, btype = orf_tbl[oid]
        if gname in gene_set:
            continue
        elif gname in gsyn_tbl:
            continue
        else:
            missing_gene_set.add(gname)
    return missing_gene_set


# written for sanity check
def compare_orf_to_gene_strand(orf_tbl, gene_tbl, gsyn_tbl, db_fn):
    target_db = gffutils.FeatureDB(db_fn)
    for oid in orf_tbl:
        ctg, strand, gname, btype = orf_tbl[oid]
        if gname in gsyn_tbl:
            target_gname = gsyn_tbl[gname]
        else:
            target_gname = gname
        for gid in gene_tbl:
            if target_gname in gene_tbl[gid]:
                gene = target_db[gid]
                if gene.strand != strand:
                    print(oid)
                    print("gene strand: " + gene.strand)
                    print("orf strand: " + strand)
                    print(gid)
                    print(gname)


def proc_tx_omic_matches(args, match_tbl, orf_tbl, gene_tbl, tx_tbl, cds_tbl, db_fn, gsyn_tbl, is_partial=False):
    if is_partial:
        fn = os.path.join(args.out_dir, "partial_transcriptomic_matches.gff")
    else:
        fn = os.path.join(args.out_dir, "full_transcriptomic_matches.gff")
    fh = open(fn, 'w')
    target_db = gffutils.FeatureDB(db_fn)
    categ_tbl = dict()
    for oid in match_tbl:
        if not is_partial:
            is_categ_1 = True
        else:
            is_categ_1 = False
        is_categ_2 = False
        is_categ_3 = False
        is_nc = False
        match_lst = match_tbl[oid]
        ctg, strand, gname, btype = orf_tbl[oid]

        proc_match_lst = list()
        for hit in match_lst:
            match_tid = hit.ctg
            gid, tx = tx_tbl[match_tid]
            exon_lst = list(target_db.children(match_tid, featuretype='exon', order_by='start'))
            if strand == '-':
                exon_lst.reverse()
            hit_shift = hit.r_st
            hit_span = hit.r_en - hit.r_st
            st_coord = -1
            en_coord = -1
            splice_sites = list()
            for exon in exon_lst:
                exon_len = exon.stop - exon.start + 1
                if st_coord == -1:
                    if exon_len < hit_shift:
                        hit_shift -= exon_len
                    else:
                        # 1-based, closed (i.e., inclusive) start coordinate
                        st_coord = exon.start + hit_shift
                        # includes st_coord base
                        valid_exon_len = exon_len - hit_shift
                        if hit_span <= valid_exon_len:
                            # 1-based, closed (i.e., exclusive) end coordinate
                            en_coord = st_coord + hit_span - 1
                            break
                        else:
                            hit_span -= valid_exon_len
                            if strand == '-':
                                sd = exon.start
                            else:
                                sd = exon.stop
                else:
                    if strand == '-':
                        sa = exon.stop
                    else:
                        sa = exon.start
                    if strand == '-':
                        splice_sites.append((sa, sd))
                    else:
                        splice_sites.append((sd, sa))
                    if exon_len < hit_span:
                        hit_span -= exon_len
                        sd = exon.stop
                    else:
                        # 1-based, closed (i.e., exclusive) end coordinate
                        en_coord = exon.start + hit_span - 1
                        break
            if en_coord == -1:
                print("alignment end coordinate cannot be -1 after exons iteration.")
                sys.exit()

            # check if the hit_st and hit_en are contained in the 5' UTR region for uORFs
            # CDS coordinates were converted to 0-based, half-open-half-closed interval coordinate system
            # from 1-based, closed interval coordinate system.
            if match_tid not in cds_tbl:
                # if CDS was not recorded for this tx, then it's considered non-coding
                is_nc = True
            else:
                cds_st, cds_en = cds_tbl[match_tid]
                if btype == "uoORF":
                    # hit.r_st must be < cds_st
                    if hit.r_st < cds_st:
                        is_categ_2 = True
                else:
                    if hit.r_en <= cds_st:
                        is_categ_2 = True

            # check if the hit is contained in the gene locus annotated as associated with this orf

            target_gname_lst = [gname]
            if gname in gsyn_tbl:
                tmp = gsyn_tbl[gname]
                target_gname_lst.append(tmp)

            match_gname_lst = list(gene_tbl[gid])
            for item in target_gname_lst:
                if item in match_gname_lst:
                    is_categ_3 = True
                    break

            if strand == '-':
                proc_match_obj = (tx.chrom, tx.strand, gid, match_tid, en_coord, st_coord, splice_sites,
                                  is_categ_1, is_categ_2, is_categ_3, is_nc)
            else:
                proc_match_obj = (tx.chrom, tx.strand, gid, match_tid, st_coord, en_coord, splice_sites,
                                  is_categ_1, is_categ_2, is_categ_3, is_nc)
            if proc_match_obj in proc_match_lst:
                print("FATAL ERROR: duplicate entries in transcriptomic match table")
                sys.exit()

            print_tx_omic_match(fh, oid, proc_match_obj, target_gname_lst, match_gname_lst, hit.cigar_str)
            proc_match_lst.append(proc_match_obj)
        categ_tbl[oid] = (gname, target_gname_lst, proc_match_lst)
    fh.close()
    return categ_tbl


def annotate_gene_omic_matches(match_tbl, orf_tbl, db_fn):
    target_db = gffutils.FeatureDB(db_fn)
    annot_match_tbl = dict()
    for oid in match_tbl:
        ctg, strand, gname, btype = orf_tbl[oid]
        match_lst = match_tbl[oid]
        aug_match_lst = list()
        for hit in match_lst:
            use_tx = False
            skip = False
            feature_lst = list(target_db.region(seqid=hit.ctg, start=hit.r_st, end=hit.r_en,
                                                featuretype="gene"))
            if len(feature_lst) == 0:
                use_tx = True
                feature_lst = list(target_db.region(seqid=hit.ctg, start=hit.r_st, end=hit.r_en,
                                                    featuretype="transcript"))
            if len(feature_lst) == 0:
                skip = True
            gid_lst = list()
            if not skip:
                if len(feature_lst) == 1:
                    feature = feature_lst[0]
                    if feature.strand != strand:
                        skip = True
                        # continue
                    if not skip:
                        if not use_tx:
                            gid = feature.id
                        else:
                            gid = feature.attributes["Parent"][0]
                        gid_lst.append(gid)
                else:
                    if not use_tx:
                        for feature in feature_lst:
                            skip = False
                            if feature.strand != strand:
                                skip = True
                            if not skip:
                                if feature.start <= hit.r_st <= feature.stop and feature.start <= hit.r_en <= feature.stop:
                                    gid = feature.id
                                    gid_lst.append(gid)
                    else:
                        gid = None
                        for feature in feature_lst:
                            if gid is None:
                                gid = feature.attributes["Parent"][0]
                                gid_lst.append(gid)
                            else:
                                tmp = feature.attributes["Parent"][0]
                                # TODO: is there an edge case here?
                                if gid != tmp:
                                    print("WARNING: potential overlapping gene locus detected.")
                                    sys.exit()
            # only alignment to protein coding gene locus, on the same strand as uORF added
            if len(gid_lst) > 0:
                aug_match_lst.append((hit, gid_lst))
            else:
                aug_match_lst.append((hit, None))
        annot_match_tbl[oid] = aug_match_lst
    return annot_match_tbl


def proc_gene_omic_matches(args, match_tbl, orf_tbl, gene_tbl, gsyn_tbl, db_fn, is_partial=False):
    if is_partial:
        fn = os.path.join(args.out_dir, "partial_genomic_matches.gff")
    else:
        fn = os.path.join(args.out_dir, "full_genomic_matches.gff")
    fh = open(fn, 'w')
    target_db = gffutils.FeatureDB(db_fn)
    categ_tbl = dict()
    for oid in match_tbl:
        if not is_partial:
            is_categ_1 = True
        else:
            is_categ_1 = False
        match_lst = match_tbl[oid]
        ctg, strand, gname, btype = orf_tbl[oid]
        hit_tbl = dict()

        target_gname_lst = [gname]
        if gname in gsyn_tbl:
            tmp = gsyn_tbl[gname]
            target_gname_lst.append(tmp)
        # gid_lst needed because I take all features overlapping the aligned coordinates
        for hit, gid_lst in match_lst:
            categ_2_lst = list()
            categ_3_lst = list()
            if gid_lst is None:
                hit_tbl[hit] = (None, None)
                continue
            for gid in gid_lst:
                is_categ_2 = False
                is_nc = False
                gene = target_db[gid]
                try:
                    if gene.attributes['gene_biotype'][0] != "protein_coding":
                        is_nc = True
                except:
                    if gene.attributes['gene_type'][0] != "protein_coding":
                        is_nc = True
                source_tx = None
                if not is_nc:
                    tx_lst = list(target_db.children(gid, featuretype='transcript'))
                    if len(tx_lst) == 0:
                        print("something might have gone wrong...")
                    for tx in tx_lst:
                        if strand != tx.strand:
                            print("something might have gone wrong with: " + oid)
                        cds_lst = list(target_db.children(tx.id, featuretype='CDS', order_by='start'))
                        exon_lst = list(target_db.children(tx.id, featuretype='exon', order_by='start'))
                        fpu = None
                        if len(cds_lst) != 0:
                            for cds in cds_lst:
                                # 0-based, half-open, half-closed interval coordinate system
                                if strand == '+':
                                    cand_1 = (tx.start - 1, cds.start)
                                    break
                                else:
                                    cand_1 = (cds.end, tx.end)
                            for exon in exon_lst:
                                # 0-based, half-open, half-closed interval coordinate system
                                if strand == '+':
                                    cand_2 = (tx.start - 1, exon.end)
                                    break
                                else:
                                    cand_2 = (exon.start, tx.end)
                            if strand == '+':
                                fpu = (tx.start - 1, min(cand_1[1], cand_2[1]))
                            else:
                                fpu = (max(cand_1[0], cand_2[0]), tx.end)

                        # OLD STUB
                        # if len(cds_lst) == 0:
                        #     fpu = None
                        # for cds in cds_lst:
                        #     # 0-based, half-open, half-closed interval coordinate system
                        #     if strand == '+':
                        #         fpu = (tx.start - 1, cds.start)
                        #         break
                        #     else:
                        #         fpu = (cds.end, tx.end)

                        if fpu is not None:
                            fpu_st, fpu_en = fpu
                            if btype == 'uORF':
                                if fpu_st <= hit.r_st < hit.r_en <= fpu_en:
                                    is_categ_2 = True
                                    # this is different from how tx_omic matches are processed
                                    # for tx_omic matches, even if the match isn't contained within
                                    # the reference tx, its id will be outputted bc alignment is
                                    # specifically to that tx; genomic alignment is to a genomic region,
                                    # on the contrary.
                                    source_tx = tx
                                    break
                            else:
                                if fpu_st <= hit.r_st:
                                    is_categ_2 = True
                                    source_tx = tx
                                    break
                # categ_2 = is the match properly positioned w.r.t a 5' UTR?
                categ_2_lst.append((gid, is_nc, is_categ_2, source_tx))
                is_categ_3 = False

                match_gname_lst = list(gene_tbl[gid])

                for item in target_gname_lst:
                    if item in match_gname_lst:
                        is_categ_3 = True
                        break
                # categ_3 = is the downstream gene same as the annotated associate gene?
                categ_3_lst.append((gid, is_categ_3, match_gname_lst))
            hit_tbl[hit] = (categ_2_lst, categ_3_lst)
        print_gene_omic_match(fh, oid, hit_tbl, target_gname_lst, is_categ_1)
        categ_tbl[oid] = (is_categ_1, target_gname_lst, hit_tbl)
    return categ_tbl


def print_gene_omic_match(fh, oid, hit_tbl, src_gname_lst, is_categ_1):
    for hit in hit_tbl:
        if hit.strand == -1:
            strand = '-'
        else:
            strand = '+'
        categ_2_lst, categ_3_lst = hit_tbl[hit]
        if categ_2_lst is None:
            if categ_3_lst is not None:
                print("something is weird")
                sys.exit()
            fh.write(hit.ctg + "\tuORFinder\taligned_segment\t" + str(hit.r_st + 1) + "\t" + str(hit.r_en) +
                     "\t.\t" + strand + "\t.\t" + "ID=" + oid + ";source_gene_name=")
            for i in range(len(src_gname_lst)):
                gname = src_gname_lst[i]
                if i == len(src_gname_lst) - 1:
                    fh.write(gname)
                else:
                    fh.write(gname + ",")
            fh.write(";category_1=" + str(is_categ_1) + ";category_2=False;category_3=False\n")
            if not is_categ_1:
                fh.write(";partial_cigar_string=" + hit.cigar_str)
            continue
        for i in range(len(categ_2_lst)):
            gid, is_nc, is_categ_2, src_tx = categ_2_lst[i]
            if is_categ_2 and src_tx is None:
                print("this cannot happen")
            gid_alt, is_categ_3, tgt_gname_lst = categ_3_lst[i]
            fh.write(hit.ctg + "\tuORFinder\taligned_segment\t" + str(hit.r_st + 1) + "\t" + str(hit.r_en) +
                     "\t.\t" + strand + "\t.\t" + "ID=" + oid + ";source_gene_name=")
            for j in range(len(src_gname_lst)):
                gname = src_gname_lst[j]
                if j == len(src_gname_lst) - 1:
                    fh.write(gname)
                else:
                    fh.write(gname + ",")
            fh.write(";category_1=" + str(is_categ_1))
            fh.write(";target_gene_id=" + gid + ";non_coding=" + str(is_nc) + ";category_2=" + str(is_categ_2))
            if src_tx is not None:
                fh.write(";target_transcript_id=" + src_tx.id)
            if gid != gid_alt:
                print("this also cannot happen")
            fh.write(";category_3=" + str(is_categ_3) + ";target_gene_name=")

            for j in range(len(tgt_gname_lst)):
                gname = tgt_gname_lst[j]
                if j == len(tgt_gname_lst) - 1:
                    fh.write(gname)
                else:
                    fh.write(gname + ",")

            if not is_categ_1:
                fh.write(";partial_cigar_string=" + hit.cigar_str)
            fh.write("\n")


def print_tx_omic_match(fh, oid, proc_match_obj, src_gname_lst, match_gname_lst, cigar):
    ctg, strand, match_gid, match_tid, match_st, match_en, splice_sites, \
        is_categ_1, is_categ_2, is_categ_3, is_nc = proc_match_obj
    fh.write(ctg + "\tuORFinder\taligned_segment\t" + str(match_st) + "\t" + str(match_en) + "\t.\t" + strand +
             "\t.\t" + "ID=" + oid + ";source_gene_name=")
    for i in range(len(src_gname_lst)):
        gname = src_gname_lst[i]
        if i == len(src_gname_lst) - 1:
            fh.write(gname)
        else:
            fh.write(gname + ",")
    fh.write(";category_1=" + str(is_categ_1))
    fh.write(";target_gene_id=" + match_gid + ";target_transcript_id=" + match_tid + ";non_coding=" + str(is_nc))
    fh.write(";category_2=" + str(is_categ_2) + ";category_3=" + str(is_categ_3) + ";target_gene_name=")
    for i in range(len(match_gname_lst)):
        gname = match_gname_lst[i]
        if i == len(match_gname_lst) - 1:
            fh.write(gname)
        else:
            fh.write(gname + ",")
    if not is_categ_1:
        fh.write(";partial_cigar_string=" + cigar)
    fh.write("\n")
    ctr = 0
    for st, en in splice_sites:
        fh.write(ctg + "\tuORFinder\tjunction\t" + str(st) + "\t" + str(en) + "\t.\t" + strand +
                 "\t.\t" + "ID=" + oid + "_exon_" + str(ctr) + ";Parent=" + match_tid + "\n")
        ctr += 1


def compile_lvl(categ_tbl, max_lvl_tbl, full_lvl_tbl, use_tx=False):
    for oid in categ_tbl:
        if use_tx:
            _, _, proc_match_lst = categ_tbl[oid]
            for match in proc_match_lst:
                _, _, _, tid, _, _, _, is_1, is_2, is_3, _ = match
                tmp = combine_categ_scores(is_1, is_2, is_3)
                if oid in max_lvl_tbl:
                    if oid not in full_lvl_tbl:
                        print("something could have gone wrong...")
                    curr_max, curr_ref, did_use_tx = max_lvl_tbl[oid]
                    max_lvl = max(curr_max, tmp)
                    if max_lvl != curr_max:
                        max_lvl_tbl[oid] = (tmp, tid, use_tx)
                    full_lvl_tbl[oid].append((tmp, tid, use_tx))
                else:
                    max_lvl_tbl[oid] = (tmp, tid, use_tx)
                    full_lvl_tbl[oid] = [(tmp, tid, use_tx)]
        else:
            is_1, src_gname_lst, hit_tbl = categ_tbl[oid]
            for hit in hit_tbl:
                categ_2_lst, categ_3_lst = hit_tbl[hit]
                if categ_2_lst is None:
                    # sanity check
                    if categ_3_lst is not None:
                        print("something is wrong")
                        sys.exit()
                    if is_1:
                        tmp = 4
                    else:
                        tmp = 0
                    if oid in max_lvl_tbl:
                        # sanity check
                        if oid not in full_lvl_tbl:
                            print("something could have gone wrong...")
                            sys.exit()
                        curr_max, curr_ref, did_use_tx = max_lvl_tbl[oid]
                        max_lvl = max(curr_max, tmp)
                        if max_lvl != curr_max:
                            max_lvl_tbl[oid] = (tmp, None, use_tx)
                        full_lvl_tbl[oid].append((tmp, None, use_tx))
                    else:
                        max_lvl_tbl[oid] = (tmp, None, use_tx)
                        full_lvl_tbl[oid] = [(tmp, None, use_tx)]
                    continue
                for i in range(len(categ_2_lst)):
                    gid, _, is_2, _ = categ_2_lst[i]
                    _, is_3, _ = categ_3_lst[i]
                    tmp = combine_categ_scores(is_1, is_2, is_3)
                    if oid in max_lvl_tbl:
                        if oid not in full_lvl_tbl:
                            print("something could have gone wrong...")
                        curr_max, curr_ref, did_use_tx = max_lvl_tbl[oid]
                        max_lvl = max(curr_max, tmp)
                        if max_lvl != curr_max:
                            max_lvl_tbl[oid] = (tmp, gid, use_tx)
                        full_lvl_tbl[oid].append((tmp, gid, use_tx))
                    else:
                        max_lvl_tbl[oid] = (tmp, gid, use_tx)
                        full_lvl_tbl[oid] = [(tmp, gid, use_tx)]


def combine_categ_scores(is_1, is_2, is_3):
    if is_1 and is_2 and is_3:
        score = 4 + 2 + 1
    elif is_1 and is_2:
        score = 4 + 2
    elif is_1 and is_3:
        score = 4 + 1
    elif is_2 and is_3:
        score = 2 + 1
    elif is_1:
        score = 4
    elif is_2:
        score = 2
    elif is_3:
        score = 1
    else:
        score = 0
    return score


def print_lvl(fh, oid, lvl, g_or_tid, use_tx):
    if g_or_tid is None and use_tx is None:
        fh.write(oid + "\t" + str(lvl) + "\t.\t.\n")
    elif g_or_tid is None:
        fh.write(oid + "\t" + str(lvl) + "\t.\t" + str(use_tx) + "\n")
    else:
        fh.write(oid + "\t" + str(lvl) + "\t" + g_or_tid + "\t" + str(use_tx) + "\n")


def main():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-ogtf", '--orf-gtf', type=str, help="", required=True)
    parser.add_argument("-ofa", '--orf-fasta', type=str, help="", required=True)
    parser.add_argument("-g", '--genome', type=str, help="", required=True)
    parser.add_argument("-tgtf", '--target-gtf', type=str, help="", required=True)
    parser.add_argument("-tfa", '--target-fasta', type=str, help="", required=True)
    parser.add_argument("-tdb", '--target-db', type=str, help="", required=False, default=None)
    parser.add_argument('-tmp', '--tmp-dir', type=str, help="", required=False, default="./tmp")
    parser.add_argument('-ti', '--transcript_index', type=str, help="", required=False, default=None)
    parser.add_argument('-gi', '--genome_index', type=str, help="", required=False, default=None)
    parser.add_argument('-syn', type=str, help="", required=False, default=None)
    parser.add_argument('--save', default=False, help="", required=False, action='store_true')
    parser.add_argument('--resume', default=False, help="", required=False, action='store_true')
    parser.add_argument('-t', '--threads', type=int, required=False, default=1)
    parser.add_argument('-o', '--out-dir', type=str, required=True)
    parser.add_argument('-k', '--kmer-size', type=int, required=False, default=15)
    parser.add_argument('--rematch', default=False, help="", required=False, action='store_true')
    args = parser.parse_args()

    if args.target_db is None:
        db_fn = os.path.join(args.tmp_dir, "target_db")
        print("target db argument not provided. will be looking for it in tmp location: %s" % db_fn)
        if not os.path.exists(args.tmp_dir):
            print("tmp directory wasn't detected. creating it...")
            os.makedirs(args.tmp_dir)
    else:
        db_fn = args.target_db

    if not os.path.isfile(db_fn):
        print("target db not found. creating it...")
        success = False
        try:
            gffutils.create_db(args.target_gtf, db_fn, id_spec="ID", merge_strategy="create_unique", keep_order=True)
            success = True
        except:
            print("initial gffutils strategy failed... trying an alternative")
        if not success:
            rm_cmd = "rm -rf " + db_fn
            call(rm_cmd, shell=True)
            print("removed previously constructed target db")
            try:
                gffutils.create_db(args.target_gtf, db_fn, id_spec={"gene": "ID", "transcript": "ID"},
                                   merge_strategy="create_unique", keep_order=True)
            except:
                print("exhausted all gffutils strategies... contact developers for a solution.")
        print("successfully created target db.")
    else:
        print("target db found.")

    if args.syn is not None:
        print("loading gene synonyms...")
        gsyn_tbl = load_gene_syn(args)

    if args.resume:
        print("resuming previous run by loading target tables from: %s" % args.tmp_dir)
        if not os.path.exists(args.tmp_dir):
            print("target index must be present tmp dir; tmp dir was not detected")
            sys.exit()
        save_fn = os.path.join(args.tmp_dir, "target.index")
        if not os.path.isfile(save_fn):
            print("target index was not detected. please check that target.index is present in your tmp dir")
            sys.exit()
        with open(save_fn, "rb") as f:
            gene_tbl, tx_tbl, gene_set, cds_tbl = pickle.load(f)
    else:
        print("building target transcriptome tables")
        gene_tbl, tx_tbl, gene_set = build_target_tbls(db_fn)

        fn = os.path.join(args.out_dir, "no_gene_names.txt")
        fh = open(fn, 'w')
        for gene in no_gname:
            fh.write(gene + "\n")
        fh.close()

        print("obtaining cds information from target transcriptome")
        cds_tbl = get_cds(args)

        if args.save:
            print("saving target transcriptome tables to tmp dir: %s" % args.tmp_dir)
            if not os.path.exists(args.tmp_dir):
                print("tmp directory wasn't detected. creating it...")
                os.makedirs(args.tmp_dir)
            save_fn = os.path.join(args.tmp_dir, "target.index")
            with open(save_fn, "wb") as f:
                pickle.dump([gene_tbl, tx_tbl, gene_set, cds_tbl], f)

    print("loading query orfs")
    orf_tbl = load_orfs(args)

    print("finding gene names not matched to a gene id in target db")
    # you can still run this stub when --resume is set
    if args.rematch:
        print("gene set size BEFORE re-matching: " + str(len(gene_set)))
        missing_gene_set = check_missing_gene(orf_tbl, gene_set, gsyn_tbl, gene_tbl, db_fn)
        print("gene set size AFTER re-matching: " + str(len(gene_set)))
    else:
        missing_gene_set = check_missing_gene(orf_tbl, gene_set, gsyn_tbl)

    fn = os.path.join(args.out_dir, "missing_genes.txt")
    fh = open(fn, 'w')
    for gene in missing_gene_set:
        fh.write(gene + "\n")
    fh.close()

    # sanity check: if ORFs are located on the same strand as the associated genes
    # compare_orf_to_gene_strand(orf_tbl, gene_tbl, gsyn_tbl, db_fn)
    # sys.exit()

    print("1st pass - mapping query orf sequences to target genome sequences")
    full_match_tbl_1, partial_match_tbl_1 = map_to_gene_ome(args, orf_tbl)

    print("annotating genomic alignments with gene locus information")
    annot_full_match_tbl_1 = annotate_gene_omic_matches(full_match_tbl_1, orf_tbl, db_fn)
    annot_partial_match_tbl_1 = annotate_gene_omic_matches(partial_match_tbl_1, orf_tbl, db_fn)

    print("processing genomic alignments")
    full_categ_tbl_1 = proc_gene_omic_matches(args, annot_full_match_tbl_1, orf_tbl, gene_tbl, gsyn_tbl, db_fn)
    partial_categ_tbl_1 = proc_gene_omic_matches(args, annot_partial_match_tbl_1, orf_tbl, gene_tbl, gsyn_tbl,
                                                 db_fn, is_partial=True)

    print("2nd pass - mapping query orf sequences to target transcriptome sequences")
    full_match_tbl_2, partial_match_tbl_2 = map_to_tx_ome(args)
    # print(len(full_match_tbl_1) + len(full_match_tbl_2))

    print("processing transcriptomic alignments")
    full_categ_tbl_2 = proc_tx_omic_matches(args, full_match_tbl_2, orf_tbl, gene_tbl, tx_tbl, cds_tbl, db_fn, gsyn_tbl)
    partial_categ_tbl_2 = proc_tx_omic_matches(args, partial_match_tbl_2, orf_tbl, gene_tbl, tx_tbl, cds_tbl, db_fn,
                                               gsyn_tbl, is_partial=True)

    print("combining genomic + transcriptomic alignment results to output ")

    max_lvl_tbl = dict()
    full_lvl_tbl = dict()
    compile_lvl(full_categ_tbl_1, max_lvl_tbl, full_lvl_tbl)
    compile_lvl(full_categ_tbl_2, max_lvl_tbl, full_lvl_tbl, use_tx=True)
    compile_lvl(partial_categ_tbl_1, max_lvl_tbl, full_lvl_tbl)
    compile_lvl(partial_categ_tbl_2, max_lvl_tbl, full_lvl_tbl, use_tx=True)

    max_lvl_fn = os.path.join(args.out_dir, "max_lvls.tsv")
    max_lvl_fh = open(max_lvl_fn, 'w')
    full_lvl_fn = os.path.join(args.out_dir, "full_lvls.tsv")
    full_lvl_fh = open(full_lvl_fn, 'w')
    synonymous_mut_fn = os.path.join(args.out_dir, "syn_mut_info.txt")
    synonymous_mut_fh = open(synonymous_mut_fn, 'w')

    seven = 0
    six = 0
    five = 0
    four = 0
    three = 0
    two = 0
    one = 0
    zero = 0

    for oid in max_lvl_tbl:
        if oid not in full_lvl_tbl:
            print("something might have gone wrong somewhere...")
        lvl, g_or_tid, use_tx = max_lvl_tbl[oid]

        if lvl == 7:
            seven += 1
        elif lvl == 6:
            six += 1
        elif lvl == 5:
            five += 1
        elif lvl == 4:
            four += 1
        elif lvl == 3:
            three += 1
        elif lvl == 2:
            two += 1
        elif lvl == 1:
            one += 1
        elif lvl == 0:
            zero += 1
        else:
            print("FATAL ERROR: unrecognized level detected for: " + oid)
            sys.exit()
        print_lvl(max_lvl_fh, oid, lvl, g_or_tid, use_tx)
        full_lvl_lst = full_lvl_tbl[oid]

        # sanity check
        if oid in synonymous_mut:
            if lvl < 4:
                print("can't be true")
                sys.exit()
            k = 0
            for lvl, _, _ in full_lvl_lst:
                if lvl >= 4:
                    k += 1

            if len(synonymous_mut[oid]) != k:
                for rec in synonymous_mut[oid]:
                    synonymous_mut_fh.write(oid + "\t" + rec.ctg + "\t" + str(rec.r_st + 1) + "\t" + str(rec.r_en) +
                                            "\t" + str(rec.strand) + "\n")
                # print(oid + "\t" + str(len(synonymous_mut[oid])) + "\t" + str(k))

        for lvl, g_or_tid, use_tx in full_lvl_lst:
            print_lvl(full_lvl_fh, oid, lvl, g_or_tid, use_tx)

    print("# of orfs with protein-level sequence match: " + str(len(synonymous_mut)))

    found_oids = set(max_lvl_tbl.keys())
    tot_oids = set(orf_tbl.keys())
    not_found_oids = tot_oids.difference(found_oids)

    for oid in not_found_oids:
        zero += 1
        print_lvl(max_lvl_fh, oid, 0, None, None)
        print_lvl(full_lvl_fh, oid, 0, None, None)

    max_lvl_fh.close()
    full_lvl_fh.close()
    synonymous_mut_fh.close()

    print("# of level 7: " + str(seven))
    print("# of level 6: " + str(six))
    print("# of level 5: " + str(five))
    print("# of level 4: " + str(four))
    print("# of level 3: " + str(three))
    print("# of level 2: " + str(two))
    print("# of level 1: " + str(one))
    print("# of level 0: " + str(zero))


if __name__ == "__main__":
    main()

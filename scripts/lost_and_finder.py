#!/usr/bin/env python

import argparse
import mappy as mp
import pyfastx
import os
import sys

USE_EQX = 0x4000000


def get_args():
    parser = argparse.ArgumentParser(description="")
    parser.add_argument("-g", "--genome", type=str, help="", required=True)
    parser.add_argument("-gi", "--genome-index", type=str, help="", required=False, default=None)
    parser.add_argument("-tx", "--transcripts", type=str, help="", required=True)
    parser.add_argument("-gname", "--gene-name", type=str, help="", required=True)
    parser.add_argument("-o", "--out-dir", type=str, help="", required=False, default="./")
    # parser.add_argument("-k", "--kmer-size", type=int, help="", required=False)
    parser.add_argument("-t", "--threads", type=int, help="", required=False, default=1)
    parser.add_argument("--strand", type=int, help="", required=False, default=None)
    args = parser.parse_args()
    return args


def map_to_genome(args):
    if args.genome_index is None:
        index_fn = args.genome + ".mmi"
        print("mm2 genome index argument not provided")
    else:
        index_fn = args.genome_index
    print("minimap index: " + index_fn)

    if not os.path.isfile(index_fn):
        print("mm2 genome index not found. creating it...")
        a = mp.Aligner(fn_idx_in=args.genome, fn_idx_out=index_fn, preset="splice", n_threads=args.threads)
    else:
        a = mp.Aligner(fn_idx_in=index_fn, preset="splice", n_threads=args.threads)

    aln_fn = args.out_dir + args.gene_name + "_aln.tsv"
    aln_fh = open(aln_fn, 'w')
    aln_fh.write("transcript_id\tchr\tstrand\tref_start\tref_end\tCIGAR\n")

    aln_tbl = dict()
    print("aligning transcript sequences to the target genome...")

    for name, seq, qual in mp.fastx_read(args.transcripts):
        mapq = None
        out_ln = None
        keep_looking = False
        for hit in a.map(seq, cs=True, extra_flags=USE_EQX):
            # just using the primary alignments (might not be correct / accurate)
            if hit.is_primary:
                if args.strand is not None:
                    if hit.strand != args.strand:
                        keep_looking = True
                        continue
                aln_fh.write(name + "\t" + hit.ctg + "\t" + str(hit.strand) + "\t" + str(hit.r_st) + "\t"
                             + str(hit.r_en) + "\t" + hit.cigar_str + "\t")
                aln_tbl[name] = (hit.ctg, hit.strand, hit.r_st, hit.r_en)
            elif keep_looking:
                if mapq is None:
                    out_ln = name + "\t" + hit.ctg + "\t" + str(hit.strand) + "\t" + str(hit.r_st) + "\t" \
                             + str(hit.r_en) + "\t" + hit.cigar_str + "\t"
                    aln_tbl[name] = (hit.ctg, hit.strand, hit.r_st, hit.r_en)
                else:
                    if hit.mapq > mapq:
                        out_ln = name + "\t" + hit.ctg + "\t" + str(hit.strand) + "\t" + str(hit.r_st) + "\t" \
                                 + str(hit.r_en) + "\t" + hit.cigar_str + "\t"
                        aln_tbl[name] = (hit.ctg, hit.strand, hit.r_st, hit.r_en)
        if keep_looking:
            aln_fh.write(out_ln)

    aln_fh.close()
    return aln_tbl


def define_gene_coord(aln_tbl, args):
    g_st = None
    g_en = None
    g_ctg = None
    g_strand = None

    print("iterating through the alignment records to finalized gene coordinates...")
    for tid in aln_tbl:
        ctg, strand, r_st, r_en = aln_tbl[tid]
        if g_st is None:
            # sanity check
            if g_en is not None or g_ctg is not None or g_strand is not None:
                print("something went wrong1")
                sys.exit()
            g_ctg = ctg
            g_strand = strand
            g_st = r_st
            g_en = r_en
        else:
            # sanity check
            if ctg != g_ctg or strand != g_strand:
                print("something went wrong2")
                print(g_ctg)
                print(g_strand)
                print(ctg)
                print(strand)
                print(r_st)
                print(r_en)
                print(tid)
                sys.exit()
            g_st = min(g_st, r_st)
            g_en = max(g_en, r_en)

    gene_fn = args.out_dir + args.gene_name + "_coords.txt"
    gene_fh = open(gene_fn, 'w')
    gene_fh.write("chr\tstrand\tstart\tend\n")
    gene_fh.write(g_ctg + "\t" + str(g_strand) + "\t" + str(g_st) + "\t" + str(g_en) + "\n")
    gene_fh.close()


def main():
    args = get_args()
    if not os.path.exists(args.out_dir):
        print("output directory not found. making it...")
        os.makedirs(args.out_dir)

    aln_tbl = map_to_genome(args)
    define_gene_coord(aln_tbl, args)


if __name__ == "__main__":
    main()





#!/usr/bin/env python

import sys


def write_gtf_line(gtf_fh, ctg, feature, start, end, strand, gid, oid, tid, ref_id, splice=False, sites=None):
    if splice and sites is None:
        print("Provide splice sites when writing out a gtf line with a candidate splice site.")
        sys.exit()
    gtf_fh.write(str(ctg) + "\tuORFexplorer\t" + feature + "\t" + str(start) + "\t" + str(end) + "\t.\t" + strand +
                 "\t.\t" + 'gene_id "' + gid + '"; ' + 'orf_id "' + oid + '"; ' + 'transcript_id "' + tid + '"; ' +
                 'reference_id "' + ref_id + '"; ')
    if splice:
        gtf_fh.write("splice sites (" + str(sites[0]) + ", " + str(sites[1]) + ");\n")
    else:
        gtf_fh.write("\n")


def load_gene_syn(fn):
    fh = open(fn, 'r')
    gsyn_tbl = dict()
    for line in fh:
        k = line.split("\t")[0]
        v = line.split("\t")[1].strip()
        gsyn_tbl[k] = v
    fh.close()
    return gsyn_tbl


def translate(seq):
    aa_seq = ""
    for i in range(0, len(seq), 3):
        triplet = seq[i:i+3]
        amino = codons[triplet][0]
        aa_seq += amino
    return aa_seq


# truncated = 1, elongated = 2, none = 0
def compare(query, target):
    end = min(len(query), len(target))
    match = 0
    res = None
    for i in range(0, end):
        match += 1
        if query[i] != target[i]:
            res = 0
            break
    if res != 0:
        if len(query) < len(target):
            res = 1
        else:
            res = 2
    return match, res


# start = 1; stop = 2; neither = 0
codons = {
    "TTT": ("F", 0),
    "TTC": ("F", 0),
    "TTA": ("L", 0),
    "TTG": ("L", 0),
    "TCT": ("S", 0),
    "TCC": ("S", 0),
    "TCA": ("S", 0),
    "TCG": ("S", 0),
    "TAT": ("Y", 0),
    "TAC": ("Y", 0),
    "TAA": ("*", 2),
    "TAG": ("*", 2),
    "TGT": ("C", 0),
    "TGC": ("C", 0),
    "TGA": ("*", 2),
    "TGG": ("W", 0),
    "CTT": ("L", 0),
    "CTC": ("L", 0),
    "CTA": ("L", 0),
    "CTG": ("L", 0),
    "CCT": ("P", 0),
    "CCC": ("P", 0),
    "CCA": ("P", 0),
    "CCG": ("P", 0),
    "CAT": ("H", 0),
    "CAC": ("H", 0),
    "CAA": ("Q", 0),
    "CAG": ("Q", 0),
    "CGT": ("R", 0),
    "CGC": ("R", 0),
    "CGA": ("R", 0),
    "CGG": ("R", 0),
    "ATT": ("I", 0),
    "ATC": ("I", 0),
    "ATA": ("I", 0),
    "ATG": ("M", 1),
    "ACT": ("T", 0),
    "ACC": ("T", 0),
    "ACA": ("T", 0),
    "ACG": ("T", 0),
    "AAT": ("N", 0),
    "AAC": ("N", 0),
    "AAA": ("K", 0),
    "AAG": ("K", 0),
    "AGT": ("S", 0),
    "AGC": ("S", 0),
    "AGA": ("R", 0),
    "AGG": ("R", 0),
    "GTT": ("V", 0),
    "GTC": ("V", 0),
    "GTA": ("V", 0),
    "GTG": ("V", 0),
    "GCT": ("A", 0),
    "GCC": ("A", 0),
    "GCA": ("A", 0),
    "GCG": ("A", 0),
    "GAT": ("D", 0),
    "GAC": ("D", 0),
    "GAA": ("E", 0),
    "GAG": ("E", 0),
    "GGT": ("G", 0),
    "GGC": ("G", 0),
    "GGA": ("G", 0),
    "GGG": ("G", 0)
}
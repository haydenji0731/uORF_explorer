#!/usr/bin/env python

import sys
from tqdm import tqdm


def get_pairs(ln):
    kv_pairs = ln.split(";")
    out = list()
    for p in kv_pairs:
        tmp = p.strip().split(" ")
        if len(tmp) < 2:
            continue
        k = tmp[0].strip()
        v = tmp[1].strip().replace('"', '')
        out.append((k, v))
    return out


def load_utr(fn):
    fh = open(fn, 'r')
    cds_tbl = dict()
    exon_tbl = dict()
    for ln in tqdm(fh):
        if ln[0] == '#':
            continue
        clean_ln = ln.strip()
        fields = clean_ln.split("\t")
        st_pos = int(fields[3])
        en_pos = int(fields[4])
        strand = fields[6]
        if fields[2] == "transcript":
            ps = get_pairs(fields[8])
            for k, v in ps:
                if k == "transcript_id":
                    tid = v
                    cds_tbl[tid] = list()
                    exon_tbl[tid] = list()
        elif fields[2] == "exon":
            exon_tbl[tid].append((st_pos, en_pos, strand))
        elif fields[2] == "CDS":
            cds_tbl[tid].append((st_pos, en_pos, strand))
    fh.close()

    utr_tbl = dict()
    # no cds, not protein coding --> ignore
    for tid in cds_tbl:
        if len(cds_tbl[tid]) == 0:
            continue
        exon_lst = sorted(exon_tbl[tid], key=lambda x: x[0])
        cds_lst = sorted(cds_tbl[tid], key=lambda x: x[0])
        fc_st, fc_en, strand = cds_lst[0]
        lc_st, lc_en, _ = cds_lst[len(cds_lst) - 1]
        fc_idx = None
        lc_idx = None
        lc_spliced = False
        fc_spliced = False
        for i in range(len(exon_lst)):
            e_st, e_en, _ = exon_lst[i]
            if e_st <= fc_st and fc_en <= e_en:
                if e_st == fc_st and fc_en == e_en:
                    fc_spliced = True
                fc_idx = i
                break
        exon_lst_rev = list(reversed(exon_lst))
        for j in range(len(exon_lst_rev)):
            e_st, e_en, _ = exon_lst_rev[j]
            if e_st <= lc_st and lc_en <= e_en:
                if e_st == lc_st and lc_en == e_en:
                    lc_spliced = True
                lc_idx = j
                break

        if fc_idx is None:
            print("shouldn't happen")
            sys.exit()
        f_utr_lst = exon_lst[:fc_idx]

        if not fc_spliced:
            f_utr_lst.append((exon_lst[fc_idx][0], fc_st - 1, strand))

        if lc_idx is None:
            print("shouldn't happen")
            sys.exit()
        l_utr_lst = exon_lst_rev[:lc_idx]

        if not lc_spliced:
            l_utr_lst.append((lc_en + 1, exon_lst_rev[lc_idx][1], strand))
        if strand == '+':
            l_utr_lst.reverse()
            # five prime, three prime
            # increasing coords
            utr_tbl[tid] = (f_utr_lst, l_utr_lst)
        else:
            f_utr_lst.reverse()
            # five prime, three prime
            # decreasing coords
            utr_tbl[tid] = (l_utr_lst, f_utr_lst)
        # if tid == "ENST00000226319.11":
        #     print(lc_idx)
        #     print(utr_tbl[tid])
    return utr_tbl


def load_qry(fn):
    fh = open(fn, 'r')
    q2r_tbl = dict()
    exon_tbl = dict()
    for ln in fh:
        clean_ln = ln.strip()
        fields = clean_ln.split("\t")
        if fields[2] == "transcript":
            ps = get_pairs(fields[8])
            strand = fields[6]
            for k, v in ps:
                if k == "transcript_id":
                    tid = v
                elif k == "reference_id":
                    rid = v
                    q2r_tbl[tid] = rid, strand
                    exon_tbl[tid] = list()
                    break
        elif fields[2] == "exon":
            e_st = int(fields[3])
            e_en = int(fields[4])
            strand = fields[6]
            exon_tbl[tid].append((e_st, e_en, strand))
    fh.close()
    return q2r_tbl, exon_tbl


def add_utr(utr_tbl, q2r_tbl, exon_tbl):
    add_utr_tbl = dict()
    for tid in exon_tbl:
        rid, strand = q2r_tbl[tid]
        fp_utr_lst, tp_utr_lst = utr_tbl[rid]
        if strand == '-':
            exon_lst = sorted(exon_tbl[tid], key=lambda x: x[0], reverse=True)
        else:
            exon_lst = sorted(exon_tbl[tid], key=lambda x: x[0])
        fe_st, fe_en, _ = exon_lst[0]
        te_st, te_en, _ = exon_lst[len(exon_lst) - 1]

        fp_idx = None
        tp_idx = None
        fp_spliced = False

        for i in range(len(fp_utr_lst)):
            u_st, u_en, _ = fp_utr_lst[i]
            if u_st <= fe_st and fe_en <= u_en:
                # print((tid, rid, strand, u_st, fe_st, fe_en, u_en))
                if strand == '+':
                    if u_st == fe_st:
                        fp_spliced = True
                else:
                    if fe_en == u_en:
                        fp_spliced = True
                fp_idx = i
                break

        for j in range(len(tp_utr_lst)):
            u_st, u_en, _ = tp_utr_lst[j]
            if strand == '+':
                if u_st == te_en - 2:
                    tp_idx = j
                    break
            else:
                if u_en == te_st + 2:
                    tp_idx = j
                    break

        if fp_idx is None:
            add_fp_utr = []
        else:
            add_fp_utr = fp_utr_lst[:fp_idx]
            if not fp_spliced:
                if strand == '+':
                    add_fp_utr.append((fp_utr_lst[fp_idx][0], fe_st - 1, strand))
                else:
                    add_fp_utr.append((fe_en + 1, fp_utr_lst[fp_idx][1], strand))
        # sanity check
        # print((tid, rid, fp_utr_lst, fp_idx, add_fp_utr, fp_spliced))

        if tp_idx is None:
            # possible that the reference tx doesn't have 3' UTR annotated
            add_tp_utr = []
        else:
            add_tp_utr = tp_utr_lst[tp_idx:]
        # sanity check
        # print((tid, rid, tp_utr_lst, tp_idx, add_tp_utr))

        add_utr_tbl[tid] = (add_fp_utr, add_tp_utr)
    return add_utr_tbl


def add_features(strand, exon_lst, add_fp_utr, add_tp_utr):
    # sanity check; exons sorted correctly?
    # if strand == '+':
    #     tmp = sorted(exon_lst, key=lambda x: x[0])
    #     assert tmp == exon_lst
    # else:
    #     tmp = sorted(exon_lst, key=lambda x: x[0], reverse=True)
    #     assert tmp == exon_lst
    feature_lst = list()
    fp_utr = [(tpl + (0,)) for tpl in add_fp_utr]
    tp_utr = [(tpl + (0,)) for tpl in add_tp_utr]
    exons = [(tpl + (1,)) for tpl in exon_lst]
    if len(fp_utr) > 0:
        # print(add_fp_utr)
        feature_lst += fp_utr[:-1]
        if strand == '+':
            feature_lst.append((fp_utr[-1][0], exons[0][1], strand, 0))
        else:
            feature_lst.append((exons[0][0], fp_utr[-1][1], strand, 0))
        feature_lst.append((exons[0][0], exons[0][1], strand, 1))
        feature_lst += exons[1:-1]
    else:
        feature_lst = exons[:-1]

    if len(tp_utr) > 0:
        # print(add_tp_utr)
        if strand == '+':
            feature_lst.append((exons[-1][0], tp_utr[0][1], strand, 0))
        else:
            feature_lst.append((tp_utr[0][0], exons[-1][1], strand, 0))
        feature_lst.append((exons[-1][0], exons[-1][1], strand, 1))
        feature_lst += tp_utr[1:]

    return feature_lst


def update_ln(st_pos, en_pos, fields):
    assert len(fields) == 9
    ln = ""
    for i in range(len(fields)):
        if i == 3:
            ln += str(st_pos) + "\t"
        elif i == 4:
            ln += str(en_pos) + "\t"
        elif i == len(fields) - 1:
            ln += fields[i] + "\n"
        else:
            ln += fields[i] + "\t"
    # print(ln)
    return ln


def build_ln(st_pos, en_pos, ft_type, fields):
    if ft_type == 0:
        ft_type_s = "exon"
    elif ft_type == 1:
        ft_type_s = "CDS"
    else:
        print("impossible")
        sys.exit()
    assert len(fields) == 9
    ln = ""
    for i in range(len(fields)):
        if i == 2:
            ln += ft_type_s + "\t"
        elif i == 3:
            ln += str(st_pos) + "\t"
        elif i == 4:
            ln += str(en_pos) + "\t"
        elif i == len(fields) - 1:
            ln += fields[i] + "\n"
        else:
            ln += fields[i] + "\t"
    # print(ln)
    return ln


def fix_rec(feature_tbl, fn, out_fn):
    fh = open(fn, 'r')
    out_fh = open(out_fn, 'w')
    for ln in fh:
        clean_ln = ln.strip()
        fields = clean_ln.split("\t")
        if fields[2] == "transcript":
            ps = get_pairs(fields[8])
            strand = fields[6]
            for k, v in ps:
                if k == "transcript_id":
                    tid = v
                    break
            features = feature_tbl[tid]
            if strand == '+':
                tx_coords = (features[0][0], features[-1][1])
            else:
                tx_coords = (features[-1][0], features[0][1])
            tx_ln = update_ln(tx_coords[0], tx_coords[1], fields)
            out_fh.write(tx_ln)
            prev_ft = None
            for ft in features:
                st_pos, en_pos, _, ft_type = ft
                if prev_ft == ft_type:
                    if ft_type == 1:
                        ft_ln = build_ln(st_pos, en_pos, 0, fields)
                        out_fh.write(ft_ln)
                ft_ln = build_ln(st_pos, en_pos, ft_type, fields)
                out_fh.write(ft_ln)
                prev_ft = ft_type
        else:
            continue
    fh.close()
    out_fh.close()


def main(ref_fn, qry_fn, out_fn):
    utr_tbl = load_utr(ref_fn)
    q2r_tbl, exon_tbl = load_qry(qry_fn)
    add_utr_tbl = add_utr(utr_tbl, q2r_tbl, exon_tbl)
    feature_tbl = dict()
    for tid in exon_tbl:
        exon_lst = exon_tbl[tid]
        _, _, strand = exon_lst[0]
        add_fp_utr, add_tp_utr = add_utr_tbl[tid]
        feature_lst = add_features(strand, exon_lst, add_fp_utr, add_tp_utr)
        feature_tbl[tid] = feature_lst
    fix_rec(feature_tbl, qry_fn, out_fn)


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2], sys.argv[3])







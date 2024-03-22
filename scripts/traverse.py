#!/usr/bin/env python

import sys
import os


def load_graph(u2dr_fn, dr2u_fn):
    u2dr_fh = open(u2dr_fn, 'r')
    dr2u_fh = open(dr2u_fn, 'r')
    tx2ref_tbl = dict()
    ref2tx_tbl = dict()
    for line in u2dr_fh:
        line = line.strip()
        fields = line.split("\t")
        for i in range(len(fields)):
            field = fields[i]
            if i == 0:
                tid = field
                tx2ref_tbl[tid] = list()
            else:
                tx2ref_tbl[tid].append(field)
    for line in dr2u_fh:
        line = line.strip()
        fields = line.split("\t")
        rid = fields[0]
        tid = fields[1]
        if rid not in ref2tx_tbl:
            ref2tx_tbl[rid] = [tid]
        else:
            ref2tx_tbl[rid].append(tid)
    return tx2ref_tbl, ref2tx_tbl


# no cycle
def find_island(tx2ref_tbl, ref2tx_tbl, src):
    island = []
    queue = [(src, True)]
    to_del = set()
    visited = set()
    queue_el = [src]
    while len(queue) > 0:
        curr, is_tx = queue.pop(0)
        queue_el.pop(0)
        island.append((curr, is_tx))
        visited.add(curr)
        if is_tx:
            to_del.add(curr)
            is_tx = False
            nxt_lvl = tx2ref_tbl[curr]
            nxt_lvl = [x for x in nxt_lvl if x not in visited and x not in queue_el]
        else:
            is_tx = True
            nxt_lvl = ref2tx_tbl[curr]
            nxt_lvl = [x for x in nxt_lvl if x not in visited and x not in queue_el]
        for node in nxt_lvl:
            queue.append((node, is_tx))
            queue_el.append(node)
    try:
        assert len(island) == len(set(island))
    except:
        print(island)
    return to_del, island


def traverse(tx2ref_tbl, ref2tx_tbl):
    src_lst = list(tx2ref_tbl)
    src = src_lst[0]
    island_map = dict()
    while True:
        to_del, island = find_island(tx2ref_tbl, ref2tx_tbl, src)
        if src in island_map:
            print("something is wrong")
            sys.exit()
        island_map[src] = island
        src_lst = [x for x in src_lst if x not in to_del]
        if len(src_lst) == 0:
            break
        src = src_lst[0]
    return island_map


def split(island_map, n):
    # file name counter
    ctr = 0
    island_fn = "island." + str(ctr) + ".txt"
    island_fh = open(island_fn, 'w')
    # number of elements in the current archipelago
    num_el = 0
    archipelago = list()
    arch_lst = list()
    rid_set = set()
    tid_set = set()
    for src in island_map:
        src_island = island_map[src]
        tmp = num_el + len(src_island)
        if tmp > n:
            for island in archipelago:
                for el, is_tx in island:
                    if is_tx:
                        tid_set.add(el)
                    else:
                        rid_set.add(el)
                    island_fh.write(el + "\t" + str(is_tx) + "\n")
            island_fh.close()
            arch_lst.append(archipelago)
            archipelago = list()
            ctr += 1
            island_fn = "island." + str(ctr) + ".txt"
            island_fh = open(island_fn, 'w')
            num_el = len(src_island)
            archipelago.append(src_island)
        else:
            archipelago.append(src_island)
            num_el += len(src_island)

    if len(archipelago) != 0:
        for island in archipelago:
            for el, is_tx in island:
                if is_tx:
                    tid_set.add(el)
                else:
                    rid_set.add(el)
                island_fh.write(el + "\t" + str(is_tx) + "\n")
        island_fh.close()
        arch_lst.append(archipelago)

    return arch_lst


def main(n, u2dr_fn, dr2u_fn):
    tx2ref_tbl, ref2tx_tbl = load_graph(u2dr_fn, dr2u_fn)
    # glob_rid_set = set(ref2tx_tbl.keys())
    # glob_tid_set = set(tx2ref_tbl.keys())
    island_map = traverse(tx2ref_tbl, ref2tx_tbl)
    arch_lst = split(island_map, n)
    print("total number of archipelagos: " + len(arch_lst))


if __name__ == "__main__":
    main(int(sys.argv[1]), sys.argv[2], sys.argv[3])

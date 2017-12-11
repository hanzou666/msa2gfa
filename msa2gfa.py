#!/usr/bin/env python
# -*- coding: utf-8 -*-
# author: Shuho Ohwada
# date: 2017/12/11

"""
Extract graph structure from multiple alignment result and output as GFA/JSON for vg

Usage:
msa2gfa.py -f msa.fa > graph.gfa (for a fasta file)
msa2gfa.py -l path_list.txt > graph.gfa (for multiple fasta files)

It supports for vg v1.5 and not <= v1.4.
"""

import argparse
import json
import re
import sys


def get_argument():
    """Parse argument"""
    parser = argparse.ArgumentParser(
        description='Extract graph structure from multiple alignment result\
                     and output as GFA/JSON')
    parser.add_argument('-f', '--fasta',
                        help='An input fasta file')
    parser.add_argument('-l', '--list',
                        help='A file with the paths of multiple fasta files')
    parser.add_argument('-j', '--json',
                        action='store_true',
                        help='Output as JSON format instead of GFA format')
    args = parser.parse_args()
    return args


def extract_graph(fasta, first_id):
    """Extract graph structure from input fasta file with msa

    input: txtfile(fasta format), int(first node id)
    output: dict(node, edge, path), int(next node id)

    Method:
    1. Parse a fasta file
    2. Get node per base with eliminating duplication in the same coordinate
    3. Merge nodes in a non-branch area in one graph into one node
    4. Add edge using node path information
    5. Transform node and edge to vg like graph
    """
    fasta_dic = parse_fasta(fasta)
    base_node_dic = get_node(fasta_dic)
    merged_node_dic, next_id = merge_nodes(base_node_dic, fasta_dic.keys(), first_id)
    new_edge_dic = add_edge(merged_node_dic, fasta_dic.keys())
    vg_like_graph = transform_to_vg_like_graph(
        merged_node_dic, new_edge_dic, fasta_dic.keys())
    return vg_like_graph, next_id


def parse_fasta(fasta):
    """Parse fasta file as dictionary"""
    fasta_dic = {}
    for tmpline in open(fasta, 'r'):
        if tmpline[0] == '>':
            seqtitle = re.match(r"^\>\s*(\S*)\s", tmpline)
            seq_name = seqtitle.group(1)
            fasta_dic[seq_name] = ''
        else:
            fasta_dic[seq_name] += tmpline[:-1]
    return fasta_dic


def get_node(fasta_dic):
    """Get nodes with base and path name from multiple sequences"""
    def get_coordinate_info(fasta_dic):
        """Divide the sequences into one base units"""
        coordinate_info_list = []
        for seqid, seq in fasta_dic.items():
            tmp_list = [(seqid, base) for base in seq]
            coordinate_info_list.append(tmp_list)
        return coordinate_info_list

    coordinate_info_list = get_coordinate_info(fasta_dic)
    node_list = []
    for tmp_coordinate_info in zip(*coordinate_info_list):
        base_set = set([i[1] for i in tmp_coordinate_info])
        if len(base_set) == 1:
            node_list.append({'base': list(base_set)[0],
                              'seq_name': [i[0] for i in tmp_coordinate_info]})
        else:
            for base in base_set:
                if base != '-':
                    node_list.append({'base': base,
                                      'seq_name': [i[0] for i in tmp_coordinate_info
                                                   if i[1] == base]})
    return {i: j for i, j in enumerate(node_list, 1)}


def merge_nodes(node_dic, seq_name_list, first_id):
    """Merge nodes in a non-branch area in one graph into one node"""
    tmp_edge_dic = add_edge(node_dic, seq_name_list)
    single_link_dic = {}
    for start_id, end_ids in tmp_edge_dic.items():
        end_id_list = [i for i in end_ids if i > start_id]
        if len(end_id_list) == 1:
            end_id = end_id_list[0]
            start_id_list = [i for i in tmp_edge_dic[end_id] if i < end_id]
            if len(start_id_list) == 1:
                single_link_dic[end_id] = start_id
    end_id_list = [i for i in sorted(single_link_dic.keys(), reverse=True)]
    for end_id in end_id_list:
        if end_id in single_link_dic.keys():
            start_id = single_link_dic[end_id]
            node_dic[start_id]['base'] += node_dic[end_id]['base']
            del node_dic[end_id]
    merged_node_dic = {i: node_dic[j]
                       for i, j in enumerate(sorted(node_dic.keys()), first_id)}
    next_id = max(merged_node_dic.keys()) + 1
    return merged_node_dic, next_id


def add_edge(node_dic, seq_name_list):
    """Connect nodes as undirected graph with tracing the original paths"""
    def get_first_node_id(node_dic, seq_name):
        """Find first node ID in the path"""
        for i, node_info in sorted(node_dic.items()):
            if seq_name in node_info['seq_name']:
                return i

    edge_dic = {i: set() for i in node_dic.keys()}
    for seq_name in seq_name_list:
        start_id = get_first_node_id(node_dic, seq_name)
        tmp_id = start_id + 1
        while tmp_id <= len(node_dic.keys()):
            if seq_name in node_dic[tmp_id]['seq_name']:
                edge_dic[tmp_id].add(start_id)
                edge_dic[start_id].add(tmp_id)
                start_id = tmp_id
            tmp_id += 1
    return edge_dic


def transform_to_vg_like_graph(node_dic, edge_dic, seq_name_list):
    """Shape graph structure as vg like graph dict"""
    vg_like_graph = {'node': [], 'edge': [], 'path': []}
    vg_like_graph['node'] = [{'name': str(i),
                              'sequence': j['base'],
                              'id': i}
                             for i, j in sorted(node_dic.items())]
    for start_id, end_ids in sorted(edge_dic.items()):
        vg_like_graph['edge'] += [{'from': start_id, 'to': end_id}
                                  for end_id in sorted(end_ids)
                                  if end_id > start_id]
    for seq_name in seq_name_list:
        tmp_mapping_list = []
        i = 0
        for node_id in sorted(node_dic.keys()):
            if seq_name in node_dic[node_id]['seq_name']:
                i += 1
                tmp_mapping_list.append({'position': {'node_id': node_id}, 'rank': i})
        vg_like_graph['path'].append(
            {'name': seq_name, 'mapping': tmp_mapping_list})
    return vg_like_graph


def output_as_json(vg_like_graph):
    """Convert vg_like_graph to JSON and output as stdout"""
    sys.stdout.write(json.dumps(vg_like_graph))
    sys.stdout.write('\n')


def output_as_gfa(vg_like_graph):
    """Convert vg_like_graph to GFA 1.0 and output as stdout

    It supports only vg v1.5.
    """
    sys.stdout.write('H\tVN:Z:1.0\n')
    for tmp_path in vg_like_graph['path']:
        seq_name = tmp_path['name']
        mapping_list = tmp_path['mapping']
        node_len_dic = {i['id']: len(i['sequence']) for i in vg_like_graph['node']}
        tmp_path_list = []
        tmp_len_list = []
        for mapping in mapping_list:
            node_id = mapping['position']['node_id']
            tmp_path_list.append('{}+'.format(node_id))
            tmp_len_list.append('{}M'.format(node_len_dic[node_id]))
        path_csv = ','.join(tmp_path_list)
        len_csv = ','.join(tmp_len_list)
        sys.stdout.write('P\t{}\t{}\t{}\n'.format(seq_name, path_csv, len_csv))
    for node in vg_like_graph['node']:
        sys.stdout.write('S\t{}\t{}\n'.format(node['id'], node['sequence']))
    for edge in vg_like_graph["edge"]:
        sys.stdout.write('L\t{}\t+\t{}\t+\t0M\n'.format(edge['from'], edge['to']))


def main():
    """Run"""
    args = get_argument()

    if args.fasta:
        vg_like_graph, tmp_counter = extract_graph(args.fasta, 1)
    elif args.list:
        first_id = 1
        vg_like_graph = {'node': [], 'edge': [], 'path': []}
        for fasta_path in args.list:
            tmp_vg_like_graph, next_id = extract_graph(fasta_path, first_id)
            first_id = next_id
            for keyname in vg_like_graph:
                vg_like_graph[keyname].append(tmp_vg_like_graph[keyname])
    else:
        sys.stderr.write(__doc__)
        sys.exit(1)

    if args.json:
        output_as_json(vg_like_graph)
    else:
        output_as_gfa(vg_like_graph)


if __name__ == '__main__':
    main()

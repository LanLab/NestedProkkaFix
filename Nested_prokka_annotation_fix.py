from time import sleep as sl
import sys
import os
import math
import pandas as pd
import numpy as np
import argparse
import glob
#import matplotlib
#matplotlib.use('pdf')
#import matplotlib.pyplot as plt
import progressbar
import json
from io import StringIO
import csv

#1 Filtering ortho groups for presence and cleaning names.
def pre_filtering(args):
    print('Filtering and cleaning Roary input.')
    # function will isolate genes found in 99% of strains for each ortholog group
    df = pd.read_csv(args.roary_csv, low_memory=False, index_col=False)

    strains_used = df.shape[1] - 14
    df.insert(loc = 6, column = 'Gene Count', value = df.iloc[:,14:].notnull().sum(axis=1))
    df = df.loc[df['Gene Count'] >= math.ceil(strains_used * args.presence_percen)]
    working_df = pd.DataFrame
    working_df = df.loc[df['Gene Count'] >= math.ceil(strains_used * args.presence_percen)]
    working_df.to_csv(args.output + '/filtered_gap.csv',sep=',', index= False)

    return working_df
def clean_naming(temp_file, args):
    # converting temp_file dataframe into list of lists for pythonic str functions etc
    temp_file = [temp_file.columns.values.tolist()] + temp_file.values.tolist()
    for line_num in range(1, len(temp_file), 1):
        for cell_num in range(15, len(temp_file[line_num]), 1):
            #fix .fasta in cells
            try:
                if '.fasta' in temp_file[line_num][cell_num]:
                    temp_file[line_num][cell_num] = temp_file[line_num][cell_num].replace('.fasta','')
            except:
                pass
            try:
                #fix cds in cells
                if 'cds-' in temp_file[line_num][cell_num]:
                    temp_file[line_num][cell_num] = temp_file[line_num][cell_num].replace('cds','')
            except:
                pass
            try:
            # fix cds in cells
                if '__' in temp_file[line_num][cell_num]:
                    temp_file[line_num][cell_num] = temp_file[line_num][cell_num].replace('__', '_')
            except:
                pass
    return temp_file

#2 Read in previoulsy created dictionary or create genome info dict.
def parsingGFFs(args):
    ###will read in previoulsy created dictionary or create genome info dict if needed
    if args.gff_dictionary:
        print('Reading in dictionary with all gff genome information.')
        with open(args.gff_dictionary) as dict_file:
            info_dict = json.loads(dict_file.read())
            return info_dict

    elif args.input_gff_folder:
        print( 'Creating dictionary with all gff genome information.')
        info_dict = creating_genome_info(args)
        return info_dict
def creating_genome_info(args):
    number_of_genomes = len(list(glob.iglob(args.input_gff_folder + '/*.gff'))) + 1
    bar = progressbar.ProgressBar(maxval=number_of_genomes, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    info_dict = {}
    genome_count = 0
    for filename in glob.iglob(args.input_gff_folder + '/*.gff'):
        bar.update(genome_count)
        genome_count = genome_count + 1
        with open(filename, 'r') as file:
            strain_name = filename.split('/')[-1].rstrip('.fasta.gff')
            strain_name = strain_name.replace('PROKKA_', '')
            strain_name =strain_name.split('.')[0]
            info_dict[strain_name] = {}
            for line in file:
                if ('CDS' in line and 'ID=' in line) or ('RNA' in line and 'ID=' in line):
                    col = line.split('\t')
                    beg_cord = col[3]
                    end_cord = col[4]
                    strand = col[6]
                    size = int(end_cord) - int(beg_cord)
                    gene_id = col[8].split(';')[0].split('=')[-1]
                    # if gene_id[-2] == '.':
                    #     gene_id = gene_id.split('.')[0]
                    for q in col[8].split(';'):
                        if "locus_tag" in q:
                            genbank = q.split('=')[-1]
                    if '_' not in col[0]:
                        contig = col[0]
                    elif '_' in col[0]:
                        if len(col[0].split('_')) > 2:
                            contig = col[0].split('_')[0] + '_' + col[0].split('_')[1]
                        else:
                            contig = col[0]
                    info_dict[strain_name][gene_id] = {'size': size, 'start': beg_cord, 'end': end_cord, 'contig': contig, 'locus_tag':genbank, 'strand':strand}
    bar.finish()
    #writing out the genome dictionary to file to save time
    print('\n')
    with open(args.output + '/GFF_annotation.json', 'w') as dict_file:
        dict_file.write(json.dumps(info_dict))
    args.gff_dictionary = args.output + '/GFF_annotation.json'
    return info_dict

#3 Gathering information core gene groups
def gathering_core_gene_information(temp_file, info_dict, args):

    ###script will go over earch ortholog, first find average size of non paralogs, then find combined sizes of each paralog.
    print('Gathering core gene information.')
    sl(1)
    percen_len = sum(1 for row in temp_file)
    bar = progressbar.ProgressBar(maxval=percen_len, widgets=[progressbar.Bar('=', '[', ']'), ' ', progressbar.Percentage()])
    bar.start()
    group_count = 1
    for line in range(1, len(temp_file)):
    # for line in range(1, 6):
        ##creating progress var
        bar.update(group_count)
        group_count = group_count + 1
        # for each strain for the core gene
        for j in range(15, len(temp_file[line])):
            parseAcc = temp_file[0][j].replace('.fasta', '').replace('PROKKA_', '').split('.')[0]
            # print(parseAcc)
            info = temp_file[line][j]
            # will leave empty strains
            if 'nan' in str(info):
                pass
            # handling cells which DONT have orthologs
            elif '\t' not in str(info):
                acc, contig, gene_id, beg_cord, end_cord, size, strand = handling_roary_single_annoations_strings(str(info), info_dict, parseAcc)
                temp_file[line][j] = (acc + '_' + str(contig) + '_' + str(gene_id) + '_' + str(beg_cord) + '_' + str(end_cord) + '_' + str(size) + '(' + str(strand) + ')')
            # handling cells which have paralogs
            elif '\t' in str(info):
                frag_cell, frag_length_list = handling_roary_paralogs_annotations_strings(str(info), info_dict, parseAcc)
                temp_file[line][j] = frag_cell

    with open(args.output + '/gap_with_details.csv', 'w') as out:
        for line in temp_file:
            for j in line:
                out.write(str(j) + ',')
            out.write('\n')

    return temp_file
def load_roary_with_detais(args):
    infile = open(args.roary_with_details).read().splitlines()
    temp_file = []
    for i in infile:
        tmp_line = i.split(',')
        temp_file.append(tmp_line)
    return temp_file

#4 Handling Roary fasely split orthologs and fasely joined paralogs
def process_roary_input(temp_file, args):
    print('Handling erroneous core genes.')
    core_gene_list = []
    excluded_list = []
    split_count = 0
    ##for each line of the Roary input
    print( 'Identifying fasely split orthologs and fasely joined paralogs.')
    for line in range(1, len(temp_file),1):
        keep_core_gene = True
        split_core_gene = False
        # print(temp_file[line])
        ##handle fasely split orthologs
        keep_core_gene, reference_core_gene_fso = fasely_split_ortholgs_filter(keep_core_gene, temp_file, line, args)
        ##handle fasely joined orthologs
        split_core_gene, reference_core_gene_fjp = fasely_joined_orthologs_filter(split_core_gene, temp_file, line, args)
        #add line to new_temp_file
        if keep_core_gene == True and split_core_gene == False:
            core_gene_list.append(reference_core_gene_fso)
        ##spliting fasely joined cells and appending
        if split_core_gene == True:
            split_count = split_count + 1
            single_gene = reference_core_gene_fjp.split('\t')
            for j in single_gene:
                core_gene_list.append(j)
        if keep_core_gene == False and split_core_gene == False:
            excluded_list.append(temp_file[line])
        # print(keep_core_gene, split_core_gene)
    return core_gene_list, excluded_list
def fasely_split_ortholgs_filter(keep_core_gene, temp_file, line, args):
    reference_core_gene = ""
    cell_in_size = 0
    true_para = 0
    false_para = 0
    for j in range(15, len(temp_file[line]), 1):
        j = str(temp_file[line][j])
        #find reference gene size if not fragmented
        if (args.reference in j) and ('\t' not in j) and ('nan' not in j):
            reference_core_gene = j
            reference_size = int(j.split('_')[-1].split('(')[0])
            # find size of other cells in the same core group
            for l in temp_file[line][15:]:
                l = str(l)
                # if other cells are a single fragment
                if '\t' not in l and 'nan' not in l and l != '':
                    other_size = int(l.split('_')[-1].split('(')[0])
                    # compare size of reference cell against other cells
                    ref_lower_range = int(0.8 * reference_size)
                    ref_upper_range = int(1.2 * reference_size)
                    if ref_lower_range <= other_size <= ref_upper_range:
                        cell_in_size = cell_in_size + 1
                    # else:
                    #     print('single', l, other_size, reference_core_gene, reference_size, ref_lower_range, ref_upper_range)
                # if other cells are fragmented
                if '\t' in l and 'nan' not in l:
                    frag_sizes_list = []
                    fragments = l.split('\t')
                    for frag in fragments:
                        frag_size = int(frag.split('_')[-1].split('(')[0])
                        frag_sizes_list.append(frag_size)
                    total_frag_size = sum(frag_sizes_list)
                    # compare size of reference cell against other cells
                    ref_lower_range = int(0.8 * reference_size)
                    ref_upper_range = int(1.2 * reference_size)
                    if ref_lower_range <= total_frag_size <= ref_upper_range:
                        cell_in_size = cell_in_size + 1
                    # else:
                    #     print('multi', fragments, total_frag_size, reference_core_gene, reference_size, ref_lower_range, ref_upper_range)

    if float(args.presence_percen) < 1.0:
        if cell_in_size <= (args.presence_percen * len(temp_file[line][15:])):
            keep_core_gene = False
    else:
        if cell_in_size != (args.presence_percen * len(temp_file[line][15:])):
            keep_core_gene = False

    return keep_core_gene, reference_core_gene
def fasely_joined_orthologs_filter(split_core_gene, temp_file, line, args):
    reference_core_gene = ""
    for j in range(15, len(temp_file[line]), 1):
        j = str(temp_file[line][j])
        if args.reference in j and '\t' in j and 'nan' not in j:
            reference_core_gene = j
            ##find combined length of fragmented cells
            frag_sizes_list = []
            fragments = j.split('\t')
            for frag in fragments:
                frag_size = int(frag.split('_')[-1].split('(')[0])
                frag_sizes_list.append(frag_size)
            total_frag_size = sum(frag_sizes_list)
            ##check if single gene is in line and if that gene is similar size to combined fragments
            single_gene_present = False
            for q in temp_file[line][15:]:
                q = str(q)
                if '\t' not in q and 'nan' not in q and q != '':
                    single_gene_size = int(q.split('_')[-1].split('(')[0])
                    lower_limit_ref = int(0.8 * total_frag_size)
                    upper_limit_ref = int(1.2 * total_frag_size)
                    if lower_limit_ref <= single_gene_size <= upper_limit_ref:
                        single_gene_present = True
            # if single gene present, then treat core gene as fasely joined paralogs
            if single_gene_present == True:
                fragments = j.split('\t')
                # comparing the sizes of fasely joined paralogs to sizes in other genome cells
                # taking reference fragment size
                both_frags_in_size = 0
                for frag in fragments:
                    frag_in_size = 0
                    frag_size = int(frag.split('_')[-1].split('(')[0])
                    # taking other cells fragment size
                    for other_cells in temp_file[line][15:]:
                        other_cells = str(other_cells)
                        other_cells_fragments = other_cells.split('\t')
                        for other_cells_frags in other_cells_fragments:
                            if 'nan' not in other_cells_frags and other_cells_frags != '':
                                other_cell_size = int(other_cells_frags.split('_')[-1].split('(')[0])
                                # creating upper and lower limits for variation in size between ref and other cell fragments
                                lower_limit_ref = int(0.8 * frag_size)
                                upper_limit_ref = int(1.2 * frag_size)
                                # comparing each reference fragment against the size of other cell fragments
                                if lower_limit_ref <= other_cell_size <= upper_limit_ref:
                                    frag_in_size = frag_in_size + 1
                        # checking if the fragment is found across >= 99% of the test genomes

                    if float(args.presence_percen) < 1.0:
                        if frag_in_size >= (args.presence_percen * len(temp_file[line][15:])):
                            both_frags_in_size = both_frags_in_size + 1
                    else:
                        if frag_in_size != (args.presence_percen * len(temp_file[line][15:])):
                            both_frags_in_size = both_frags_in_size + 1

                    # if frag_in_size >= (args.presence_percen * len(temp_file[line][15:])):
                    #     both_frags_in_size = both_frags_in_size + 1
                # check if should split core genes
                if both_frags_in_size == len(fragments):
                    split_core_gene = True
    return split_core_gene, reference_core_gene
def handling_roary_single_annoations_strings(j, info_dict, acc):
    #will handle all the different roary annotations to extract gene_id, start, end and size
    gene_id = j.rstrip('"')
    if gene_id not in info_dict[acc]:
        gene_id = "cds" + gene_id
    contig = info_dict[acc][gene_id]['contig']
    beg_cord = info_dict[acc][gene_id]['start']
    end_cord = info_dict[acc][gene_id]['end']
    size = info_dict[acc][gene_id]['size']
    strand = info_dict[acc][gene_id]['strand']
    return acc, contig, gene_id, beg_cord, end_cord, size, strand
def handling_roary_paralogs_annotations_strings(j, info_dict,acc):
    frag_length_list = []
    frag_count = 0
    para_count = j.count('\t') + 1
    frag_list = []
    fragments = j.split('\t')
    for frag in range(para_count):
        gene_id = fragments[frag].rstrip('"')
        entry = info_dict[acc]
        if gene_id not in info_dict[acc]:
            gene_id = "cds"+gene_id

        contig = info_dict[acc][gene_id]['contig']
        beg_cord = int(info_dict[acc][gene_id]['start'])
        end_cord = int(info_dict[acc][gene_id]['end'])
        size = info_dict[acc][gene_id]['size']
        strand = info_dict[acc][gene_id]['strand']
        frag_list.append(acc + '_' + contig + '_' + gene_id + '_' + str(beg_cord) + '_' + str(end_cord) + '_' + str(size) + '(' + str(strand) + ')')
        frag_length_list.append((beg_cord, end_cord))
    frag_cell = '\t'.join(frag_list)
    return frag_cell, frag_length_list

#5 Writing out core genes
def core_gene_out(core_gene_list, excluded_list, args):
    print( 'Writing out list of core genes to: ' + args.output + '/core_genes.csv')
    gene_count = 0
    with open(args.output + '/core_genes.csv','w') as outfile:
        outfile.write(','.join(["Accession", "CDS", "Locus Tag", "Start", "End", "Length"]))
        outfile.write('\n')
        for line in core_gene_list:
            col = line.split('_')
            acc = args.reference
            cds = col[-4]
            start = col[-3]
            end = col[-2]
            length = col[-1]
            outfile.write(','.join([acc,cds,start,end,length]))
            outfile.write('\n')
            gene_count += 1
    with open(args.output + '/excluded_list.csv','w') as outfile:
        for line in excluded_list:
            for cell in line:
                outfile.write(str(cell) + ',')
            outfile.write('\n')
    print('Total Core Genes:' + str(gene_count))

#aux
def merge(intervals):
   if not intervals:
       return []
   data = []
   for interval in intervals:
       data.append((interval[0], 0))
       data.append((interval[1], 1))
   data.sort()
   merged = []
   stack = [data[0]]
   for i in range(1, len(data)):
       d = data[i]
       if d[1] == 0:
           # this is a lower bound, push this onto the stack
           stack.append(d)
       elif d[1] == 1:
           if stack:
               start = stack.pop()
           if len(stack) == 0:
               # we have found our merged interval
               merged.append((start[0], d[0]))
   return merged

#main
def parseargs():
    parser = argparse.ArgumentParser(formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("-o", "--output", help="Output direcotry.", required=True)
    parser.add_argument("-p", "--presence_percen", help="Percentage of isolates a gene must be present.", default=0.95, type=float)
    parser.add_argument("-r", "--reference", help="Reference accession.", required=True)

    roary_parser_group = parser.add_mutually_exclusive_group(required=True)
    roary_parser_group.add_argument("-i", "--roary_csv", help="Roary presence and absence matrix.")
    roary_parser_group.add_argument("-id", "--roary_with_details", help="Roary presence and absence matrix with annotationd details added.")

    gff_parser_group = parser.add_mutually_exclusive_group(required=True)
    gff_parser_group.add_argument("-gff", "--input_gff_folder", help="path to folder containing all the prokka gff file")
    gff_parser_group.add_argument("-dict", "--gff_dictionary", help="Preprocessed .JSON of GFF annotations.")

    args = parser.parse_args()



    return args
def main():
    args = parseargs()

    if args.roary_csv:
        temp_file = pre_filtering(args)
        temp_file = clean_naming(temp_file, args)
        info_dict = parsingGFFs(args)
        temp_file = gathering_core_gene_information(temp_file, info_dict, args)
        core_gene_list, excluded_list = process_roary_input(temp_file, args)
        core_gene_out(core_gene_list, excluded_list, args)

    elif args.roary_with_details:
        print(args.roary_with_details)
        temp_file = load_roary_with_detais(args)
        core_gene_list, excluded_list = process_roary_input(temp_file, args)
        core_gene_out(core_gene_list, excluded_list, args)

if __name__ == '__main__':
    main()














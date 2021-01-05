#!/usr/bin/env python3
# -*- coding: utf-8 -*-

'''This script is able to align reads to a database with genomes of reference.
The result obtained includes a fasta file  with the generated scaffolds,
coverage images for the mapped references, non-redundant sequence file
and a summary of what was found'''

from PIL import Image
from PIL import ImageTk
from tkinter import ttk
from statistics import mean
from operator import itemgetter
from pandas.plotting import table
from tkinter import Tk, INSIDE, filedialog, Toplevel, IntVar, messagebox, DISABLED, NORMAL, SUNKEN, CENTER
import os
import re
import time
import gzip
import math
import queue
import threading
import subprocess
import pandas as pd
import tkinter as tk
import seaborn as sns
import matplotlib.pyplot as plt
__author__ = "Ary Rivillas"
__copyright__ = "Copyright 2020, Biotecnologia Microbiana Research Group, Universidad Nacional de Colombia - Sede Medellín"
__credits__ = [
    "Pablo Gutierrez",
    "Mauricio Marin",
    "Yuliana Gallo",
    "Daniel Tejada",
    "Andrea Restrepo",
    "Susana Giraldo"]
__license__ = "MIT"
__version__ = "1.0"
__maintainer__ = "Ary Rivillas"
__email__ = "amrivillast@unal.edu.co"
__status__ = "Development"
__date__ = '2020/06/19'


'''Required libraries Pillow, Seaborn, matplotlib, pandas'''


start_time = time.time()
dir_path = os.path.dirname(os.path.realpath(__file__))
result_path = "Results"
db_path = "Databases"
script_path = "Scripts"
db_name = "viral_refseqs"
db_name_tsv = "database_viral.tsv"
base_name_mgb_out = "magicblast_out.sam"
base_name_short_sam = "short_sam.tsv"
base_name_temp_query = "tempQuery.fasta"
base_name_temp_subject = "tempSubject.fasta"
base_name_temp_result_blastn = "tempOut.txt"
path_to_results = os.path.join(dir_path, result_path)
path_to_results_magicblast = os.path.join(path_to_results, base_name_mgb_out)
path_to_results_short_sam = os.path.join(path_to_results, base_name_short_sam)
path_to_db_tsv = os.path.join(dir_path, db_path, db_name_tsv)
path_to_db = os.path.join(dir_path, db_path, db_name)
path_to_scripts = os.path.join(dir_path, script_path)
path_to_temp_query = os.path.join(path_to_results, base_name_temp_query)
path_to_temp_subject = os.path.join(path_to_results, base_name_temp_subject)
path_to_temp_query_result_blastn = os.path.join(
    path_to_results, base_name_temp_result_blastn)
check_result_folder = os.path.isdir(path_to_results)

if not check_result_folder:

    os.mkdir(path_to_results)

list_of_images = []
list_of_graphs_to_pdf = []
checker_time = []
global flag_sam_proc
flag_sam_proc = 0
global check_mb_proc
check_mb_proc = 0
global limit
limit = 0.1
global pros
pros = 0
global ver
ver = 0


def main():

    root = Tk()
    VirusMapperUI(root)
    root.mainloop()


def rna_seq_aligner(self, name_file1):
    '''Align with Magic-BLAST remotely'''

    global ver
    name_file_to_process = self.name_main + "_nr.fasta"
    path_to_nr_file = os.path.join(self.path_to_sr, name_file_to_process)
    command1 = 'magicblast -db ' + path_to_db + ' -query ' + path_to_nr_file + \
        ' -reftype transcriptome -infmt fasta -outfmt sam -no_unaligned -out ' + \
        path_to_results_magicblast
    subprocess.Popen(command1, shell=False)
    time.sleep(5)
    ver = 1


def non_red_seqs(self):
    '''Generate a non-redundant sequences file'''

    def complement_fnc(x):
        '''Complementarity'''

        answer = []
        base = {'A': 'T',
                'T': 'A',
                'G': 'C',
                'C': 'G',
                'N': 'N'}

        for nucleotide in x:

            answer.append(base[nucleotide])

        x_rev = ''.join(answer[::])

        return x_rev

    def convert_fnc(tuple, di):

        for a, b in tuple:

            di.setdefault(a, []).append(b)

        return di

    self.statusbar.config(text="Generating non-redundant sequences")
    print("\nGenerating non-redundant sequences...\n")

    count = 0
    r_name = []
    r_seq = []
    r_seqs_com = []
    r_seqs_com2 = []
    ph_t = 0
    update_seq = ""
    seq_name_dict = {}
    seq_mult_dict = {}
    seq_ph_dict = {}
    seq_ph_a_dict = {}

    sequence_file1 = self.route_file1
    sequence_file2 = self.route_file2
    name = sequence_file1[::-1].split("/")[0][::-1]

    name_ext = sequence_file1[::-1].split(".")[0][::-1]
    name_ext_checker = ["gz"]

    self.root.update_idletasks()

    if name_ext in name_ext_checker:

        if sequence_file2 != "":

            filenames = [sequence_file1, sequence_file2]

            for fname in filenames:

                with gzip.open(fname, "rt") as infile:

                    for line in infile:

                        count += 1
                        line_mod = line.strip()
                        self.root.update()

                        if count == 1:

                            r_name.append(line_mod)
                            r_name_dict = line_mod

                        if count == 2:

                            r_seq.append(line_mod)
                            r_comp = complement_fnc(line_mod)
                            r_seqs_com.append(line_mod)
                            r_seqs_com.append(r_comp)
                            r_seqs_com.sort()

                            if r_seqs_com[0] not in seq_name_dict:

                                update_seq = r_seqs_com[0]
                                number = 1
                                new_element2 = {update_seq: number}
                                seq_mult_dict.update(new_element2)

                            if r_seqs_com[0] in seq_name_dict:

                                number = seq_mult_dict.get(r_seqs_com[0])
                                number += 1
                                update_seq = r_seqs_com[0]
                                new_element3 = {update_seq: number}
                                seq_mult_dict.update(new_element3)

                        new_element = {update_seq: r_name_dict}

                        if update_seq != "":

                            seq_name_dict.update(new_element)

                        if count == 4:

                            for x in line_mod:

                                phred = int(ord(x) - 33)
                                ph_t = ph_t + int(phred)

                            ph_prom = round(ph_t / int(len(line_mod)))
                            r_comp2 = complement_fnc(r_seq[0])
                            r_seqs_com2.append(r_seq[0])
                            r_seqs_com2.append(r_comp2)
                            r_seqs_com2.sort()

                            if r_seqs_com2[0] not in seq_ph_dict:

                                new_element4 = {r_seqs_com2[0]: ph_prom}
                                seq_ph_dict.update(new_element4)

                            n = seq_mult_dict.get(r_seqs_com2[0])

                            if r_seqs_com2[0] in seq_ph_dict and n > 1:

                                number_2 = seq_ph_dict.get(r_seqs_com[0])
                                number_2 += ph_prom
                                new_element6 = {r_seqs_com2[0]: number_2}
                                seq_ph_dict.update(new_element6)

                            ph_prom = 0
                            ph_t = 0
                            count = 0
                            r_name = []
                            r_seq = []
                            r_seqs_com = []
                            r_seqs_com2 = []
                            update_seq = ""

        else:

            with gzip.open(sequence_file1, "rt") as in_file:

                for line in in_file:

                    count += 1
                    line_mod = line.strip()
                    self.root.update()

                    if count == 1:

                        r_name.append(line_mod)
                        r_name_dict = line_mod

                    if count == 2:

                        r_seq.append(line_mod)
                        r_comp = complement_fnc(line_mod)
                        r_seqs_com.append(line_mod)
                        r_seqs_com.append(r_comp)
                        r_seqs_com.sort()

                        if r_seqs_com[0] not in seq_name_dict:

                            update_seq = r_seqs_com[0]
                            number = 1
                            new_element2 = {update_seq: number}
                            seq_mult_dict.update(new_element2)

                        if r_seqs_com[0] in seq_name_dict:

                            number = seq_mult_dict.get(r_seqs_com[0])
                            number += 1
                            update_seq = r_seqs_com[0]
                            new_element3 = {update_seq: number}
                            seq_mult_dict.update(new_element3)

                    new_element = {update_seq: r_name_dict}

                    if update_seq != "":

                        seq_name_dict.update(new_element)

                    if count == 4:

                        for x in line_mod:

                            phred = int(ord(x) - 33)
                            ph_t = ph_t + int(phred)

                        ph_prom = round(ph_t / int(len(line_mod)))
                        r_comp2 = complement_fnc(r_seq[0])
                        r_seqs_com2.append(r_seq[0])
                        r_seqs_com2.append(r_comp2)
                        r_seqs_com2.sort()

                        if r_seqs_com2[0] not in seq_ph_dict:

                            new_element4 = {r_seqs_com2[0]: ph_prom}
                            seq_ph_dict.update(new_element4)

                        n = seq_mult_dict.get(r_seqs_com2[0])

                        if r_seqs_com2[0] in seq_ph_dict and n > 1:

                            number_2 = seq_ph_dict.get(r_seqs_com[0])
                            number_2 += ph_prom
                            new_element6 = {r_seqs_com2[0]: number_2}
                            seq_ph_dict.update(new_element6)

                        ph_prom = 0
                        ph_t = 0
                        count = 0
                        r_name = []
                        r_seq = []
                        r_seqs_com = []
                        r_seqs_com2 = []
                        update_seq = ""

    else:

        if sequence_file2 != "":

            filenames = [sequence_file1, sequence_file2]

            for fname in filenames:

                with open(fname, "r") as infile:

                    for line in infile:

                        count += 1
                        line_mod = line.strip()
                        self.root.update()

                        if count == 1:

                            r_name.append(line_mod)
                            r_name_dict = line_mod

                        if count == 2:

                            r_seq.append(line_mod)
                            r_comp = complement_fnc(line_mod)
                            r_seqs_com.append(line_mod)
                            r_seqs_com.append(r_comp)
                            r_seqs_com.sort()

                            if r_seqs_com[0] not in seq_name_dict:

                                update_seq = r_seqs_com[0]
                                number = 1
                                new_element2 = {update_seq: number}
                                seq_mult_dict.update(new_element2)

                            if r_seqs_com[0] in seq_name_dict:

                                number = seq_mult_dict.get(r_seqs_com[0])
                                number += 1
                                update_seq = r_seqs_com[0]
                                new_element3 = {update_seq: number}
                                seq_mult_dict.update(new_element3)

                        new_element = {update_seq: r_name_dict}

                        if update_seq != "":

                            seq_name_dict.update(new_element)

                        if count == 4:

                            for x in line_mod:

                                phred = int(ord(x) - 33)
                                ph_t = ph_t + int(phred)

                            ph_prom = round(ph_t / int(len(line_mod)))
                            r_comp2 = complement_fnc(r_seq[0])
                            r_seqs_com2.append(r_seq[0])
                            r_seqs_com2.append(r_comp2)
                            r_seqs_com2.sort()

                            if r_seqs_com2[0] not in seq_ph_dict:

                                new_element4 = {r_seqs_com2[0]: ph_prom}
                                seq_ph_dict.update(new_element4)

                            n = seq_mult_dict.get(r_seqs_com2[0])

                            if r_seqs_com2[0] in seq_ph_dict and n > 1:

                                number_2 = seq_ph_dict.get(r_seqs_com[0])
                                number_2 += ph_prom
                                new_element6 = {r_seqs_com2[0]: number_2}
                                seq_ph_dict.update(new_element6)

                            ph_prom = 0
                            ph_t = 0
                            count = 0
                            r_name = []
                            r_seq = []
                            r_seqs_com = []
                            r_seqs_com2 = []
                            update_seq = ""

        else:

            with open(sequence_file1, "r") as in_file:

                for line in in_file:

                    count += 1
                    line_mod = line.strip()
                    self.root.update()

                    if count == 1:

                        r_name.append(line_mod)
                        r_name_dict = line_mod

                    if count == 2:

                        r_seq.append(line_mod)
                        r_comp = complement_fnc(line_mod)
                        r_seqs_com.append(line_mod)
                        r_seqs_com.append(r_comp)
                        r_seqs_com.sort()

                        if r_seqs_com[0] not in seq_name_dict:

                            update_seq = r_seqs_com[0]
                            number = 1
                            new_element2 = {update_seq: number}
                            seq_mult_dict.update(new_element2)

                        if r_seqs_com[0] in seq_name_dict:

                            number = seq_mult_dict.get(r_seqs_com[0])
                            number += 1
                            update_seq = r_seqs_com[0]
                            new_element3 = {update_seq: number}
                            seq_mult_dict.update(new_element3)

                    new_element = {update_seq: r_name_dict}

                    if update_seq != "":

                        seq_name_dict.update(new_element)

                    if count == 4:

                        for x in line_mod:

                            phred = int(ord(x) - 33)
                            ph_t = ph_t + int(phred)

                        ph_prom = round(ph_t / int(len(line_mod)))
                        r_comp2 = complement_fnc(r_seq[0])
                        r_seqs_com2.append(r_seq[0])
                        r_seqs_com2.append(r_comp2)
                        r_seqs_com2.sort()

                        if r_seqs_com2[0] not in seq_ph_dict:

                            new_element4 = {r_seqs_com2[0]: ph_prom}
                            seq_ph_dict.update(new_element4)

                        n = seq_mult_dict.get(r_seqs_com2[0])

                        if r_seqs_com2[0] in seq_ph_dict and n > 1:

                            number_2 = seq_ph_dict.get(r_seqs_com[0])
                            number_2 += ph_prom
                            new_element6 = {r_seqs_com2[0]: number_2}
                            seq_ph_dict.update(new_element6)

                        ph_prom = 0
                        ph_t = 0
                        count = 0
                        r_name = []
                        r_seq = []
                        r_seqs_com = []
                        r_seqs_com2 = []
                        update_seq = ""

    seqs = seq_ph_dict.keys()

    for i in seqs:

        self.root.update()
        phred = seq_ph_dict.get(i)
        count3 = seq_mult_dict.get(i)
        mean = float(round(phred / count3, 2))
        new_element7 = {i: mean}
        seq_ph_a_dict.update(new_element7)

    seqs__ = sorted(seq_mult_dict.items(), key=itemgetter(1), reverse=True)

    dictionary = {}
    convert_fnc(seqs__, dictionary)
    seqs5 = dictionary.keys()
    count5 = 0

    nr_name = name.split(".")[0]
    self.name_main = nr_name
    result_folder = "Results_" + self.name_main
    path_to_specific_results = os.path.join(
        dir_path, result_path, result_folder)
    self.path_to_sr = path_to_specific_results

    check_result_folder_spec = os.path.isdir(path_to_specific_results)

    if not check_result_folder_spec:

        os.mkdir(path_to_specific_results)

    path_to_nr = os.path.join(dir_path, result_path, self.path_to_sr, nr_name)

    self.statusbar.config(text="Saving non redundant sequences")
    print("Saving non redundant sequences...\n")

    for i in seqs5:

        self.root.update()
        count5 += 1

        ph_s = seq_ph_a_dict.get(i)
        ss = dictionary.get(i)
        ss_m = str(ss)
        ss_mod = ss_m.replace("[", "", 1)
        ss_mod2 = ss_mod.replace("]", "", 1)

        with open(path_to_nr + "_nr.fasta", "a+") as fasta_file:

            print(
                ">SEQ_" +
                str(count5) +
                '\t' +
                str(ss_mod2) +
                '\t' +
                str(ph_s),
                file=fasta_file)
            print(i + '\n', file=fasta_file)

    file1_r = path_to_nr + "_nr.fasta"

    t = threading.Thread(target=rna_seq_aligner(self, file1_r))
    t.start()


def verify(self):
    '''Is Magic-BLAST running?'''

    global check_mb_proc
    global ver
    time_mod = os.path.getmtime(path_to_results_magicblast)
    checker_time.append(time_mod)

    if len(checker_time) > 40:

        last = checker_time[-1]
        last_two_minutes = checker_time[-36]

        if last == last_two_minutes:

            check_mb_proc = 3

        if check_mb_proc == 3:

            ver = 3  # Magic-BLAST process is finished


def sam_processor(self):

    def string2list(string):
        '''Convert string into to list'''

        return [char for char in string]

    def list2string(list):
        '''Convert list into to string'''

        string = ""

        for elements in list:

            string += str(elements)

        return string

    def translate_cigar_code(untranslate, refSeq):
        '''Translate numeric scaffolds'''

        cont = 0
        scaffold_l = []

        for i in untranslate:

            if i == 1:

                scaffold_l.append(refSeq[cont].capitalize())

            if i == 2:

                scaffold_l.append("a")

            if i == 3:

                scaffold_l.append("t")

            if i == 5:

                scaffold_l.append("c")

            if i == 7:

                scaffold_l.append("g")

            if i == 9:

                scaffold_l.append("-")

            if i == 0:

                scaffold_l.append("n")

            cont += 1

        return (list2string(scaffold_l))

    global flag_sam_proc
    global limit
    global pros

    if flag_sam_proc == 0:

        flag_sam_proc = 1
        '''Relevant information is extracted from the database'''
        col_headers2 = [
            'code',
            'description',
            'family',
            'genus',
            'pubmed',
            'long',
            'host',
            'isolate',
            'country',
            'segment',
            'status',
            'gc',
            'seq']
        data_db = pd.read_csv(
            path_to_db_tsv,
            delimiter='\t',
            names=col_headers2)

        '''TSV database information is saved in dictionaries (detected sequences only)'''
        descritp = data_db["description"].tolist()
        countries = data_db["country"].tolist()
        families = data_db["family"].tolist()
        segment = data_db["segment"].tolist()
        genus_l = data_db["genus"].tolist()
        pubmed = data_db["pubmed"].tolist()
        length = data_db["long"].tolist()
        host = data_db["host"].tolist()
        code = data_db["code"].tolist()
        seq = data_db["seq"].tolist()
        dict_code_description = dict(zip(code, descritp))
        dict_code_countries = dict(zip(code, countries))
        dict_code_segment = dict(zip(code, segment))
        dict_code_family = dict(zip(code, families))
        dict_code_genus = dict(zip(code, genus_l))
        dict_code_pubmed = dict(zip(code, pubmed))
        dict_code_length = dict(zip(code, length))
        dict_code_host = dict(zip(code, host))
        dict_code_seq = dict(zip(code, seq))

        ''''Dictionaries to use are initialized'''
        dict_code_multiplicity = {}
        dict_code_2average = {}
        dict_code_phred = {}

        if os.path.isfile(path_to_results_short_sam):

            os.remove(path_to_results_short_sam)

        name_nr = self.name_main
        name_file_consensus = name_nr + "_mapping_assembly.fasta"
        name_file_summary = name_nr + "_summary.tsv"
        path_to_RGMA = os.path.join(
            dir_path,
            result_path,
            self.path_to_sr,
            name_file_consensus)
        path_to_summary = os.path.join(
            dir_path,
            result_path,
            self.path_to_sr,
            name_file_summary)

        '''Indexing'''
        with open(path_to_results_magicblast, "r") as lines:

            for line in lines:

                locus = '^@'
                locus_hit = re.match(locus, line)

                if locus_hit:

                    pass

                else:

                    x = line.split()

                    '''Filter by FLAG and PHRED. If not primary out. If not PQ < 30 out'''
                    if x[3] != '256' and float(x[2]) > 30:

                        '''SHORT-SAM file is created'''
                        with open(path_to_results_short_sam, "a+") as fasta_file:

                            print(x[0]
                                  + '\t'
                                  + x[1]
                                  + '\t'
                                  + x[2]
                                  + '\t'
                                  + x[4]
                                  + '\t'
                                  + x[5]
                                  + '\t'
                                  + x[11], file=fasta_file)
                    else:

                        pass

        col_headers = [
            'ID',
            'MULTIPLICITY',
            'PHRED',
            'REFERENCE',
            'POSITION',
            'READ']

        if os.path.isfile(path_to_results_short_sam):

            self.statusbar.config(text="Processing alignment")
            print("Processing alignment...\n")

            '''SHORT-SAM is sorted by reference and position'''
            sam = pd.read_csv(
                path_to_results_short_sam,
                delimiter='\t',
                header=None,
                names=col_headers)
            sam_ordered = sam.sort_values(by=['REFERENCE', 'POSITION'])
            codes = sam_ordered.REFERENCE.unique()

            if pros == 0:

                for i in codes:
                    '''Memory array is created'''

                    ncol = dict_code_length.get(i)
                    temp_df = pd.DataFrame(
                        0, index=range(12), columns=range(
                            ncol + 10))
                    seqid = sam_ordered.loc[sam_ordered['REFERENCE'] == i]
                    check_size = dict_code_length.get(i)
                    ref_seq = dict_code_seq.get(i)
                    how_much_rows = len(seqid.index)

                    if how_much_rows == 1:

                        pass

                    else:

                        for index, row in seqid.iterrows():
                            '''Memory matrix is completed for the respective unique code'''

                            self.root.update()
                            multiplicity = row["MULTIPLICITY"]
                            phred_average = row["PHRED"]
                            start = int(row["POSITION"])
                            read_i = row["READ"]

                            if (start + len(read_i)) > check_size:

                                pass

                            else:

                                subject_var = dict_code_seq.get(i)

                                '''Files are created to run BLASTn'''
                                with open(path_to_temp_subject, "w") as subject_file:

                                    print(">" + i, file=subject_file)
                                    print(subject_var, file=subject_file)

                                with open(path_to_temp_query, "w") as query_file:

                                    print(">read_nr", file=query_file)
                                    print(read_i, file=query_file)

                                '''Quick BLASTn'''
                                command2 = ("blastn -subject "
                                            + path_to_temp_subject
                                            + " -max_hsps 1 "
                                            + " -outfmt 3 -query "
                                            + path_to_temp_query
                                            + "  -out  "
                                            + path_to_temp_query_result_blastn)
                                blast_n = subprocess.Popen(
                                    command2, shell=False, stderr=subprocess.DEVNULL)
                                blast_n.wait()

                                with open(path_to_temp_query_result_blastn, 'r') as result_blastn:

                                    start_subject = 2000000
                                    full_subject = ""
                                    full_query = ""
                                    end_subject = 0

                                    for line in result_blastn:
                                        '''BLASTn result file is processed with regular expressions'''

                                        query_ = '^Query_1'
                                        query_hit = re.match(query_, line)

                                        if query_hit:

                                            full_query += line.strip().split()[
                                                2]

                                        subject_ = '^Subject_1'
                                        subject_hit = re.match(subject_, line)

                                        if subject_hit:

                                            full_subject += line.strip().split()[
                                                2]

                                            if int(line.strip().split()
                                                   [1]) < start_subject:

                                                start_subject = int(
                                                    line.strip().split()[1])

                                            if int(line.strip().split()
                                                   [3]) > end_subject:

                                                end_subject = int(
                                                    line.strip().split()[3])

                                if i not in dict_code_multiplicity:

                                    dict_code_multiplicity[i] = multiplicity
                                    dict_code_phred[i] = phred_average
                                    dict_code_2average[i] = 1

                                else:

                                    counter = int(
                                        dict_code_multiplicity.get(i))
                                    counter = counter + int(multiplicity)
                                    element_mult = {i: counter}
                                    dict_code_multiplicity.update(element_mult)

                                    counter2 = float(dict_code_phred.get(i))
                                    counter2 = counter2 + float(phred_average)
                                    element_phred = {i: counter2}
                                    dict_code_phred.update(element_phred)

                                    counter3 = float(dict_code_2average.get(i))
                                    counter3 = counter3 + 1
                                    element_2ave = {i: counter3}
                                    dict_code_2average.update(element_2ave)

                                list_to_actualize_df = string2list(
                                    full_subject)
                                list_to_detect_changes = string2list(
                                    full_query)
                                pos_ini = start_subject - 1
                                counter4 = 0

                                for ic in list_to_actualize_df:
                                    '''The memory array is updated'''

                                    if ic == ".":

                                        temp_df.iat[0, pos_ini] = 1
                                        mod_multiplicity_1 = temp_df.loc[1].iat[pos_ini]
                                        mod_multiplicity_1 = mod_multiplicity_1 + multiplicity
                                        temp_df.iat[1,
                                                    pos_ini] = mod_multiplicity_1

                                    else:

                                        if list_to_detect_changes[counter4] == "A":

                                            temp_df.iat[2, pos_ini] = 2
                                            mod_multiplicity_2 = temp_df.loc[3].iat[pos_ini]
                                            mod_multiplicity_2 = mod_multiplicity_2 + multiplicity
                                            temp_df.iat[3,
                                                        pos_ini] = mod_multiplicity_2

                                        if list_to_detect_changes[counter4] == "T":

                                            temp_df.iat[4, pos_ini] = 3
                                            mod_multiplicity_3 = temp_df.loc[5].iat[pos_ini]
                                            mod_multiplicity_3 = mod_multiplicity_3 + multiplicity
                                            temp_df.iat[5,
                                                        pos_ini] = mod_multiplicity_3

                                        if list_to_detect_changes[counter4] == "C":

                                            temp_df.iat[6, pos_ini] = 5
                                            mod_multiplicity_4 = temp_df.loc[7].iat[pos_ini]
                                            mod_multiplicity_4 = mod_multiplicity_4 + multiplicity
                                            temp_df.iat[7,
                                                        pos_ini] = mod_multiplicity_4

                                        if list_to_detect_changes[counter4] == "G":

                                            temp_df.iat[8, pos_ini] = 7
                                            mod_multiplicity_5 = temp_df.loc[9].iat[pos_ini]
                                            mod_multiplicity_5 = mod_multiplicity_5 + multiplicity
                                            temp_df.iat[9,
                                                        pos_ini] = mod_multiplicity_5

                                        if list_to_detect_changes[counter4] == "-":

                                            temp_df.iat[10, pos_ini] = 9
                                            mod_multiplicity_5 = temp_df.loc[11].iat[pos_ini]
                                            mod_multiplicity_5 = mod_multiplicity_5 + multiplicity
                                            temp_df.iat[11,
                                                        pos_ini] = mod_multiplicity_5

                                    pos_ini += 1
                                    counter4 += 1

                                if os.path.isfile(
                                        path_to_temp_query_result_blastn):

                                    os.remove(path_to_temp_query_result_blastn)

                                if os.path.isfile(path_to_temp_subject):

                                    os.remove(path_to_temp_subject)

                                if os.path.isfile(path_to_temp_query):

                                    os.remove(path_to_temp_query)

                        if i in dict_code_multiplicity:

                            transpose_temp_df = temp_df.transpose()
                            scaffold_over_march = []
                            list_of_multi = []

                            for index2, row2 in transpose_temp_df.iterrows():
                                '''Select the nucleotide most frequently in a particular position'''

                                multiplicity_of_1 = int(row2[1])
                                multiplicity_of_A = int(row2[3])
                                multiplicity_of_T = int(row2[5])
                                multiplicity_of_C = int(row2[7])
                                multiplicity_of_G = int(row2[9])
                                multiplicity_of_d = int(row2[11])

                                if multiplicity_of_1 > multiplicity_of_A and multiplicity_of_1 > multiplicity_of_T and multiplicity_of_1 > multiplicity_of_C and multiplicity_of_1 > multiplicity_of_G\
                                        and multiplicity_of_1 > multiplicity_of_d:

                                    scaffold_over_march.append(int(row2[0]))
                                    list_of_multi.append(int(row2[1]))

                                if multiplicity_of_A > multiplicity_of_1 and multiplicity_of_A > multiplicity_of_T and multiplicity_of_A > multiplicity_of_C and multiplicity_of_A > multiplicity_of_G\
                                        and multiplicity_of_A > multiplicity_of_d:

                                    scaffold_over_march.append(int(row2[2]))
                                    list_of_multi.append(int(row2[3]))

                                if multiplicity_of_T > multiplicity_of_A and multiplicity_of_T > multiplicity_of_1 and multiplicity_of_T > multiplicity_of_C and multiplicity_of_T > multiplicity_of_G\
                                        and multiplicity_of_T > multiplicity_of_d:

                                    scaffold_over_march.append(int(row2[4]))
                                    list_of_multi.append(int(row2[5]))

                                if multiplicity_of_C > multiplicity_of_A and multiplicity_of_C > multiplicity_of_T and multiplicity_of_C > multiplicity_of_1 and multiplicity_of_C > multiplicity_of_G\
                                        and multiplicity_of_C > multiplicity_of_d:

                                    scaffold_over_march.append(int(row2[6]))
                                    list_of_multi.append(int(row2[7]))

                                if multiplicity_of_G > multiplicity_of_A and multiplicity_of_G > multiplicity_of_T and multiplicity_of_G > multiplicity_of_C and multiplicity_of_G > multiplicity_of_1\
                                        and multiplicity_of_G > multiplicity_of_d:

                                    scaffold_over_march.append(int(row2[8]))
                                    list_of_multi.append(int(row2[9]))

                                if multiplicity_of_d > multiplicity_of_1 and multiplicity_of_d > multiplicity_of_A and multiplicity_of_d > multiplicity_of_T and multiplicity_of_d > multiplicity_of_C\
                                        and multiplicity_of_d > multiplicity_of_G:

                                    scaffold_over_march.append(int(row2[10]))
                                    list_of_multi.append(int(row2[11]))

                                if multiplicity_of_1 == multiplicity_of_A and multiplicity_of_1 == multiplicity_of_T and multiplicity_of_1 == multiplicity_of_C and multiplicity_of_1 == multiplicity_of_G\
                                        and multiplicity_of_1 == multiplicity_of_d:

                                    scaffold_over_march.append(0)
                                    list_of_multi.append(int(row2[1]))

                            scaffolddd = translate_cigar_code(
                                scaffold_over_march, ref_seq)
                            mx_list = max(list_of_multi)
                            range_gen_max = int(len(list_of_multi))
                            positions_of_genome = []

                            for ixz in range(1, range_gen_max + 1, 1):

                                positions_of_genome.append(ixz)

                            d_graph = {
                                'Genome position (nt)': positions_of_genome,
                                'No. of reads': list_of_multi}
                            df_graph = pd.DataFrame(d_graph)
                            zeros = list_of_multi.count(0)
                            perct = (check_size - zeros) / check_size

                            if perct > limit:
                                '''If the scaffold coverage percentage is greater than 10 percent, it is added to the results'''

                                self.root.update()
                                title_graph = dict_code_description.get(
                                    i) + " [" + i + "]"
                                sns.relplot(
                                    x="Genome position (nt)",
                                    y="No. of reads",
                                    data=df_graph,
                                    kind='line',
                                    height=2,
                                    aspect=4.5,
                                    linewidth=1)
                                plt.title(title_graph)
                                plt.tight_layout()
                                plt.ylim(0, mx_list)
                                plt.fill_between(
                                    df_graph["Genome position (nt)"], df_graph["No. of reads"])
                                name_of_file = self.name_main + "_" + i + '.png'
                                path_to_image = os.path.join(
                                    self.path_to_sr, name_of_file)
                                plt.savefig(path_to_image)
                                list_of_images.append(path_to_image)
                                plt.clf()
                                plt.cla()
                                plt.close()
                                self.root.update()

                                with open(path_to_RGMA, "a+") as fasta_file:
                                    '''Save scaffold in file'''

                                    limit_perc = round(limit * 100)
                                    self.root.update()
                                    print(
                                        dict_code_description.get(i) +
                                        " [" +
                                        i +
                                        "]" +
                                        " scaffold, has a coverage greater than " +
                                        str(limit_perc) +
                                        " percent, based on the reference\n")
                                    phred_sum = float(dict_code_phred.get(i))
                                    phred_quantity = float(
                                        dict_code_2average.get(i))
                                    phred_average_total = round(
                                        (phred_sum / phred_quantity), 2)

                                    '''Undetermined nucleotides from the start and end of the genome are removed'''
                                    seq_ns = scaffolddd
                                    n_sf = '^[n]+'
                                    nsf_hit = re.match(n_sf, seq_ns)

                                    if nsf_hit:

                                        clearf = nsf_hit.group(0)
                                        cleanf = seq_ns.replace(clearf, "", 1)

                                    else:

                                        cleanf = seq_ns

                                    cadenaInvertida = cleanf[::-1]
                                    n_sr = '^[n]+'
                                    nsr_hit = re.match(n_sr, cadenaInvertida)

                                    if nsr_hit:

                                        clearr = nsr_hit.group(0)
                                        cleanr = cadenaInvertida.replace(
                                            clearr, "", 1)

                                    else:

                                        cleanr = cadenaInvertida

                                    seq_clear = cleanr[::-1]

                                    '''Scaffold coverage percentage is calculated'''
                                    count_n = seq_clear.count("n")
                                    all_chars = len(seq_clear)
                                    coverage = int(
                                        (((all_chars - count_n) / all_chars) * 100))
                                    exact_matches = len(
                                        re.findall(r'[A-Z]', seq_clear))
                                    all_non_n = all_chars - count_n
                                    identity_of_contigs = int(
                                        (((1 - (all_non_n - exact_matches) / all_non_n)) * 100))

                                    print(">" +
                                          str(i) +
                                          "\t" +
                                          str(dict_code_description.get(i)) +
                                          "\t" +
                                          "Segment: " +
                                          str(dict_code_segment.get(i)) +
                                          "\t" +
                                          "Length: " +
                                          str(dict_code_length.get(i)) +
                                          "\t" +
                                          "#Reads: " +
                                          str(dict_code_multiplicity.get(i)) +
                                          "\t" +
                                          "Phred: " +
                                          str(phred_average_total) +
                                          "\t" +
                                          "Scaffold length(Coverage): " +
                                          str(all_chars) +
                                          "(" +
                                          str(coverage) +
                                          "%)" +
                                          "\t" +
                                          "Avg. Identity: " +
                                          str(identity_of_contigs) +
                                          "%" +
                                          "\t" +
                                          "Family: " +
                                          str(dict_code_family.get(i)) +
                                          "\t" +
                                          "Genus: " +
                                          str(dict_code_genus.get(i)) +
                                          "\t" +
                                          "Host: " +
                                          str(dict_code_host.get(i)) +
                                          "\t" +
                                          "Country: " +
                                          str(dict_code_countries.get(i)) +
                                          "\t" +
                                          "Pubmed: " +
                                          str(dict_code_pubmed.get(i)), file=fasta_file)

                                    print(
                                        str(seq_clear) + '\n', file=fasta_file)

                                table_tsv = str(title_graph) + '\t' + str(dict_code_length.get(i)) + '\t' + str(all_chars) + " (" + str(
                                    coverage) + "%)" + '\t' + str(identity_of_contigs) + "%" + '\t' + str(dict_code_multiplicity.get(i))

                                with open(path_to_summary, "a+") as summary_file:

                                    print(table_tsv, file=summary_file)

                time.sleep(5)

                if os.path.isfile(path_to_summary):
                    '''Summary table image is created'''

                    data2 = pd.read_csv(
                        path_to_summary,
                        delimiter='\t',
                        names=[
                            'Name [ID]',
                            'ID length',
                            'Scaffold length(coverage)',
                            'Avg. Identity',
                            'Associate reads'])
                    data2.sort_values(
                        by=['Associate reads'],
                        inplace=True,
                        ascending=False)
                    fig, ax = plt.subplots(figsize=(9, 2))
                    ax.xaxis.set_visible(False)
                    ax.yaxis.set_visible(False)
                    ax.set_frame_on(False)
                    data2 = data2.reset_index(drop=True)
                    tabla = table(
                        ax,
                        data2,
                        loc='center',
                        cellLoc='center',
                        rowLoc='center',
                        rowLabels=[''] *
                        data2.shape[0])  # where df is your data frame
                    tabla.auto_set_font_size(False)
                    tabla.auto_set_column_width(
                        col=list(range(len(data2.columns))))
                    tabla.set_fontsize(8)
                    tabla.scale(1, 1)
                    name_of_summary = self.name_main + "_summary.png"
                    path_to_img_summary = os.path.join(
                        self.path_to_sr, name_of_summary)
                    plt.savefig(path_to_img_summary)
                    list_of_images.append(path_to_img_summary)
                    plt.clf()
                    plt.cla()
                    plt.close()
                    pros = 3
                    flag_sam_proc = 1
                    end = time.time()

                else:

                    self.statusbar.config(
                        text="No significant similarity found")
                    print("No significant similarity found!")
                    self.progress.stop()
                    tk.messagebox.showwarning(
                        "Process finished", "No significant similarity found.")

        else:

            self.statusbar.config(text="No significant similarity found")
            print("No significant similarity found!")
            self.progress.stop()
            tk.messagebox.showwarning(
                "Process finished",
                "No significant similarity found.")


def active_results_button(self):
    '''Process finished?'''

    global ver
    global pros
    global flag_sam_proc

    if ver == 3:  # If Magic-BLAST process is finished

        if flag_sam_proc == 0:

            t3 = threading.Thread(target=sam_processor(self))
            t3.start()

        if pros == 3:  # sam_processor is done

            ver = 4  # Pass this function
            pros = 1  # No process again

            if os.path.isfile(path_to_results_short_sam):

                os.remove(path_to_results_short_sam)

            if os.path.isfile(path_to_results_magicblast):

                os.remove(path_to_results_magicblast)

            end = time.time()
            minutes = round((end - start_time) / 60, 2)
            self.progress.stop()
            infooo = "Time of execution: " + str(minutes) + " minutes"
            print(infooo)
            self.statusbar.config(text=infooo)

            def ress(self):
                '''Results button window'''

                self.top_level = Toplevel(self.root)
                self.top_level.wait_visibility()
                self.top_level.geometry("900x540")
                self.top_level.title("VirusMAPPER - Results")
                self.top_level['bg'] = '#ccc2bd'
                self.top_level.resizable(False, False)
                top_level = self.top_level
                list_of_grapgs = []

                for iasx in list_of_images:

                    temp_img = ImageTk.PhotoImage(Image.open(iasx))
                    list_of_grapgs.append(temp_img)

                stop_next = len(list_of_grapgs)

                self.img_lbl = ttk.Label(top_level, image=list_of_grapgs[0])
                self.img_lbl.place(x=0, y=0, width=900, height=540)

                self.btn7_r = ttk.Button(top_level,
                                         text="<<",
                                         state="disabled")
                self.btn7_r.place(x=360, y=510)

                self.btn8_r = ttk.Button(top_level,
                                         text=">>",
                                         command=lambda: next_(self, 2))
                self.btn8_r.place(x=460, y=510)

                def next_(self, image_number):
                    '''Results window, next button'''

                    self.img_lbl.place_forget()
                    self.img_lbl.config(image=list_of_grapgs[image_number - 1])
                    self.img_lbl.place(x=0, y=0, width=900, height=560)

                    self.btn7_r.config(
                        command=lambda: _back(
                            self, image_number - 1))
                    self.btn8_r.config(
                        command=lambda: next_(
                            self, image_number + 1))

                    if image_number == stop_next:

                        self.btn8_r.config(state="disabled")

                    if image_number != 1:

                        self.btn7_r.config(state="normal")

                    if image_number == 1:

                        self.btn8_r.config(state="normal")

                def _back(self, image_number):
                    '''Results window, back button'''

                    self.img_lbl.place_forget()
                    self.img_lbl.config(image=list_of_grapgs[image_number - 1])
                    self.img_lbl.place(x=0, y=0, width=900, height=560)

                    self.btn7_r.config(
                        command=lambda: _back(
                            self, image_number - 1))
                    self.btn8_r.config(
                        command=lambda: next_(
                            self, image_number + 1))

                    if image_number == 1:

                        self.btn7_r.config(state="disabled")

                    if image_number != stop_next:

                        self.btn8_r.config(state="normal")

                    if image_number == stop_next:

                        self.btn7_r.config(state="normal")

            self.btn_4.config(state=NORMAL, command=lambda: ress(self))
            self.root.update()
            tk.messagebox.showinfo(
                message="The results were saved in the working directory",
                title="Process finished")


def open_query(self):
    '''Open file 1 button'''

    filename1 = filedialog.askopenfilename(
        title="Select a file", filetypes=[
            ("fastq files", "*.fastq *.fastq.gz *.fq")])
    self.route_file1 = filename1
    routeInv = filename1[::-1]
    routeFrag = routeInv.split("/")
    name_neg = routeFrag[0]
    name_fix = name_neg[::-1]
    self.lbl1.config(text=name_fix)
    self.name_file1 = name_fix
    self.btn_3.config(state=NORMAL)


def open_querymate(self):
    '''Open file 2 button'''

    filename2 = filedialog.askopenfilename(
        title="Select a file", filetypes=[
            ("fastq files", "*.fastq *.fastq.gz *.fq")])
    self.route_file2 = filename2
    routeInv2 = filename2[::-1]
    routeFrag2 = routeInv2.split("/")
    name_neg2 = routeFrag2[0]
    name_fix2 = name_neg2[::-1]
    self.lbl2.config(text=name_fix2)
    self.name_file2 = name_fix2


def process(self):
    '''Process button'''

    check_file = self.lbl1.cget("text")

    if len(check_file) == 0:

        tk.messagebox.showerror("Error", "No files")

    elif self.route_file2 == "" and self.value == 2:

        tk.messagebox.showerror("Error", "Select both files")

    else:

        self.progress = ttk.Progressbar(self.root, mode="indeterminate")
        self.progress.place(x=20, y=240, width=660)
        self.progress.start()
        self.btn_1.config(state=DISABLED)
        self.btn_2.config(state=DISABLED)
        self.btn_3.config(state=DISABLED)
        self.radio_btn1.config(state=DISABLED)
        self.radio_btn2.config(state=DISABLED)
        t1 = threading.Thread(target=non_red_seqs(self))
        t1.start()
        self.root.update_idletasks()


def click_quantity(self, value):
    '''Single or paired-end reads'''

    self.value = value

    if value == 1:

        self.btn_1.config(state="normal")
        self.btn_2.config(state="disabled")

    if value == 2:

        self.btn_1.config(state="normal")
        self.btn_2.config(state="normal")


class VirusMapperUI:
    '''GUI of VirusMAPPER'''

    def __init__(self, root):

        s = ttk.Style()
        s.configure('.', font=("Book Antiqua", 11))
        self.root = root
        self.root.geometry("700x315")
        self.root.title("VirusMAPPER")
        self.root['bg'] = '#ccc2bd'
        self.root.resizable(False, False)
        ico_name = "ico_virusmapper_.ico"
        path_to_ico = os.path.join(dir_path, db_path, ico_name)
        self.root.iconbitmap(path_to_ico)
        self.value = 0
        self.name_file1 = ""
        self.name_file2 = ""
        self.route_file1 = ""
        self.route_file2 = ""
        self.message = "message"
        self.name_main = ""
        self.path_to_sr = ""
        self.statusbar = ttk.Label(
            root,
            text="Welcome",
            relief=SUNKEN,
            borderwidth=1,
            anchor=CENTER)
        self.statusbar.place(x=20, y=270, width=660)
        self.tabs_of_gui()
        self.widgets_tab_1()
        self.widgets_tab_2()
        self.widgets_tab_3()
        self.periodic_call()

    def periodic_call(self):
        '''Execute this every 5 seconds'''

        self.root.after(5000, self.periodic_call)
        self.root.update_idletasks()

        if ver == 1:

            verify(self)

        if ver == 3 and flag_sam_proc == 0:

            active_results_button(self)

    def tabs_of_gui(self):

        self.tab_control = ttk.Notebook(self.root)
        self.tab = ttk.Frame(self.tab_control)
        self.tab_1 = ttk.Frame(self.tab_control)
        self.tab_2 = ttk.Frame(self.tab_control)
        self.tab_control.add(self.tab, text="Home")
        self.tab_control.add(self.tab_1, text="RNA-seq aligner")
        self.tab_control.add(self.tab_2, text="About us")
        self.tab_control.place(x=20,
                               y=20,
                               width=660,
                               height=210,
                               bordermode=INSIDE)

    def widgets_tab_1(self):

        global r
        r = IntVar()
        r.set("0")

        self.information = ttk.Label(
            self.tab_1,
            text='''This program support single or paired-end reads from Illumina\
 sequencing technology \nin FASTQ format (or fastq.gz):''')
        self.information.place(x=10, y=10)

        self.radio_btn1 = ttk.Radiobutton(
            self.tab_1,
            text="Single reads",
            variable=r,
            value=1,
            command=lambda: click_quantity(
                self,
                r.get()))
        self.radio_btn1.place(x=130, y=60)

        self.btn_1 = ttk.Button(
            self.tab_1,
            text="Open file 1",
            state="disabled",
            command=lambda: open_query(self))
        self.btn_1.place(x=310, y=60)

        self.lbl1 = ttk.Label(self.tab_1,
                              text="No file selected")
        self.lbl1.place(x=425, y=60)

        self.radio_btn2 = ttk.Radiobutton(
            self.tab_1,
            text="Paired-end reads",
            variable=r,
            value=2,
            command=lambda: click_quantity(
                self,
                r.get()))
        self.radio_btn2.place(x=130, y=90)

        self.btn_2 = ttk.Button(
            self.tab_1,
            text="Open file 2",
            state="disabled",
            command=lambda: open_querymate(self))
        self.btn_2.place(x=310, y=90)

        self.lbl2 = ttk.Label(self.tab_1,
                              text="No file selected")
        self.lbl2.place(x=425, y=90)

        self.btn_3 = ttk.Button(
            self.tab_1,
            text="Process",
            state="disabled",
            command=lambda: process(self))
        self.btn_3.place(x=210, y=140)

        self.btn_4 = ttk.Button(self.tab_1, text="Results", state="disabled")
        self.btn_4.place(x=350, y=140)

    def widgets_tab_2(self):

        self.information2 = ttk.Label(
            self.tab,
            text='''\nVirusMAPPER is able to align reads to a database with genomes of reference.\n\n The result obtained includes a fasta file \
with the generated scaffolds, coverage images for the \n mapped references, non-redundant sequence file and a summary of what was found.\
    \n\n If you have any trouble please feel free to contact us.''')
        self.information2.place(x=10, y=10)

    def widgets_tab_3(self):

        self.information2 = ttk.Label(
            self.tab_2, text='''Written by Ary Rivillas [amrivillast@unal.edu.co] [github.com/ARYrt]
Biotecnologia Microbiana Research Group (COL0001119)
Universidad Nacional de Colombia - Sede Medellín

                                                                                               Financed by
Secretaria de Agricultura y Desarrollo Rural de Antioquia. MinCiencias. CES and Universidad Nacional de Colombia
Sede Medellín. Convenio No. 4600007658-779. Project: "Desarrollo de una plataforma molecular y bioinformática para
el diagnóstico de virus en cultivos y material de siembra de papa (Solanum tuberosum y S. phureja) en Antioquia"
(Código 1101-805-62787).''', font=("Book Antiqua", 9))
        self.information2.place(x=10, y=10)


if __name__ == "__main__":

    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-
#
# @created: 08.09.2021
# @author: Danil
# @contact: Zilov

import argparse
import os


def mkblastdb(assembly_fasta, debug=False):
    blast_db_file = assembly_fasta + ".ndb"
    command = f"makeblastdb -in {assembly_fasta} -dbtype nucl -parse_seqids"
    if os.path.exists(blast_db_file):
        print("Database already created, continue!\n")
        return True
    else:
        if debug:
            print(f"\nCommand to make blast_db\n\n{command}\n")
            return False
        else:
            print(f"\nRunning {command}\n")
            os.system(command)
            print("mtDNA blast_db id created!\n")
            return True


def run_blast_nt(reference_mtdna, blast_db, threads, output_file, debug=False):
    command = f"blastn -query {reference_mtdna} -db {blast_db} -max_target_seqs 50 " \
              f"-outfmt '6 qseqid sseqid pident qstart qend length evalue sscinames staxids' " \
              f"-evalue 1e-5 -num_threads {threads} -out {output_file}"
    if os.path.exists(output_file):
        print("Already blasted, continue!\n")
        return True
    else:
        if debug:
            print(f"\nCommand to find mtDNA in assembly\n\n{command}\n")
            return False
        else:
            print(f"\nRunning {command}\n")
            os.system(command)
            print("Blast of assembly content is done!\n")
            return True


def fasta_contig_finder(genome, header_to_find):
    header_to_find = f">{header_to_find}"
    with open(genome) as fh:
        for i, line in enumerate(fh):
            line = line.strip()
            if line.startswith(">"):
                header = line.split()[0]
                if not header == header_to_find:
                    header = None
                else:
                    seq = []
            else:
                if header:
                    seq.append(line)
    return "".join(seq)


def sort_blast_results(blast_results):
    results = []
    with open(blast_results) as fh:
        for line in fh:
            line = line.strip().split("\t")
            results.append([line[1], line[2], line[3], line[4], line[5], line[6]])
        results.sort(key=lambda x: int(x[-2]), reverse=True)
        return results


def write_top_matches(sorted_results_list, results_summary):
    with open(results_summary, "w") as fw:
        header = "\t".join(["#contig", "identity", "start", "stop", "length", "p_value"])
        fw.write(f'{header}\n')
        for match in sorted_results_list[:15]:
            string_to_write = "\t".join(match)
            fw.write(f'{string_to_write}\n')
        print(f"\nTop matches was written in {results_summary}!")


def write_mt_fasta(assembly_fasta, sorted_blast_results, mtdna_file):
    top_match_header = sorted_blast_results[0][0]
    sequence = fasta_contig_finder(assembly_fasta, top_match_header)
    with open(mtdna_file, "w") as fw:
        fw.write(f">{top_match_header}\n{sequence}\n")
    print(f"\nmtDNA was written in {mtdna_file}!")


def run_mtdna_search(assembly_fasta, mtdna_fasta, blast_results, results_summary, mtdna_results, threads, debug):
    database = mkblastdb(assembly_fasta, debug)
    if database:
        blast = run_blast_nt(mtdna_fasta, assembly_fasta, threads, blast_results, debug)
        if blast:
            sorted_results = sort_blast_results(blast_results)
            write_top_matches(sorted_results, results_summary)
            write_mt_fasta(assembly_fasta, sorted_results, mtdna_results)
            print("\nDONE!\n")
            return True
        else:
            print("Error in blast running!")
    else:
        print("Error in makeblastdb running!")

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='mt_finder -  a tool for mtDNA extracting from genome')
    parser.add_argument('-a', '--assembly', help="path to genome asssembly in FASTA format", required=True)
    parser.add_argument('-m', '--mtdna', help="path to reference mtDNA in FASTA format", required=True)
    parser.add_argument('-p', '--prefix', help="prefix for output files", default=None)
    parser.add_argument('-o', '--outdir', help='output directory', required=True)
    parser.add_argument('-t', '--threads', help='number of threads [default == 8]', default="8")
    parser.add_argument('-d', '--debug', help='do not run, print commands only', action='store_true', default=False)

    args = vars(parser.parse_args())

    assembly_fasta = os.path.abspath(args["assembly"])
    mtdna_fasta = os.path.abspath(args["mtdna"])
    threads = args["threads"]
    debug = args["debug"]
    prefix = args["prefix"]
    if not prefix:
        prefix = os.path.splitext(os.path.basename(assembly_fasta))[0]
    outdir = os.path.abspath(args["outdir"])
    if not os.path.exists(outdir):
        os.makedirs(outdir)


    # define output files
    blast_results = os.path.join(outdir, f"{prefix}_mtdna_blast.outfmt6")
    top_results_summary = os.path.join(outdir, f"{prefix}_mtdna_top_blast.tsv")
    mtdna_results = os.path.join(outdir, f"{prefix}_mtdna_contig.fasta")

    # run pipeline
    run_mtdna_search(assembly_fasta, mtdna_fasta, blast_results, top_results_summary, mtdna_results, threads, debug)

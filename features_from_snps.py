# -*- coding: utf-8 -*-
"""
By J. Lucas Boatwright

Find any genes near variant positions
"""

import os
import sys
import sqlite3
import argparse
import gffutils
from tqdm import tqdm
from collections import OrderedDict

SQLITE_PATH = "Sbicolor_454_v3.1.1.gene.gff3.sqlite3"
ANNOTATION_PATH = "Sbicolor_454_v3.1.1.annotation_info.txt"


def parse_arguments():
    """Parse arguments passed to script"""
    parser = argparse.ArgumentParser(
        description=f"This script is designed to find any genes or features "
                    f"near variants or within genomic ranges.\n\n")

    required_named = parser.add_argument_group('required arguments')

    required_named.add_argument(
        "--input",
        type=str,
        required=True,
        help="The name of the input file (CSV).",
        action="store")

    required_named.add_argument(
        "--input-type",
        type=str,
        required=True,
        help=f"The type of input file: gapit, coordinate, range. " 
             f"Note: headers are **required** - or first line lost. " 
             f"*gapit - the unedited GAPIT output with header* " 
             f"*coordinate - markerName,chromosome,position *** " 
             f"*range - markerName,chromosome,position,start,end* ",
        action="store")

    required_named.add_argument(
        "--output",
        type=str,
        required=True,
        help="The output file to be created (TSV).",
        action="store")

    parser.add_argument(
        "--gff",
        type=str,
        required=False,
        default=SQLITE_PATH,
        help="The GFF3 file for the genome or a gffutils database.",
        action="store")

    parser.add_argument(
        "--distance",
        type=float,
        required=False,
        default=10.0,
        help="The distance in kb from a SNP to search for genes (default: 10).",
        action="store")

    parser.add_argument(
        "--fdr",
        type=float,
        required=False,
        default=0.05,
        help=f"The significance threshold for associated SNPs. Only " 
             f"applicable with input_type gapit.",
        action="store")

    parser.add_argument(
        "--info",
        type=str,
        required=False,
        default=ANNOTATION_PATH,
        help="The gene info file as obtained from Phytozome (TSV).",
        action="store")

    parser.add_argument(
        "--feature-type",
        type=str,
        required=False,
        help="The type of features desired from the annotation (i.e. gene,CDS).",
        default="gene",
        action="store")

    return parser.parse_args()


def file_to_list(filename):
    with open(filename) as f:
        output = f.read().splitlines()
    return output


def get_gene_info_dictionary(info):
    if not os.path.exists(ANNOTATION_PATH):
        sys.stderr.write(f"Annotation file {ANNOTATION_PATH} does not exist and is required.\n")
        sys.exit(1)
    gene_info = {}
    with open(info) as f:
        for line in f.read().splitlines():
            split_line = line.split()
            try:
                gene_info[split_line[1]] = gene_info[split_line[1]] + [line]
            except KeyError:
                gene_info[split_line[1]] = [line]
    return gene_info


def open_gff_db(gff):
    try:
        db = gffutils.FeatureDB(gff)
    except sqlite3.DatabaseError:
        if os.path.exists(f"{gff}.sqlite3"):
            sys.stderr.write(f"File {gff}.sqlite3 exists. Using existing database.\n")
            db = gffutils.FeatureDB(f"{gff}.sqlite3")
        else:
            db = gffutils.create_db(gff, f"{gff}.sqlite3")
    return db


def get_significant_snps(gwas, max_fdr):
    significant_snps = []
    for line in gwas:
        if not line.startswith("SNP"):
            split_line = line.split(',')
            current_fdr = float(split_line[8])
            if current_fdr <= max_fdr:
                significant_snps.append(line)
    return significant_snps


def generate_snp_range(sig_gwas, input_type, distance):
    snp_groups = OrderedDict()

    for x, line in enumerate(sig_gwas):
        split_line = line.split(',')
        snp_group = split_line[0]
        chromosome = split_line[1]
        position = int(split_line[2])
        if input_type == "range":
            snp_groups[snp_group] = [
                chromosome,
                position,
                int(split_line[3]),
                int(split_line[4])]
        elif input_type in ["gapit", "coordinate"]:
            snp_groups[snp_group] = [
                chromosome,
                position,
                position - distance,
                position + distance]
    return snp_groups


def coordinates_to_chromosome_range(chromosome, minimum, maximum):
    chromosome_range = "Chr{0:02d}:{1}-{2}".format(
        int(chromosome.lower().replace("chr", "")),
        minimum,
        maximum)
    return chromosome_range


def find_overlapping_features(snp_groups, gff_db, gene_info, feature_type, outfile):
    """For each SNP group, search the gff3 (gene) file for any features that
    overlap the range. Find the feature in the info file, and calculate the
    distance from the actual SNP"""
    with open(outfile, 'w') as output:
        output.write("SNP_GROUP\tCHR\tSNP_POS\tGENE_ID\t" +
                     "GENE_Start\tGENE_End\tGeneINFO\n")
        for key, value in tqdm(snp_groups.items()):
            chromosome, position, minimum, maximum = value
            region = coordinates_to_chromosome_range(chromosome, minimum, maximum)
            features = []
            try:
                for line in gff_db.region(region, featuretype=feature_type):
                    # sys.stderr.write(str(line) + "\n")
                    gene_id = ".".join(line.id.split('.')[0:2])
                    gene_start = line.start
                    gene_end = line.end
                    for feature in gene_info[gene_id]:
                        if feature not in features:
                            output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                                key,
                                chromosome,
                                position,
                                gene_id,
                                gene_start,
                                gene_end,
                                ",".join(feature.split()[3:])
                            )
                            )
                        features.append(feature)
            except ValueError:
                continue


if __name__ == "__main__":
    args = parse_arguments()
    input_file = file_to_list(args.input)
    radial_distance = int(args.distance * 1000)
    input_format = args.input_type.lower()
    if input_format == "gapit":
        variant_positions = get_significant_snps(input_file, args.fdr)
    elif input_format in ["coordinate", "range"]:
        variant_positions = input_file[1:]
    else:
        sys.stderr.write(f"Input type {args.input_type} not recognized.\n")
        sys.exit(1)

    snp_ranges = generate_snp_range(variant_positions, input_format, radial_distance)

    gff_database = open_gff_db(args.gff)
    gene_information = get_gene_info_dictionary(args.info)

    find_overlapping_features(snp_ranges,
                              gff_database,
                              gene_information,
                              args.feature_type,
                              args.output)

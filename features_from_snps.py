#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
By J. Lucas Boatwright 

Find any genes near SNPs found to be significant in a GWAS
"""

import sys
import argparse
import gffutils
from collections import OrderedDict

def parse_arguments():
    """Parse arguments passed to script"""
    parser = argparse.ArgumentParser(description="This script was " + 
            "originally designed to find any genes near SNPs found " +
            "to be significant in a GWAS analysis performed in GAPIT " + 
            "but has since been generalize for any coordinates.\n\n")

    requiredNamed = parser.add_argument_group('required arguments')

    requiredNamed.add_argument(
            "--input", 
            type=str, 
            required=True, 
            help="The name of the input file (CSV).",
            action="store")

    requiredNamed.add_argument(
            "--input_type",
            type=str,
            required=True,
            help="The type of input file: gapit, coordinate, range." + 
            " Note: headers are **required** - or first line lost." +
            " *gapit - the unedited GAPIT output with header*" +
            " *coordinate - markerName,chromosome,position ***" +
            " *range - markerName,chromosome,position,start,end*",
            action="store")

    requiredNamed.add_argument(
            "--output", 
            type=str, 
            required=True, 
            help="The output file to be created (TSV).", 
            action="store")

    parser.add_argument(
            "--gff", 
            type=str, 
            required=False, 
            default="/zfs/tillers/Reference_Genomes/BTx623/v3.1.1/" + 
            "annotation/Sbicolor_454_v3.1.1.gene.gff3.sqlite3",
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
            help="The significance threshold for significant SNPs. Only " +
            "applicable with input_type gapit.", 
            action="store")

    parser.add_argument(
            "--info", 
            type=str, 
            required=False, 
            default="/zfs/tillers/Reference_Genomes/BTx623/v3.1.1/" + 
            "annotation/Sbicolor_454_v3.1.1.annotation_info.txt",
            help="The gene info file as obtained from Phytozome (TSV).", 
            action="store")

    parser.add_argument(
            "--feature_type",
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
    except ValueError:
        db = gffutils.create_db(gff, gff + ".sqlite3")
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
        #snp_group = x + 1
        split_line = line.split(',')
        snp_group = split_line[0]
        chromosome = split_line[1]
        position = int(split_line[2])
        if input_type == "range":
            snp_groups[snp_group] = [
                    chromosome,
                    position,
                    int(split_line[3]),
                    int(split_line[4]),
                    1]
        elif input_type in ["gapit", "coordinate"]:
            snp_groups[snp_group] = [
                    chromosome,
                    position,
                    position - distance,
                    position + distance,
                    1]
    return snp_groups


def find_overlapping_features(snp_groups, gff_db, gene_info, feature_type, outfile):
    """For each SNP group, search the gff3 (gene) file for any features that 
    overlap the range. Find the feature in the info file, and calculate the 
    distance from the actual SNP"""
    with open(outfile, 'w') as output:
        output.write("SNP_GROUP\tCHR\tSNP_POS\tNUM_SNPs\tGENE_ID\t" + 
                "GENE_Start\tGENE_End\tGeneINFO\n")
        for key, value in snp_groups.items():
            chromosome, position, minimum, maximum, snp_count = value
            region = "Chr{0:02d}:{1}-{2}".format(int(chromosome.lower().replace("chr","")), 
                    minimum, 
                    maximum) 
            features = []
            for line in gff_db.region(region, featuretype=feature_type):
                # sys.stderr.write(str(line) + "\n")
                gene_id = ".".join(line.id.split('.')[0:2])
                gene_start = line.start
                gene_end = line.end
                for feature in gene_info[gene_id]:
                    if feature not in features:
                        output.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(key, 
                            chromosome, 
                            position, 
                            snp_count, 
                            gene_id, 
                            gene_start,
                            gene_end,
                            ",".join(feature.split()[3:])
                            )
                            )
                    features.append(feature)


def main(args):
    """Main exection function"""
    input_file = file_to_list(args.input)
    distance   = int(args.distance * 1000)
    input_type = (args.input_type).lower()
    if input_type == "gapit":
        sig_gwas = get_significant_snps(input_file, args.fdr)
    elif input_type in ["coordinate", "range"]:
        sig_gwas = input_file[1:]
    else:
        sys.stderr.write("Input type {0} not recognized.".format(args.input_type))
        sys.exit(1)

    snp_groups = generate_snp_range(sig_gwas, input_type, distance)

    gff_db     = open_gff_db(args.gff)
    info       = get_gene_info_dictionary(args.info)

    find_overlapping_features(snp_groups, 
            gff_db, 
            info, 
            args.feature_type, 
            args.output)


if __name__ == "__main__":
    args = parse_arguments()
    main(args)


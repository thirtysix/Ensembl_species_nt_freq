#!/usr/bin/python
# -*- coding: utf-8 -*-
# Python vers. 3.8.0 ###########################################################
# Libraries ####################################################################
import os
import pandas as pd
import pathlib

import sys
import ftplib
from ftplib import FTP

import gzip
from Bio import SeqIO

import time
################################################################################
# Description/Notes ############################################################
################################################################################
"""

"""

################################################################################
# Base-level Functions #########################################################
################################################################################
def load_json(object_filename):
    """Load a json file to object."""

    import os
    import json
    

    if os.path.exists(object_filename):
        with open(object_filename, 'r') as object_file:
            return json.load(object_file)

def save_json(json_dict, json_filename):
    """Save json object to file."""
    import json
    
    with open(json_filename, 'w') as json_file:
        json.dump(json_dict, json_file)

################################################################################
# Task-specific Functions ######################################################
################################################################################


################################################################################
# Initiating Variables #########################################################
################################################################################
# name - data dn
data_dn = "data"
# name - genome dn
genome_dn = os.path.join(data_dn, "Ensembl_genomes")
# make - data/genome dir
pathlib.Path(genome_dn).mkdir(parents=True, exist_ok=True)

# name - species nt count fn
all_species_nt_count_d_fn = "species_nt_count.json"
# name - species nt freq fn
all_species_nt_freq_d_fn = "species_nt_freq.json"


################################################################################
# Execution ####################################################################
################################################################################

ftp=FTP("ftp.ensembl.org")
ftp.login("anonymous","")
ftp_target_dn = "pub/release-110/fasta/"
species_dns = ftp.nlst(ftp_target_dn)

all_species_nt_count_d = {}
all_species_nt_count_ls = []

if os.path.exists(all_species_nt_count_d_fn):
    all_species_nt_count_d = load_json(all_species_nt_count_d_fn)

# iterate - ftp directories within ftp folder
for species_dn in species_dns:
    start_time = time.time()

    # initiate - dict to hold nt counts for this genome
    species_nt_counts_d = {"A":0, "C":0, "G":0, "T":0}

    # name - current species
    species = os.path.basename(species_dn)
    print(species)

    # bool - check for calculated nt counts for species in all_species_nt_count_d
    if species not in all_species_nt_count_d:

        # set - empty genome filename
        species_genome_fn = ""

        # set - species out dir
        species_odn = os.path.join(genome_dn, species)

        # make - species out dir
        pathlib.Path(species_odn).mkdir(parents=True, exist_ok=True)

        # get - files in species dir
        species_fns = [os.path.join(species_odn, x) for x in os.listdir(species_odn) if os.path.isfile(os.path.join(species_odn, x))]

        # bool
        for species_fn in species_fns:
            if species.lower() == os.path.basename(species_fn).split(".")[0].lower() and os.path.getsize(species_fn)>0:
                species_genome_fn = species_fn

        if species_genome_fn == "":
        
            # name - ftp dir species + dna
            dna_dn = os.path.join(species_dn, "dna")

            # get - contents of ftp dir species + dna
            dna_dn_contents = ftp.nlst(dna_dn)

            # set - fn to empty
            dna_rm_fn = ""

            # iterate - contents of ftp dir species + dna
            for fn in dna_dn_contents:
            
                # bool - found dna_rm fasta file
                if ".dna_rm.toplevel.fa.gz" in fn:
                    dna_rm_fn = fn

            if dna_rm_fn == "":
                print(dna_rm_fn, "not found")

            else:
                # name - ftp file basename
                out_basename = os.path.basename(dna_rm_fn)
                # set - local out filename
                species_genome_fn = os.path.join(species_odn, out_basename)

                # bool - if out fn does not exist
                if not os.path.exists(species_genome_fn):
                    print("downloading:", dna_rm_fn)

                    # save - download and save local file 
                    with open(species_genome_fn, 'wb') as out_h:
                        ftp.retrbinary("RETR " + dna_rm_fn, out_h.write, 1024)

        # bool - if species_genome_fn does exist
        if os.path.exists(species_genome_fn) and os.path.getsize(species_genome_fn)>0:
            print("file exists:", species_genome_fn)

            # load - genome file
            with gzip.open(species_genome_fn, "rt") as handle:
                for record in SeqIO.parse(handle, "fasta"):
                    for nt in species_nt_counts_d.keys():
                        species_nt_counts_d[nt] += record.seq.count(nt)
            
            # copy - species_nt_counts_d
            species_nt_counts_d_dup = species_nt_counts_d.copy()

            # calculate - add reverse strand counts
            species_nt_counts_d['A']+=species_nt_counts_d_dup['T']
            species_nt_counts_d['T']+=species_nt_counts_d_dup['A']
            species_nt_counts_d['C']+=species_nt_counts_d_dup['G']
            species_nt_counts_d['G']+=species_nt_counts_d_dup['C']

            # add - species genome nt counts d to all species'
            all_species_nt_count_d[species] = species_nt_counts_d

            # save - all species' json d to file
            save_json(all_species_nt_count_d, all_species_nt_count_d_fn)

    # bool - species nt count data already in all_species_nt_count_d
    if species in all_species_nt_count_d:
        # get - species nt count data
        species_nt_counts_d = all_species_nt_count_d[species]

    print(species_nt_counts_d)
    
    # initiate - entry ls for species
    species_entry_ls = [species] + [species_nt_counts_d[nt] for nt in 'ACGT']

    # calculate - sum of all nucleotides
    nt_sum = sum(list(species_nt_counts_d.values()))

    # set - starting ls nt frequency
    nt_freqs = [0,0,0,0]

    # bool - check if sum greater than zero
    if nt_sum >0:
    
        # calculate - nt frequencies for genome
        for i, nt in enumerate("ACGT"):
            nt_count = species_nt_counts_d[nt]
            nt_freq = nt_count/nt_sum
            nt_freqs[i] = nt_freq
            #species_entry_ls.append(nt_freq)

    # add - species data to running list
    all_species_nt_count_ls.append(species_entry_ls+nt_freqs)

    print("time:", round(time.time() - start_time, 3))

# calculate - nt freqs for all species; as dict
all_species_nt_freq_d = {species:{nt:count/sum(nt_count_d.values()) for nt,count in nt_count_d.items()} for species, nt_count_d in all_species_nt_count_d.items()}

# output - nt freqs for all species dict
save_json(all_species_nt_freq_d, all_species_nt_freq_d_fn)

# initiate - df of species data
all_species_nt_count_df = pd.DataFrame(all_species_nt_count_ls, columns = ["species", "A_count", "C_count", "G_count", "T_count", "A_freq", "C_freq", "G_freq", "T_freq"])

# output - species data
all_species_nt_count_df.to_csv("Ensembl_species_nt_data.tsv", sep="\t", index=False)






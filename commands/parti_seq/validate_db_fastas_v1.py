#!/usr/bin/python

import os
import sys

import modules as m


cmd = m.Command(__file__)

data_dpath = os.path.join(cmd.dpath, "data")

db_fpath = "/media/sl/T7-Shield/PaRTISeq/blast_db"

fasta_json = {
  # Human Genome Project (HGP)
  "human_hgp": {
    "fname": "human.fa"
  },
  # Telomere-to-Telomere (T2T) Consortium
  "human_t2t": {
    "fname": "CHM13v2.0.fa",
  },

  "step1_db": {
    "fname": "Step1_ref_230828.fa",
  },

  "step2_db": {
    "fname": "Step2_ref_230902.fa",
  },

  "step3_db": {
    "fname": "Step3_ref_230901.fa"
  }
}


def validate_db_fastas(arg_fasta_fpath=None):

  deflines_fpath = os.path.join(data_dpath, "deflines.tsv")
  deflines_summary_fpath = os.path.join(data_dpath, "deflines_summary.tsv")
  dup_deflines_fpath = os.path.join(data_dpath, "dup_deflines.tsv")
  dup_defline_ids_fpath = os.path.join(data_dpath, "dup_defline_ids.tsv")
  dup_defline_ids_fa_fpath = os.path.join(data_dpath, "dup_defline_ids.fa")
  id_map_fpath = os.path.join(db_fpath, "ID_map_230511.txt")
  plasmid_id_map_fpath = os.path.join(db_fpath, "ID_ref_plasmid_map.txt")
  dup_ids_in_map_fpath = os.path.join(data_dpath, "dup_ids_in_map.tsv")
  missing_ids_fpath = os.path.join(data_dpath, "missing_ids.txt")
  missing_ids_summary_fpath = os.path.join(data_dpath, "missing_ids_summary.tsv")


  ##############################################################################
  # Put all deflines from all fastas into a single file

  print("Scanning for deflines...")
  f_out = open(deflines_fpath, 'w')
  for db_name, db in fasta_json.items():
    fasta_fname = db["fname"]
    fasta_fpath = os.path.join(db_fpath, fasta_fname)

    with open(fasta_fpath, 'r') as f_in:
      print(f" Reading {fasta_fname}...")
      defline_qty = 0
      for i, line in enumerate(f_in):
        sys.stdout.write("\r")
        if line.startswith(">"):
          defline_qty += 1
          f_out.write(f"{db_name}\t{line}")
        if (i + 1) % 1_000_000 == 0:
          sys.stdout.write(f"  {i+1} lines read")
      db["line_qty"] = i + 1
      db["defline_qty"] = defline_qty
      sys.stdout.write(f"  {i+1} lines read")
      print(f"\n  defline_qty: {defline_qty}")
  f_out.close()

  # Write defline summary to file
  total_defline_qty = 0
  with open(deflines_summary_fpath, 'w') as f:
    for db_name, db in fasta_json.items():
      total_defline_qty += db['defline_qty']
      f.write(
        f"{db_name}:\n"
        f"  line_qty: {db['line_qty']}\n"
        f"  defline_qty: {db['defline_qty']}\n"
      )
    f.write(f"total_defline_qty: {total_defline_qty}\n")
  ##############################################################################

  
  ##############################################################################
  # Create a file for duplicate deflines and a file for duplicate IDs

  deflines = set()
  dup_deflines = set()
  dup_defline_qty = 0
  defline_ids = set()
  dup_defline_ids = set()
  dup_defline_id_qty = 0

  print("Scanning for duplicate deflines and IDs...")
  with open(deflines_fpath, 'r') as f_in:

    for i, line in enumerate(f_in):
      sys.stdout.write("\r")
      ar_line = line.strip("\n").split("\t")
      
      # Check defline
      db_name = ar_line[0]
      defline = ar_line[1]
      if defline in deflines:
        dup_defline_qty += 1
        dup_deflines.add(defline)
      else:
        deflines.add(defline)

      # Check ID in defline
      ar_defline = defline.split(" ")
      defline_id = ar_defline[0]
      if defline_id in defline_ids:
        dup_defline_id_qty += 1
        dup_defline_ids.add(defline_id)
      else:
        defline_ids.add(defline_id)
      if (i+1) % 1_000_000 == 0:
        sys.stdout.write(f"  {i+1} lines read")
    print(f"\n Duplicate defline qty: {dup_defline_qty}")
    print(f" Duplicate IDs in deflines qty: {dup_defline_id_qty}")
  
  # Write duplicate deflines to file
  print("Writing duplicate deflines and IDs to file...")
  with open(deflines_fpath, 'r') as f_in, \
    open(dup_deflines_fpath, 'w') as f_dup_deflines, \
    open(dup_defline_ids_fpath, 'w') as f_dup_defline_ids:
    
    for i, line in enumerate(f_in):
      sys.stdout.write("\r")
      ar_line = line.strip("\n").split("\t")

      # Check defline
      db_name = ar_line[0]
      defline = ar_line[1]
      if defline in dup_deflines:
        f_dup_deflines.write(line)

      # Check ID in defline
      ar_defline = defline.split(" ")
      defline_id = ar_defline[0]
      if defline_id in dup_defline_ids:
        f_dup_defline_ids.write(line)
      
      sys.stdout.write(f" {i+1} lines read")
    print("\n Done writing lines.") 
  ##############################################################################
  

  ##############################################################################
  # Write fasta of duplicate defline IDs to file      
  # This will account for deuplicate deflines as well

  print("Getting sequences from duplicate IDs in deflines...")
  
  # Get the saved duplicate deflines
  dup_deflines = set()
  with open(dup_defline_ids_fpath, 'r') as f:
    for line in f:
      ar_line = line.strip("\n").split("\t")
      defline = ar_line[1]
      dup_deflines.add(defline)
  
  # Scan for duplicate deflines
  # This is extra thurough in case anything was missed
  f_out = open(dup_defline_ids_fa_fpath, 'w')
  for db_name, db in fasta_json.items():
    fasta_fname = db["fname"]
    fasta_fpath = os.path.join(db_fpath, fasta_fname)

    with open(fasta_fpath, 'r') as f_in:
      
      print(f" Reading {fasta_fname}...")
      is_duplicate = False
      
      for i, line in enumerate(f_in):
        sys.stdout.write("\r")
        
        if line.startswith(">") and (line.strip("\n") in dup_deflines):
          is_duplicate = True
          f_out.write(f">{db_name}{line}")

        elif line.startswith(">") and (line.strip("\n") not in dup_deflines):
          is_duplicate = False
        
        elif is_duplicate:
          f_out.write(line)
          
        if (i+1) % 1_000_000 == 0:
          sys.stdout.write(f"  {i+1} lines read")
      
      sys.stdout.write(f"  {i+1} lines read")
      print()
  f_out.close()
  ##############################################################################


  ##############################################################################
  # Check for duplicates in the ID Map
  
  id_map_ids = set()
  dup_id_map_ids = set()
  with open(id_map_fpath, 'r') as f:
    for line in f:
      ar_line = line.strip("\n").split("\t")
      id = ar_line[0]
      if id in id_map_ids:
        dup_id_map_ids.add(id)
      else:
        id_map_ids.add(id)
  print(f"IDs in Map qty: {len(id_map_ids)}")
  
  # Get the duplicate lines
  dup_id_map_id_lines = []
  with open(id_map_fpath, 'r') as f:
    for line in f:
      ar_line = line.strip("\n").split("\t")
      id = ar_line[0]
      if id in dup_id_map_ids:
        dup_id_map_id_lines.append(line)
  print(f"Duplicate IDs in Map qty: {len(dup_id_map_id_lines)}")

  # Write lines with duplicate IDs to file
  dup_id_map_id_lines.sort()
  with open(dup_ids_in_map_fpath, 'w') as f:
    for line in dup_id_map_id_lines:
      f.write(line)
  ##############################################################################


  ##############################################################################
  # Check if there are IDs extra or missing between the ID Map file and
  # the fasta files.

  f_summary = open(missing_ids_summary_fpath, 'w')

  id_map_ids = set()
  with open(id_map_fpath, 'r') as f:
    for line in f:
      ar_line = line.strip("\n").split("\t")
      id = ar_line[0]
      id_map_ids.add(id)
  print(f"IDs in Map qty: {len(id_map_ids)}")
  f_summary.write(f"IDs in Map qty: {len(id_map_ids)}\n")

  defline_ids = set()
  with open(deflines_fpath, 'r') as f:
    for line in f:
      ar_line = line.strip("\n").split("\t")
      defline = ar_line[1]
      ar_defline = defline.split(" ")
      id = ar_defline[0][1:]
      defline_ids.add(id)
  
  fasta_missing_ids = set()
  id_map_missing_ids = set()
  with open(missing_ids_fpath, 'w') as f:
    
    f.write("IDs missing from fastas:\n")
    for id in id_map_ids:
      if id not in defline_ids:
        fasta_missing_ids.add(id)
        f.write(f"{id}\n")
    if not len(fasta_missing_ids) > 0:
      f.write("none\n")
    
    f.write("IDs missing from ID map:\n")
    for id in defline_ids:
      if id not in id_map_ids:
        id_map_missing_ids.add(id)
        f.write(f"{id}\n")
    if not len(id_map_missing_ids) > 0:
      f.write("none\n")
  
  print(f"IDs missing from fastas: {len(fasta_missing_ids)}")
  f_summary.write(f"IDs missing from fastas: {len(fasta_missing_ids)}\n")
  print(f"IDs missing from ID map: {len(id_map_missing_ids)}")
  f_summary.write(f"IDs missing from ID map: {len(id_map_missing_ids)}\n")

  # Check if the IDs missing from the ID map are in the plasmid ID map
  plasmid_id_map_ids = set()
  with open(plasmid_id_map_fpath, 'r') as f:
    for line in f:
      ar_line = line.strip("\n").split("\t")
      id = ar_line[0]
      plasmid_id_map_ids.add(id)
  print(f"IDs in plasmid Map qty: {len(plasmid_id_map_ids)}")
  f_summary.write(f"IDs in plasmid Map qty: {len(plasmid_id_map_ids)}\n")

  fasta_missing_plasmid_ids = set()
  plasmid_id_map_missing_ids = set()
  
  for id in plasmid_id_map_ids:
    if id not in defline_ids:
      fasta_missing_plasmid_ids.add(id)
  
  for id in id_map_missing_ids:
    if id not in plasmid_id_map_ids:
      plasmid_id_map_missing_ids.add(id)
  
  print(f"Plasmid IDs missing from fastas: {len(fasta_missing_plasmid_ids)}")
  f_summary.write(f"Plasmid IDs missing from fastas: {len(fasta_missing_plasmid_ids)}\n")
  print(f"IDs missing from Plasmid ID map: {len(plasmid_id_map_missing_ids)}")
  f_summary.write(f"IDs missing from Plasmid ID map: {len(plasmid_id_map_missing_ids)}\n")

  f_summary.close()


  ##############################################################################


  # # Get ids and warn if ID is already in set 

  # fasta_fpath = m.Fasta(arg_fasta_fpath)
  
  # ids = []
  # duplicate_ids = []
  # with open(fasta_fpath, 'r') as f:
  #   for line in f:
  #     if line.startswith(">"):
  #       id = line[1:].strip("\n")
  #       if id in ids and id not in duplicate_ids:
  #         duplicate_ids.append(id)
  #       ids.append(id)
  # if len(duplicate_ids) > 0:
  #   print(f"Warning: Duplicate IDs found in fasta file.")
  
  # # with open("get_lineages.duplicate_ids.tsv", 'w') as f:
  # #   for duplicate_id in duplicate_ids:
  # #     f.write(f"{duplicate_id}\n")

  # # with open("get_lineages.duplicates.fasta", 'w') as f_out:
  # #   print("Copying duplicates:")
  # #   for duplicate_id in duplicate_ids:
  # #     print(f" {duplicate_id}")
  # #     with open(fasta.fpath, 'r') as f_in:
  # #       is_duplicate = False
  # #       for line in f_in:
  # #         if line.startswith(">") and duplicate_id in line:
  # #           is_duplicate = True
  # #         elif line.startswith(">") and duplicate_id not in line:
  # #           is_duplicate = False
  # #         if is_duplicate:
  # #           f_out.write(line)


  # lineages = {
  #   # ids: [
  #   #   species # from ID_map_230511
  #   # ]
  # }
  
  # id_to_species_map = []
  # missing_ids = [] # IDs that are in the ID map but not in the fasta
  # species_list = []
  # taxon_id_map_fpath = "/media/sl/T7-Shield/PaRTISeq/blast_db/ID_map_230511.txt"
  # print("Reading ID map...")
  # with open(taxon_id_map_fpath, 'r') as f:
  #   for line in f:
  #     ar_line = line.strip("\n").split("\t")
  #     id = ar_line[0]
  #     species = ar_line[1].replace("_", " ")
      
  #     if species not in species_list:
  #       species_list.append(species)
  #     i_species_list = species_list.index(species)
      
  #     if id in ids:
  #       i_id = ids.index(id)
  #     else:
  #       i_id = None
  #       missing_ids.append([id, species])
      
  #     id_to_species_map.append(
  #       [i_id, i_species_list]
  #     )

  #     sys.stdout.write(f" {len(id_to_species_map)} lines read\r")
  #   print(f" {len(id_to_species_map)} lines read")

  # if len(missing_ids) > 0:
  #   print(f"Warning: Some IDs are in the ID map but not in the fasta.")
  # with open("get_lineages.missing_ids.tsv", 'w') as f:
  #   for missing_id in missing_ids:
  #     f.write("\t".join(missing_id)+"\n")
      
  # # Dump unique species list for visual check
  # with open("get_lineages.species.txt", 'w') as f:
  #   for species in species_list:
  #     f.write(f"{species}\n")

#!/usr/bin/python

import json
import os
from pathlib import Path
import sys

import modules as m


debug = m.Debugger.devtools_debug

cmd = m.Command(__file__)
data_dpath = cmd.data_dpath

'''
  python -m fire \
    ./commands/parti_seq/validate_db_fastas_v2.py \
    validate_db_fastas \
    --fastas_fpath=/media/sl/T7-Shield/PaRTISeq/blast_db_genomes \
    --fasta_json_fpath=commands/parti_seq/data/fasta-meta.json \
    --db_dpath=/media/sl/T7-Shield/PaRTISeq/blast_db
'''
def validate_db_fastas(
  fastas_fpath,
  fasta_json_fpath,
  db_dpath # eg, /media/sl/T7-Shield/PaRTISeq/blast_db
):
  '''
  {
    "fastas": {
      "human.fa": {
        "comment": "Human Genome Project (HGP)",
        "name": "human_hgp",
        "fname": "human.fa"
      },
      ...
    }
  }
  '''
  
  fasta_json = m.Json.load_json_as_dot_notation(fasta_json_fpath)
  all_deflines_fpath = os.path.join(data_dpath, "all_deflines.tsv")
  dup_deflines_fpath = os.path.join(data_dpath, "dup_deflines.tsv")
  dup_defline_ids_fpath = os.path.join(data_dpath, "dup_defline_ids.tsv")
  dup_defline_ids_fa_fpath = os.path.join(data_dpath, "dup_defline_ids.fa")
  id_map_fpath = os.path.join(db_dpath, "ID_map_230511.txt")
  plasmid_id_map_fpath = os.path.join(db_dpath, "ID_ref_plasmid_map.txt")
  dup_ids_in_map_fpath = os.path.join(data_dpath, "dup_ids_in_map.tsv")
  fasta_missing_ids_fpath = os.path.join(data_dpath, "fasta_missing_ids.txt")
  id_map_missing_ids_fpath = os.path.join(data_dpath, "id_map_missing_ids.txt")
  missing_ids_summary_fpath = os.path.join(data_dpath, "missing_ids_summary.tsv")
  organisms_fpath = os.path.join(data_dpath, "organisms.json")

  fastas = []
  for name in fasta_json.fastas:
    fasta_fname = fasta_json.fastas[name].fname
    fastas.append(m.Fasta(os.path.join(fastas_fpath, fasta_fname)))
  

  # ############################################################################
  # # Scan for and export deflines, this takes a while.
  # # Could include qty of bases in the output in the future.

  # for fasta in fastas:
  #   fasta.export_deflines(dpath=data_dpath)
    
  # defline_fpaths = []
  # for fasta in fastas:
  #   defline_fpaths.append(os.path.join(
  #     data_dpath,
  #     m.File.get_fname_without_extension(fasta.fpath) + "_deflines.tsv"
  #   ))

  # m.Fasta.concat_defline_files(
  #   fpaths = defline_fpaths,
  #   out_fpath = all_deflines_fpath,
  #   keep_originals = True
  # )
    
  # ############################################################################
  
  # ############################################################################
  # # Get the defline qty
    
  # with open(all_deflines_fpath, 'r') as f:
  #   for line in f:
  #     ar_line = line.split("\t")
  #     fasta_fname = ar_line[0]
      
  #     if "defline_qty" not in fasta_json.fastas[fasta_fname]:
  #       fasta_json.fastas[fasta_fname].defline_qty = 1
  #     else: 
  #       fasta_json.fastas[fasta_fname].defline_qty += 1
      
  # # Get total defline qty
  # for fasta_fname in fasta_json.fastas:
  #   this_defline_qty = fasta_json.fastas[fasta_fname].defline_qty
  #   if "defline_qty" not in fasta_json:
  #     fasta_json.defline_qty = this_defline_qty
  #   else:
  #     fasta_json.defline_qty += this_defline_qty

  # ############################################################################
      
  # ############################################################################
  # # Create a file for duplicate deflines and a file for duplicate IDs

  # deflines = set()
  # dup_deflines = set()
  # dup_defline_qty = 0
  # defline_ids = set()
  # dup_defline_ids = set()
  # dup_defline_id_qty = 0

  # print("Scanning for duplicate deflines and IDs...")
  # with open(all_deflines_fpath, 'r') as f_in:

  #   for i, line in enumerate(f_in):
  #     sys.stdout.write("\r")
  #     ar_line = line.strip("\n").split("\t")
      
  #     # Check defline
  #     db_name = ar_line[0]
  #     defline = ar_line[1]
  #     if defline in deflines:
  #       dup_defline_qty += 1
  #       dup_deflines.add(defline)
  #     else:
  #       deflines.add(defline)

  #     # Check ID in defline
  #     ar_defline = defline.split(" ")
  #     defline_id = ar_defline[0]
  #     if defline_id in defline_ids:
  #       dup_defline_id_qty += 1
  #       dup_defline_ids.add(defline_id)
  #     else:
  #       defline_ids.add(defline_id)
  #     if (i + 1) % 1_000_000 == 0:
  #       sys.stdout.write(f"  {i+1} lines read")
  #   sys.stdout.write(f"  {i+1} lines read")
  #   print(f"\n  Duplicate defline qty: {dup_defline_qty}")
  #   print(f"  Duplicate IDs in deflines qty: {dup_defline_id_qty}")
  
  # # Write duplicate deflines to file
  # print("Writing duplicate deflines and IDs to file...")
  # with open(all_deflines_fpath, 'r') as f_in, \
  #   open(dup_deflines_fpath, 'w') as f_dup_deflines, \
  #   open(dup_defline_ids_fpath, 'w') as f_dup_defline_ids:
    
  #   for i, line in enumerate(f_in):
  #     sys.stdout.write("\r")
  #     ar_line = line.strip("\n").split("\t")

  #     # Check defline
  #     db_name = ar_line[0]
  #     defline = ar_line[1]
  #     if defline in dup_deflines:
  #       f_dup_deflines.write(line)

  #     # Check ID in defline
  #     ar_defline = defline.split(" ")
  #     defline_id = ar_defline[0]
  #     if defline_id in dup_defline_ids:
  #       f_dup_defline_ids.write(line)
      
  #     if (i + 1) % 1_000_000 == 0:
  #       sys.stdout.write(f"  {i+1} lines read")
  #   sys.stdout.write(f"  {i+1} lines read")

  #   print()

  # ############################################################################
    
  # ############################################################################
  # # Write fasta of duplicate defline IDs to file      
  # # This will account for deuplicate deflines as well

  # print("Getting sequences from duplicate IDs in deflines...")
  
  # # Get the saved duplicate deflines
  # dup_deflines = set()
  # with open(dup_defline_ids_fpath, 'r') as f:
  #   for line in f:
  #     ar_line = line.strip("\n").split("\t")
  #     defline = ar_line[1]
  #     dup_deflines.add(defline)
  
  # # Scan for duplicate deflines
  # # This is extra thurough in case anything was missed
  # f_out = open(dup_defline_ids_fa_fpath, 'w')
  # for fasta_fname in fasta_json.fastas:
  #   fasta_fpath = os.path.join(fastas_fpath, fasta_fname)

  #   with open(fasta_fpath, 'r') as f_in:
      
  #     print(f" Reading {fasta_fname}...")
  #     is_duplicate = False
      
  #     for i, line in enumerate(f_in):
  #       sys.stdout.write("\r")
        
  #       if line.startswith(">") and (line.strip("\n") in dup_deflines):
  #         is_duplicate = True
  #         f_out.write(f">{db_name}{line}")

  #       elif line.startswith(">") and (line.strip("\n") not in dup_deflines):
  #         is_duplicate = False
        
  #       elif is_duplicate:
  #         f_out.write(line)
          
  #       if (i+1) % 1_000_000 == 0:
  #         sys.stdout.write(f"  {i+1} lines read")
      
  #     sys.stdout.write(f"  {i+1} lines read")
  #     print()
  # f_out.close()

  # ############################################################################
  
  # ############################################################################
  # # Check for duplicates in the ID Map
  
  # id_map_ids = set()
  # dup_id_map_ids = set()
  # with open(id_map_fpath, 'r') as f:
  #   for line in f:
  #     ar_line = line.strip("\n").split("\t")
  #     id = ar_line[0]
  #     if id in id_map_ids:
  #       dup_id_map_ids.add(id)
  #     else:
  #       id_map_ids.add(id)
  # print(f"IDs in Map qty: {len(id_map_ids)}")
  
  # # Get the duplicate lines
  # dup_id_map_id_lines = []
  # with open(id_map_fpath, 'r') as f:
  #   for line in f:
  #     ar_line = line.strip("\n").split("\t")
  #     id = ar_line[0]
  #     if id in dup_id_map_ids:
  #       dup_id_map_id_lines.append(line)
  # print(f"Duplicate IDs in Map qty: {len(dup_id_map_id_lines)}")

  # # Write lines with duplicate IDs to file
  # dup_id_map_id_lines.sort()
  # with open(dup_ids_in_map_fpath, 'w') as f:
  #   for line in dup_id_map_id_lines:
  #     f.write(line)

  # ############################################################################
  
  ############################################################################
  # Check if there are IDs extra or missing between the ID Map file and
  # the fasta files.

  f_summary = open(missing_ids_summary_fpath, 'w')

  id_map_ids = dict()
  with open(id_map_fpath, 'r') as f:
    for line in f:
      ar_line = line.strip("\n").split("\t")
      id = ar_line[0]
      lineage = ar_line[1]
      id_map_ids[id] = lineage
  print(f"IDs in Map qty: {len(id_map_ids)}")
  f_summary.write(f"IDs in Map qty: {len(id_map_ids)}\n")

  defline_ids = dict()
  with open(all_deflines_fpath, 'r') as f:
    for line in f:
      ar_line = line.strip("\n").split("\t")
      defline = ar_line[1]
      ar_defline = defline.split(" ")
      id = ar_defline[0][1:]
      defline_ids[id] = defline
  
  fasta_missing_ids = dict()
  with open(fasta_missing_ids_fpath, 'w') as f:
    for id in id_map_ids:
      if id not in defline_ids:
        fasta_missing_ids[id] = id_map_ids[id]
        f.write(f"{id}\t{id_map_ids[id]}\n")
  
  id_map_missing_ids = dict()
  with open(id_map_missing_ids_fpath, 'w') as f:
    for id in defline_ids:
      if id not in id_map_ids:
        id_map_missing_ids[id] = defline_ids[id]
        f.write(f"{id}\t{defline_ids[id]}\n")
  
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

  ############################################################################

  # ############################################################################
  # # Get organisms, typically genus and species

  # organisms = []
  # organism_name_to_organisms_idx = dict()
  # for id in defline_ids:
    
  #   if id not in id_map_ids:
  #     continue
    
  #   organism = id_map_ids[id].replace("_", " ")
    
  #   if organism not in organism_name_to_organisms_idx:
  #     organisms.append({
  #       'name': organism,
  #       'sequence_ids': [id]
  #     })
  #     organism_name_to_organisms_idx[organism] = len(organisms) - 1
    
  #   else:
  #     idx = organism_name_to_organisms_idx[organism]
  #     organisms[idx]['sequence_ids'].append(id)
      
  # m.Json.save_to_file(organisms, organisms_fpath)

  # ############################################################################

  ############################################################################
  # Get lineages from organism data
  
  organisms = m.Json.load_json_file(organisms_fpath)
  organism_qty = len(organisms)
  print("Requesting taxon IDs from NIH...")
  for i, organism in enumerate(organisms):
    sys.stdout.write("\r")
    name = organism['name']
    taxon_id = m.ncbi.Entrez.Taxonomy.get_taxon_id(query=name)
    # debug(esearch_json)
    organisms[i]['nih_taxon_id'] = taxon_id
    sys.stdout.write(f"  {i+1} of {organism_qty} retrieved.")
      
  print(f"  {i+1} of {organism_qty} retrieved.")
  
  m.Json.save_to_file(organisms, organisms_fpath)

  ############################################################################


  # Write meta data to file      
  with open(os.path.join(data_dpath, "fasta-meta-out.json"), 'w') as f:
    f.write(json.dumps(fasta_json, indent=2))
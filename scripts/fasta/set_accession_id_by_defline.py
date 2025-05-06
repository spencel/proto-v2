
import os

from modules import json as m_json
from modules import file_sys


def set_accession_id_by_defline(
	deflines_fpath: str,
	taxons_fpath: str,
	defline_idx: int = 0,
	has_header: bool = False,
	has_accession: bool = True,
	out_fpath: str|None = None,
	unresolved_fpath: str|None = None
) -> int:
	"""_summary_

	Args:
			fpath (str): _description_
			defline_idx (int, optional): Index of the defline in the input table. Defaults to 0.
			has_header (bool, optional): Whether the file has a header row. Defaults to False.
			has_accession (bool, optional): Whether the defline starts with an accession. Defaults to True.
			out_fpath (str | None, optional): _description_. Defaults to None.

	Returns:
			int: Number of deflines that couldn't be resolved.
	"""
	# Create list of taxons and sort by genus-species name length
	taxons: list = []
	with open(taxons_fpath) as f:
		next(f)
		for line in f:
			ar_line = line.strip('\n').split('\t')
			species_name = ar_line[0]
			taxon_id = ar_line[1]
			taxons.append([len(species_name), species_name, taxon_id])

	# Sort taxons by genus-species name and alphabetically
	taxons = sorted(taxons, key=lambda x: (-x[0], x[1]))
	m_json.save_to_file(taxons, 'logs/taxons.json')

	# Get list of deflines and sort it
	deflines: list = []
	with open(deflines_fpath) as f:
		if has_header: next(f)
		for line in f:
			ar_line = line.strip('\n').split('\t')
			defline = ar_line[defline_idx]
			accession_id: str| None = None
			if has_accession: 
				accession_id = defline[1:].split(' ')[0]
				defline = defline[len(f'>{accession_id} '):]
			deflines.append([defline, accession_id])
	deflines = sorted(deflines)
	m_json.save_to_file(deflines, 'logs/deflines.json')
	
	if not out_fpath:
		out_fpath = file_sys.File.get_fpath_without_extension(
			deflines_fpath
		) + '-taxon-ids.tsv'

	with open(out_fpath, 'w') as f:
		for taxon in taxons:
			
			# End search if there are no more deflines to scan
			if len(deflines) == 0:
				break
			
			name_length = taxon[0]
			species_name = taxon[1]
			taxon_id = taxon[2]

			# Scan the defline list for matches and pop it
			i = 0
			while i < len(deflines):
				ar_defline = deflines[i]
				defline = ar_defline[0]
				accession_id = ar_defline[1]
				if species_name == defline[:name_length]:
					f.write('\t'.join([str(taxon_id), accession_id]) + '\n')
					del deflines[i]
				else:
					i += 1
	
	if not unresolved_fpath:
		unresolved_fpath = file_sys.File.get_fpath_without_extension(
			deflines_fpath
		) + '-unresolved-taxons.tsv'
	
	with open(unresolved_fpath, 'w') as f:
		for defline in deflines:
			f.write('\t'.join(defline) + '\n')
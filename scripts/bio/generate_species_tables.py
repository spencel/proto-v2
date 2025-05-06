
import os

from modules import file_sys


NIH_TAXON_ID_COL_NAME = 'nih_taxon_id'
NIH_GI_COL_NAME = 'nih_gi'
# This value is split into 2 separate columns later
NIH_ACCESSION_VERSION_NAME= 'nih_accession_version'
NIH_TAXON_NAME_COL_NAME = 'taxon_name'
# Column names in input files related to column names in outputfiles
COLUMN_NAMES = {
	'taxid': NIH_TAXON_ID_COL_NAME,
	'gi': NIH_GI_COL_NAME,
	'accessionversion': NIH_ACCESSION_VERSION_NAME,
	'organism': NIH_TAXON_NAME_COL_NAME
}
TAXONS_FILENAME = 'taxon.tsv'
TAXON_NAMES_FILENAME = 'taxon-name.tsv'
NIH_GENBANK_ACCESSIONS_FILENAME = 'nih-genbank-accession.tsv'
DUPLICATE_ACCESSIONS_FILENAME = 'duplicate-accessions.tsv'


def generate_species_tables(
	fpaths: list[str] = [
		'/media/sl/ExtremePro/Parti-Seq/db_nanopore_exports/human-export-accession-ids-esummaries.tsv',
		'/media/sl/ExtremePro/Parti-Seq/db_nanopore_exports/CHM13v2.0-export-accession-ids-esummaries.tsv',
		'/media/sl/ExtremePro/Parti-Seq/db_nanopore_exports/GRCh38p14-export-accession-ids-esummaries.tsv',
		'/media/sl/ExtremePro/Parti-Seq/db_nanopore_exports/Step1_ref_231215-export-accession-ids-esummaries.tsv',
		'/media/sl/ExtremePro/Parti-Seq/db_nanopore_exports/Step2_ref_240115-export-accession-ids-esummaries.tsv',
		'/media/sl/ExtremePro/Parti-Seq/db_nanopore_exports/Step3_ref_231125-export-accession-ids-esummaries.tsv'
	],
	out_dpath: str = 'data/scripts/bio/generate_species_tables'
):
	
	# Generate output dirpath
	file_sys.make_dirs(out_dpath)

################################################################################
# Data structure 
	
	next_taxon_id = 0
	taxons = {}
	'''
		taxons = {
			'NIH_TAXON_ID': {
				'id': int
			},
			...
		}
	'''

	# There can be alternative names, eg, basionym, common name, etc., for the
	# same taxon ID.
	next_taxon_name_id = 0
	taxon_names = {}
	'''
		taxon_names = {
			'NIH_TAXON_NAME': {
				'id': int
				'taxon_id': int -> taxon.id
			},
			...
		}
	'''

	# Stores accessions
	# Split accession and version to separate columns in file
	next_nih_genbank_accession_id = 0
	nih_genbank_accessions	= {}
	'''
		nih_genbank_accessions	= {
			'NIH_ACCESSION_VERSION': {
				'id': int
				'nih_gi': int
				'taxon_id': int -> taxon.id
			},
			...
		}
	'''

################################################################################
# Generate data

	# Track duplicate accessions
	total_duplicate_accessions = 0
	f_duplicates = open(os.path.join(out_dpath, DUPLICATE_ACCESSIONS_FILENAME), 'w')

	# Get unique taxon IDs and make file containing unique taxon IDs
	total_lines_read = 0
	taxons_fpath = os.path.join(out_dpath, 'taxon.tsv')
	for fpath in fpaths:
		with open(fpath) as f_in:
			
			# Get column names and indices
			line = f_in.readline()
			ar_line = line.rstrip('\n').split('\t')
			col_map = {}
			for idx, col_name in enumerate(ar_line):
				if col_name in COLUMN_NAMES:
					col_map[COLUMN_NAMES[col_name]] = idx
				else:
					raise Exception(f'Column name {col_name} not recognized.')
			
			# Read rest of lines
			this_taxon_id: int
			this_taxon_name_id: int
			this_nih_genbank_accession_id: int
			for line in f_in:
				ar_line = line.rstrip('\n').split('\t')
				nih_taxon_id = ar_line[col_map[NIH_TAXON_ID_COL_NAME]]
				nih_gi = ar_line[col_map[NIH_GI_COL_NAME]]
				accession_version = ar_line[col_map[NIH_ACCESSION_VERSION_NAME]]
				taxon_name = ar_line[col_map[NIH_TAXON_NAME_COL_NAME]]
				# print(f'nih_taxon_id: {nih_taxon_id}')
				# print(f'gi: {gi}')
				# print(f'accession_version: {accession_version}')
				# print(f'taxon_name: {taxon_name}')
				# break

				# Handle NIH taxon ID
				if nih_taxon_id not in taxons:
					this_taxon_id = next_taxon_id
					taxons[nih_taxon_id] = {'id': this_taxon_id}
					next_taxon_id += 1
				else:
					this_taxon_id = taxons[nih_taxon_id]['id']

				# Handle taxon name
				if taxon_name not in taxon_names:
					this_taxon_name_id = next_taxon_name_id
					taxon_names[taxon_name] = {
						'id': this_taxon_name_id,
						'taxon_id': this_taxon_id # -> taxon.id
					}
					next_taxon_name_id += 1
				
				# Handle accession version
				if accession_version not in nih_genbank_accessions:
					this_nih_genbank_accession_id = next_nih_genbank_accession_id
					nih_genbank_accessions[accession_version] = {
						'id': this_nih_genbank_accession_id,
						'nih_gi': nih_gi,
						'taxon_id': this_taxon_id # -> taxon.id
					}
					next_nih_genbank_accession_id += 1
				else:
					# Track duplicate accessions
					f_duplicates.write(
						'\t'.join([accession_version, taxon_name, nih_gi, nih_taxon_id]) + '\n'
					)
					total_duplicate_accessions += 1
				
				total_lines_read += 1
	
	f_duplicates.close()

	# Compare total lines read to the lines in NIH_GENBANK_ACCESSIONS_FILENAME
	print(f'total_lines_read: {total_lines_read}')
	print(f'total_duplicate_accessions: {total_duplicate_accessions}')

################################################################################
# Write data to files

	# Handle taxons
	with open(os.path.join(out_dpath, TAXONS_FILENAME), 'w') as f:
		for nih_taxon_id, item in taxons.items():
			id = item['id']
			f.write('\t'.join(map(str, [id, nih_taxon_id])) + '\n')
	
	# Handle taxon names
	with open(os.path.join(out_dpath, TAXON_NAMES_FILENAME), 'w') as f:
		for taxon_name, item in taxon_names.items():
			id = item['id']
			taxon_id = item['taxon_id'] # -> taxon.id
			f.write('\t'.join(map(str, [id, taxon_name, taxon_id])) + '\n')
	
	# Handle accession versions
	with open(os.path.join(out_dpath, NIH_GENBANK_ACCESSIONS_FILENAME), 'w') as f:
		for accession_version, item in nih_genbank_accessions.items():
			id = item['id']
			ar_accession_version = accession_version.split('.')
			accession = ar_accession_version[0]
			version = ar_accession_version[1]
			nih_gi = item['nih_gi']
			taxon_id = item['taxon_id'] # -> taxon.id
			f.write('\t'.join(map(str, [id, accession, version, nih_gi, taxon_id])) + '\n')
	

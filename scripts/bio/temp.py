
import json
import os

from modules.ncbi import Nucleotide
from modules.ncbi import Taxonomy
from modules.postgresql import Psycopg2

from config import db as db_configs


def temp():

	sr_fpath = '/media/sl/T7-Shield/PaRTISeq/blast_db_genomes/all-export.tsv'
	lr_fpath ='/media/sl/ExtremePro/Parti-Seq/db_nanopore_exports/all-accession-ids.tsv'

	sr_accessions = set()
	lr_accessions = set()

	missing_from_sr = []
	missing_from_lr = []

	with open(sr_fpath) as f:
		f.readline()
		for line in f:
			accession = line.split(' ')[0][1:]
			sr_accessions.add(accession)

	with open(lr_fpath) as f:
		f.readline()
		for line in f:
			accession = line.rstrip('\n')
			lr_accessions.add(accession)
	
	with open('missing_from_lr.tsv', 'w') as f:
		for accession in sr_accessions:
			if accession not in lr_accessions:
				missing_from_lr.append(accession)
				f.write(f'{accession}\n')

	with open('missing_from_sr.tsv', 'w') as f:
		for accession in lr_accessions:
			if accession not in sr_accessions:
				missing_from_sr.append(accession)
				f.write(f'{accession}\n')

def check_for_homonyms(
	fpath = '/media/sl/ExtremePro/Parti-Seq/system-controls/db-data-staging/taxon-name.tsv'
):
	taxon_names = set()
	with open(fpath) as f:
		for line in f:
			ar_line = line.rstrip('\n').split('\t')
			taxon_name = ar_line[0]
			if taxon_name in taxon_names:
				print(f'Homonym: {taxon_name}')
			else:
				taxon_names.add(taxon_name)
		
def gen_nih_genbank_accession_data():

	# Columns
	## accession
	## version
	## gi
	## taxon.nih_id -> taxon.id

	db_staging_fpath = '/media/sl/ExtremePro/Parti-Seq/system-controls/db-data-staging/nih_genbank_accession.tsv'
	f_out = open(db_staging_fpath, 'w')
	unique_accession_versions = set()

################################################################################
# Long Reads DB

	def handle_lr_file(input_dpath, input_fname):
		input_fpath = os.path.join(input_dpath, input_fname)
		with open(input_fpath) as f:
			
			accession_version_col_names = ['accessionversion']
			accession_version_idx: int
			
			gi_col_names = ['gi']
			gi_idx: int

			nih_taxon_id_col_names = ['taxid']
			nih_taxon_id_idx: int
			
			header_line = f.readline()
			ar_header_line = header_line.rstrip('\n').split('\t')
			for idx, col_name in enumerate(ar_header_line):
				if col_name in accession_version_col_names:
					accession_version_idx = idx
				elif col_name in gi_col_names:
					gi_idx = idx
				elif col_name in nih_taxon_id_col_names:
					nih_taxon_id_idx = idx

			for line in f:
				ar_line = line.rstrip('\n').split('\t')
				accession_version = ar_line[accession_version_idx]
				accession = accession_version.split('.')[0]
				version = accession_version.split('.')[1]
				gi = ar_line[gi_idx]
				nih_taxon_id = ar_line[nih_taxon_id_idx]
				if accession_version not in unique_accession_versions:
					unique_accession_versions.add(accession_version)
					out = '\t'.join([accession, version, gi, nih_taxon_id])
					f_out.write(f'{out}\n')

	input_dpath = '/media/sl/ExtremePro/Parti-Seq/db_nanopore_exports/'

	input_fnames = [
		'CHM13v2.0-export-accession-ids-esummaries.tsv',
		'GRCh38p14-export-accession-ids-esummaries.tsv',
		'human-export-accession-ids-esummaries.tsv',
		'Step1_ref_231215-export-accession-ids-esummaries.tsv',
		'Step1_ref_231215-export-unresolved-taxons-esummaries.tsv',
		'Step2_ref_240115-export-accession-ids-esummaries.tsv',
		'Step3_ref_231125-export-accession-ids-esummaries.tsv',
	]

	for input_fname in input_fnames:
		handle_lr_file(input_dpath, input_fname)
	
################################################################################
# Short Reads DB

	additional_unique_accessions = set()

	# Skip these deflines
	spike_ins = ['Allobacillus_halotolerans', 'Imtechella_halotolerans']
	global spike_in_skipped_qty
	spike_in_skipped_qty = 0

	def handle_sr_file(input_dpath, input_fname):
		global spike_in_skipped_qty
		input_fpath = os.path.join(input_dpath, input_fname)
		with open(input_fpath) as f:
			for line in f:
				
				# Skip spike ins
				is_spine_in = False
				for spike_in in spike_ins:
					if spike_in in line:
						spike_in_skipped_qty += 1
						is_spine_in = True
				if is_spine_in:
					continue

				# Split line by spaces and get the first element which is the accession
				ar_line = line.split(' ')
				accession_version = ar_line[0][1:]

				# Raise error if there is no '.' that should separate accession and version
				if '.' not in accession_version:
					raise Exception(f"Missing '.' in accession version.\n{accession_version}")
				
				if accession_version not in unique_accession_versions:
					additional_unique_accessions.add(accession_version)

	input_dpath = '/media/sl/T7-Shield/PaRTISeq/db_short_reads_genomes'

	input_fnames = [
		'Step1_ref_230828-export.tsv',
		'Step2_ref_230902-export.tsv',
		'Step3_ref_230901-export.tsv'
	]

	for input_fname in input_fnames:
		handle_sr_file(input_dpath, input_fname)

	print(f'spike_in_skipped_qty: {spike_in_skipped_qty}')
	print(f'additional_unique_accessions: {additional_unique_accessions}')

	# Get esummaries for these additional accessions

	# Now add them
	input_dpath = '/home/sl/dev/proto-v2/data/scripts/bio'
	input_fname = 'accessions_not_in_lr_db-esummaries.tsv'

	handle_lr_file(input_dpath, input_fname)
								
	f_out.close()

def get_taxon_esummaries():

	Nucleotide.export_esummaries(
		accession_ids_fpath = '/home/sl/dev/proto-v2/data/scripts/bio/accession-ids.tsv'
	)

def get_all_taxon_names():

	out_fpath = '/media/sl/ExtremePro/Parti-Seq/system-controls/db-data-staging/taxon_name.tsv'
	f_out = open(out_fpath, 'w')
	
	all_taxon_ids_fpath = '/home/sl/dev/proto-v2/data/scripts/bio/all-taxon-ids.tsv'
	taxon_ids = {}
	'''
		taxon_ids = {
			'NIH_TAXON_ID': 'TAXON_NAME',
			...
		}
	'''
	with open(all_taxon_ids_fpath) as f:
		for line in f:
			taxon_ids[line.rstrip('\n')] = ''

################################################################################
# Long Reads DB

	def handle_lr_file(input_dpath, input_fname):
		input_fpath = os.path.join(input_dpath, input_fname)
		with open(input_fpath) as f:

			nih_taxon_id_col_names = ['taxid']
			nih_taxon_id_idx: int
			
			taxon_name_col_names = ['organism']
			taxon_name_idx: int
			
			header_line = f.readline()
			ar_header_line = header_line.rstrip('\n').split('\t')
			for idx, col_name in enumerate(ar_header_line):
				if col_name in taxon_name_col_names:
					taxon_name_idx = idx
				elif col_name in nih_taxon_id_col_names:
					nih_taxon_id_idx = idx

			for line in f:
				ar_line = line.rstrip('\n').split('\t')
				nih_taxon_id = ar_line[nih_taxon_id_idx]
				taxon_name = ar_line[taxon_name_idx]
				if nih_taxon_id == 'NULL':
					f_out.write(f"{taxon_name}\t{nih_taxon_id}\n")
					continue
				if taxon_ids[nih_taxon_id] == '':
					taxon_ids[nih_taxon_id] = taxon_name
					f_out.write(f"{taxon_name}\t{nih_taxon_id}\n")
				elif taxon_ids[nih_taxon_id] != taxon_name:
					raise Exception(
						f"Discrepancy: \n"
						f" nih_taxon_id: {nih_taxon_id}\n"
						f" taxon_name: {taxon_ids[nih_taxon_id]}\n"
						f" taxon_name: {taxon_name}"
					)
				

	input_dpath = '/media/sl/ExtremePro/Parti-Seq/db_nanopore_exports/'

	input_fnames = [
		'CHM13v2.0-export-accession-ids-esummaries.tsv',
		'GRCh38p14-export-accession-ids-esummaries.tsv',
		'human-export-accession-ids-esummaries.tsv',
		'Step1_ref_231215-export-accession-ids-esummaries.tsv',
		'Step1_ref_231215-export-unresolved-taxons-esummaries.tsv',
		'Step2_ref_240115-export-accession-ids-esummaries.tsv',
		'Step3_ref_231125-export-accession-ids-esummaries.tsv',
	]

	for input_fname in input_fnames:
		handle_lr_file(input_dpath, input_fname)

	# Taxon IDs for which there is a sub clade representing it
	input_dpath = '/home/sl/dev/proto-v2/data/scripts/bio'
	input_fname = 'additional-taxon-names-ids.tsv'
	handle_lr_file(input_dpath, input_fname)

	
	for nih_taxon_id, taxon_name in taxon_ids.items():
		if taxon_name == '':
			print(nih_taxon_id)

def get_taxon_names_for_true_and_representative_taxon_id():

	taxon_names_fpath = '/media/sl/ExtremePro/Parti-Seq/system-controls/db-data-staging/taxon-name.tsv'
	taxon_names = {}

	with open(taxon_names_fpath) as f:
		for line in f:
			ar_line = line.rstrip('\n').split('\t')
			taxon_name = ar_line[0]
			taxon_id = ar_line[1]
			taxon_names[taxon_id] = taxon_name

	taxon_ids_fpath = '/home/sl/dev/proto-v2/data/scripts/bio/taxon-id-represent-id.tsv'
	taxon_representative_names = {}
	unique_taxon_ids = set()
	
	with open(taxon_ids_fpath) as f:
		f.readline()
		for line in f:
			ar_line = line.rstrip('\n').split('\t')
			taxon_id = ar_line[0]
			representative_id = ar_line[1]
			
			unique_taxon_ids.add(taxon_id)
			unique_taxon_ids.add(representative_id)
			
			taxon_name = taxon_names[taxon_id]
			representative_name = taxon_names[representative_id]
			
			key = f'{taxon_id}_{representative_id}'
			value = f'{taxon_name}\t{representative_name}'
			if key not in taxon_representative_names:
				taxon_representative_names[key] = value
			else:
				print(f'duplicate: {key} {value}')
	
	for taxon_id in unique_taxon_ids:
		if taxon_id not in taxon_names:
			print(f'Missing from taxon_names.tsv: {taxon_id}')
	
	for taxon_id in taxon_names:
		if taxon_id not in unique_taxon_ids:
			print(f'Missing from representative names: {taxon_id}')
	
	out_fpath = '/media/sl/ExtremePro/Parti-Seq/system-controls/db-data-staging/taxon-representative.tsv'
	with open(out_fpath, 'w') as f:
		for taxon_ids, taxon_names in taxon_representative_names.items():
			f.write(f'{taxon_names}\n')

	# MANUALLY ADD THE NAMES WITH NULL IDS

def load_local_parti_cular_db(
	db_config_name = 'postgres',
	new_db_name = 'parti_cular_local',
	schema_fpath = '/media/sl/ExtremePro/Parti-Seq/system-controls/db-schema.json',
	tables_dpath = '/media/sl/ExtremePro/Parti-Seq/system-controls/db-tables-import'
):
	'''
	To all database staging table tsv files (taxon, taxon-name, taxon-representative):
		Needed to manually fix names with '[' and ']'
		Needed to manually add 
			Macrococcus caseolyticus
			Macrococcus canis
			Paeniclostridium sordellii
			Prevotella buccae'
			Check the rest of Prevotella genus
			Adeno-associated virus - 2
			Adeno-associated virus-5
			Human papillomavirus 166
			Human papillomavirus 30
			Human papillomavirus 71
			Human papillomavirus 18
			Human papillomavirus 131
			Human papillomavirus type 101
			Human papillomavirus 112
			Human papillomavirus 128
			Human papillomavirus 129
			Human papillomavirus type 101
			Memnoniella echinata
			Mycobacterium grossiae
			Nosema ceranae
			Primate bocaparvovirus 1
	'''
	'''
	CREATE MATERIALIZED VIEW taxon_representative_view AS
	SELECT 
		tn1.name  as name,
		t1.nih_id as nih_taxon_id,
		tn2.name  as representative_name,
		t2.nih_id as representative_nih_taxon_id
	FROM taxon_representative
	LEFT JOIN taxon_name tn1 ON tn1.id = taxon_representative.taxon_name_id
	LEFT JOIN taxon_name tn2 ON tn2.id = taxon_representative.representative_taxon_name_id
	LEFT JOIN taxon      t1  ON t1.id  = tn1.taxon_id
	LEFT JOIN taxon      t2  ON t2.id  = tn2.taxon_id;
	'''
	db_config = db_configs.postgresql.databases[db_config_name]
	
	conn_config = {
		'dbname': db_config.dbname,
		'user': db_config.user,
		'password': db_config.password,
		'host': db_config.host,
		'port': db_config.port
	}

	# Delete and create database
	psycopg2 = Psycopg2(conn_config)
	psycopg2.open_cursor()

	psycopg2.del_db(new_db_name)
	psycopg2.add_db(new_db_name)

	psycopg2.close_connection()


	conn_config = {
		'dbname': new_db_name,
		'user': db_config.user,
		'password': db_config.password,
		'host': db_config.host,
		'port': db_config.port
	}

	# Connect to newly created database
	psycopg2 = Psycopg2(conn_config)
	psycopg2.open_cursor()

	psycopg2.import_tables(schema_fpath, tables_dpath)

	psycopg2.close_connection()

def get_name_and_genome_length_from_parti_cular_flat_file_db(
	fpaths = [
		'/media/sl/T7-Shield/PaRTISeq/db_short_reads/Pathogen_genome_size_230511.txt',
		'/media/sl/ExtremePro/Parti-Seq/db_nanopore/Pathogen_genome_size_230511.txt'
	],
	out_fpath = 'data/scripts/bio/all_flat_db_species_genome_sizes.tsv'
):
	
	genome_sizes = {}
	'''
	{
		'ORGANISM_NAME': GENOME_SIZE
	}
	'''

	special_names = {
		'Imtechella_halotolerans': 'Imtechella halotolerans',
		'Allobacillus_halotolerans': 'Allobacillus halotolerans'
	}

	fpath = fpaths[0]
	with open(fpath) as f:
		for line in f:
			ar_line = line.strip('\n').split('\t')
			accession = ar_line[0]
			genome_size = int(ar_line[1])
			organism_name = ar_line[2]

			if organism_name in special_names:
				organism_name = special_names[organism_name]
			
			if organism_name not in genome_sizes:
				genome_sizes[organism_name] = genome_size
			else:
				genome_sizes[organism_name] += genome_size

	
	those_genome_sizes = {}
	fpath = fpaths[1]
	with open(fpath) as f:
		for line in f:
			ar_line = line.strip('\n').split('\t')
			accession = ar_line[0]
			genome_size = int(ar_line[1])
			organism_name = ar_line[2]

			if organism_name in special_names:
				organism_name = special_names[organism_name]
			
			if organism_name not in those_genome_sizes:
				those_genome_sizes[organism_name] = genome_size
			else:
				those_genome_sizes[organism_name] += genome_size
	
	genome_size_mismatches = {}

	for organism_name in those_genome_sizes:
		if organism_name not in genome_sizes:
			genome_sizes[organism_name] = those_genome_sizes[organism_name]
		elif genome_sizes[organism_name] != those_genome_sizes[organism_name]:
			print(
				f'Organism: {organism_name} size is in consistent:\n'
				f' {genome_sizes[organism_name]} != {those_genome_sizes[organism_name]}'
			)

			genome_size_mismatches[organism_name] = [
				genome_sizes[organism_name], those_genome_sizes[organism_name]
			]
	
	print(json.dumps(genome_size_mismatches, indent=2))

def export_lineages(
	existing_taxons_fpath = '/media/sl/ExtremePro/Parti-Seq/system-controls/db-data-staging/taxon.tsv',
	new_esummaries_fpath = '/home/sl/dev/proto-v2/data/scripts/bio/new-esummaries.json',
	out_fpath = '/home/sl/dev/proto-v2/data/scripts/bio/taxon-lineages.tsv'
):
	
	# Get existing taxon IDs
	taxon_ids = set()
	with open(existing_taxons_fpath) as f:
		for line in f:
			ar_line = line.rstrip('\n').split('\t')
			taxon_id = ar_line[0]
			taxon_ids.add(taxon_id)
	
	# Get any new taxon IDs
	new_esummaries = json.load(open(new_esummaries_fpath))
	new_taxon_id_qty = 0
	for esummary in new_esummaries:
		taxon_id = str(esummary['taxid'])
		if taxon_id not in taxon_ids:
			new_taxon_id_qty += 1
		taxon_ids.add(taxon_id)
	print(f'new_taxon_id_qty: {new_taxon_id_qty}')

	# Export lineages
	Taxonomy.export_lineages(
		taxon_ids = taxon_ids,
		out_fpath = out_fpath
	)

def separate_lineage_table(
	in_fpath = '/home/sl/dev/proto-v2/data/scripts/bio/taxon-lineages.tsv',
	taxon_fpath = '/media/sl/ExtremePro/Parti-Seq/system-controls/db-data-staging/taxon.tsv',
	ancestor_fpath = '/media/sl/ExtremePro/Parti-Seq/system-controls/db-data-staging/taxon-ancestor.tsv',
	rank_fpath = '/media/sl/ExtremePro/Parti-Seq/system-controls/db-data-staging/taxon-rank.tsv',
	name_fpath = '/media/sl/ExtremePro/Parti-Seq/system-controls/db-data-staging/taxon-name.tsv',
	name_type_fpath = '/media/sl/ExtremePro/Parti-Seq/system-controls/db-data-staging/taxon-name-type.tsv'
):
	"""The taxon-lineages.tsv data is generated by modules.ncbi.Taxonomy.export_lineages() from
	the old taxon.tsv database staging data, which mean the taxon IDs queried are comprehensive.
	However, need to make sure that the data from new-esummaries.json is included.

	Args:
			in_fpath (str, optional): _description_. Defaults to '/home/sl/dev/proto-v2/data/scripts/bio/taxon-lineages.tsv'.
			out_fpath (str, optional): _description_. Defaults to '/media/sl/ExtremePro/Parti-Seq/system-controls/db-data-staging/taxon.tsv'.

	Raises:
			Exception: _description_
	"""
	taxons = dict()
	# {
	#		TAXON_ID: ['NAME', 'RANK', 'PARENT_ID']
	# }
	
	# Get taxonomy data from taxon-lineages.tsv
	with open(in_fpath) as f:
		# Skip header row
		f.readline()

		for line in f:
			ar_line = line.rstrip('\n').split('\t')
			taxon_id = ar_line[0]
			name = ar_line[1]
			rank = ar_line[2]
			ancestor_id = ar_line[3]
			taxons[taxon_id] = [name, rank, ancestor_id]
	
	# Write to taxon.tsv
	#	stage:  nih_taxon_id | rank
	# import: nih_taxon_id | rank_id
	with open(taxon_fpath, 'w') as f:
		for taxon_id, values in taxons.items():
			name = values[0]
			rank = values[1]
			line = '\t'.join(map(str, [taxon_id, name, rank]))
			f.write(f'{line}\n')
	
	# Write to taxon-ancestor.tsv
	#	stage:  nih_taxon_id | nih_ancestor_taxon_id
	# import: taxon_id     | taxon_id
	with open(ancestor_fpath, 'w') as f:
		for taxon_id, values in taxons.items():
			ancestor_id = values[2]
			line = '\t'.join(map(str, [taxon_id, ancestor_id]))
			f.write(f'{line}\n')
	
	# Get unique ranks
	ranks = set()
	for taxon_id, values in taxons.items():
		rank = values[1]
		ranks.add(rank)

	# Write to taxon-rank.tsv
	# stage:  rank
	# import: rank
	with open(rank_fpath, 'w') as f:
		for rank in ranks:
			f.write(f'{rank}\n')
	
	NAME_TYPE_STANDARD = 'standard'
	NAME_TYPE_UNKNOWN = 'unknown'

	# Get existing taxon names
	existing_taxons = {}
	# {
	#		'TAXON_ID': {
	#			'NAME': 'NAME_TYPE,
	#			...
	#		}
	# }
	with open(name_fpath) as f:
		for line in f:
			ar_line = line.rstrip('\n').split('\t')
			name = ar_line[0]
			taxon_id = ar_line[1]

			# Determine whether the name is the standard scientific name
			# The standard names were requested from NIH and are in the
			# taxon-lineages.tsv file.
			standard_name = taxons[taxon_id][0]
			name_type: str
			if name == standard_name:
				name_type = NAME_TYPE_STANDARD
			else:
				name_type = NAME_TYPE_UNKNOWN

			# A taxon id can be share with more than 1 taxon name
			if taxon_id in existing_taxons:
				# Add name to set
				existing_taxons[taxon_id][name] = name_type
			else:
				existing_taxons[taxon_id] = {
					name: name_type
				}
	
	# Add any new taxon names
	added_qty = 0
	for taxon_id, values in taxons.items():
		name = values[0]

		# Determine whether the name is the standard scientific name
		# The standard names were requested from NIH and are in the
		# taxon-lineages.tsv file.
		standard_name = taxons[taxon_id][0]
		name_type: str
		if name == standard_name:
			name_type = NAME_TYPE_STANDARD
		else:
			name_type = NAME_TYPE_UNKNOWN

		# A taxon id can be share with more than 1 taxon name
		if taxon_id in existing_taxons:
			# Add name to set
			existing_taxons[taxon_id][name] = name_type
		else:
			existing_taxons[taxon_id] = {
				name: name_type
			}
			added_qty += 1
	print(f'Added {added_qty} taxon names')

	# Write taxon names to taxon-name.tsv
	# stage:  taxon_name | taxon_name_type    | nih_taxon_id
	# import: taxon_name | taxon_name_type_id | taxon_id
	# name:   name       | type_id            | taxon_id
	with open(name_fpath, 'w') as f:
		for taxon_id, names in existing_taxons.items():
			# All taxon names returned by export_lineages() are standard
			for name, name_type in names.items():
				line = '\t'.join([name, name_type, taxon_id])
				f.write(f'{line}\n')
	
	# Assign standard names (the standardized scientific names)

	# Write taxon name types to taxon-name-type.tsv
	# stage:  taxon_name_type | definition
	# import: taxon_name_type | definition
	# name:   name            | definition
	with open(name_type_fpath, 'w') as f:
		f.write(
			f'{NAME_TYPE_STANDARD}\tthe currently accepted name of a taxon.\n'
			f'{NAME_TYPE_UNKNOWN}\tis not a standard name and is either a custom name, basionym, homotypic synonym, etc.\n'
		)
	
# Alias:
#   accession = genbank_accession
#		org = organism
#		rd = research & development
#		st = staging TSVs
def update_staging_data(
	script_out_dpath = 'data/scripts/bio',
	rd_accessions_fpath = '/media/sl/T7-Shield/PaRTISeq/db_short_reads/Pathogen_genome_size_240415.txt',
	st_accessions_fpath = '/media/sl/ExtremePro/Parti-Seq/system-controls/db-data-staging/nih-genbank-accession.tsv',
	st_taxons_fpath = '/media/sl/ExtremePro/Parti-Seq/system-controls/db-data-staging/taxon.tsv',
	st_taxon_names_fpath = '/media/sl/ExtremePro/Parti-Seq/system-controls/db-data-staging/taxon-name.tsv'
):
	
	rd_spikein_names = {'Imtechella_halotolerans', 'Allobacillus_halotolerans'}

################################################################################
# Update Genbank Accessions
#		nih-genbank-accessions.tsv

	new_esums_fpath = os.path.join(script_out_dpath, 'new-esummaries.json')
	esummaries: list
	if os.path.isfile(new_esums_fpath):
		print('Checking stored esummaries...')
		esummaries = json.load(open(new_esums_fpath))
	
	# Get accessions from R&D data
	# db/Pathogen_genome_size_240415.txt
	rd_accessions = {}
	# {
	# 	'ACCESSION.VERSION': ['GENOME_LENGTH', 'ORGANISM_NAME', '', ''],
	#		...
	# }
	with open(rd_accessions_fpath) as f:
		for line in f:
			ar_line = line.rstrip('\n').split('\t')
			accession_version = ar_line[0]
			genome_length = ar_line[1]
			org_name = ar_line[2]
			# Skip spike ins
			if accession_version in rd_spikein_names:
				continue
			rd_accessions[accession_version] = [genome_length, org_name, '', '']
	
	# Get accessions from staging data
	# staging/nih-genbank-accession.tsv
	st_accessions = {}
	# {
	# 	'ACCESSION.VERSION': ['', '', 'GI_ID', 'TAXON_ID'],
	#		...
	# }
	with open(st_accessions_fpath) as f:
		for line in f:
			ar_line = line.rstrip('\n').split('\t')
			accession = ar_line[0]
			version = ar_line[1]
			accession_version = f'{accession}.{version}'
			gi_id = ar_line[2]
			taxon_id = ar_line[3]
			st_accessions[accession_version] = ['', '', gi_id, taxon_id]
	
	# Check for new accessions
	missing_accessions = {}
	missing_accessions_from_esummary = []
	# {
	# 	'ACCESSION.VERSION': ['GENOME_LENGTH', 'ORGANISM_NAME', '', ''],
	#		...
	# }
	for accession_version, items in rd_accessions.items():
		if accession_version not in st_accessions:
			missing_accessions[accession_version] = items
	if len(missing_accessions) > 1:
		print(
			'Warning: staging data is missing {len(missing_accessions)} accessions.'
		)
		# Get the esummaries and their accession versions
		if os.path.isfile(new_esums_fpath):
			print('Checking stored esummaries...')
			esummaries = json.load(open(new_esums_fpath))
			esummary_accession_versions = set()
			for esummary in esummaries:
				esummary_accession_versions.add(esummary['accessionversion'])

			# Check if accession versions are missing
			for rd_accession_version in missing_accessions:
				if rd_accession_version not in esummary_accession_versions:
					missing_accessions_from_esummary.append(rd_accession_version)
			print(f' Missing esummary qty: {len(missing_accessions_from_esummary)}')
		
		else:
			missing_accessions_from_esummary = list(missing_accessions.keys())
		print(esummaries)
		

		if len(missing_accessions_from_esummary) > 0:
			# Get esummaries of missing accessions
			esummaries.extend(Nucleotide.get_esummaries(
				accession_ids = missing_accessions_from_esummary
			))

		# Write new accessions to file
		with open(new_esums_fpath, 'w') as f:
			f.write(json.dumps(esummaries, indent=2))

		# Add new accessions to staging NIH genbank accessions table
		for esummary in esummaries:
			accession_version = esummary['accessionversion']
			gi_id = esummary['gi']
			nih_taxon_id = esummary['taxid']
			st_accessions[accession_version] = ['', '', gi_id, nih_taxon_id]
		
		# Write updated nih-genbank-accession data to file
		with open(st_accessions_fpath, 'w') as f:
			for accession_version, items in st_accessions.items():
				ar_accession_version = accession_version.split('.')
				accession = ar_accession_version[0]
				version = ar_accession_version[1]
				gi_id = items[2]
				nih_taxon_id = items[3]
				f.write('\t'.join(map(str, [
					accession, version, gi_id, nih_taxon_id
				])) + '\n')

################################################################################
# Update NIH taxon IDs
#		taxon.tsv

	# Get existing NIH taxon IDs
	nih_taxon_ids = set()
	with open(st_taxons_fpath) as f:
		for line in f:
			ar_line = line.rstrip('\n').split('\t')
			nih_taxon_id = ar_line[0]
			nih_taxon_ids.add(nih_taxon_id)
	
	# Add any new NIH taxon IDs
	for esummary in esummaries:
		nih_taxon_ids.add(str(esummary['taxid']))
	
	# Write updated taxon ID data to file
	with open(st_taxons_fpath, 'w') as f:
		for nih_taxon_id in nih_taxon_ids:
			f.write(f'{nih_taxon_id}\n')


################################################################################
# Update taxon
#		taxon-name.tsv

	# Get existing taxon names
	taxon_names = dict()
	with open(st_taxon_names_fpath) as f:
		for line in f:
			ar_Line = line.rstrip('\n').split('\t')
			taxon_name = ar_Line[0]
			taxon_id = ar_Line[1]
			taxon_names[taxon_name] = taxon_id

	# Add new taxon names
	for esummary in esummaries:
		organism_name = esummary['organism']
		nih_taxon_id = str(esummary['taxid'])
		if organism_name in taxon_names:
			existing_nih_taxon_id = taxon_names[organism_name]
			if nih_taxon_id != existing_nih_taxon_id:
				raise Exception()
		else:
			taxon_names[organism_name] = nih_taxon_id
			print(f'added {organism_name}')
	
	# Write updated organism data to file
	with open(st_taxon_names_fpath, 'w') as f:
		for taxon_name, nih_taxon_id in taxon_names.items():
			f.write('\t'.join([taxon_name, str(nih_taxon_id)]) + '\n')


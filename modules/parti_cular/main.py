
import datetime
import os

from modules.postgresql import Psycopg2

from config import db


BACFILE_ROW_NAMES_TABLE_FNAME = 'bacfile-row-names.tsv'
# Alias SYSTEM_CONTROL_SAMPLES_TABLE_FNAME
SC_SAMPLES_TABLE_FNAME = 'system-control-samples.tsv'


def convert_bacfile_to_db_files(
	in_fpath: str = '/home/sl/dev/parti-cular/parti_seq/pipelines/short_reads/background/Bacfile_bac_plasmid_M411011A.txt',
	out_dpath: str = '/media/sl/ExtremePro/Parti-Seq/system-controls/bacfiles'
):
	this_log_dpath = os.path.join('logs', 'modules', 'parti-cular')
	if not os.path.isdir(this_log_dpath): os.makedirs(this_log_dpath)
	this_log_fpath = os.path.join(this_log_dpath, 'convert_bacfile_to_db_files.log')
	# Reset this function's log
	with open(this_log_fpath, 'w') as f:
		pass

	psycopg2 = Psycopg2(
		conn_config = db.postgresql.databases.parti_cular_local
	)

	# Get NIH taxon IDs
	# The taxon name can be used as the key because it's assumed that only species
	# are considered and their should be no homonyms at the species rank.
	taxon_names = {}
	'''
		taxon_names = {
			'NAME': ORGANISM_ID,
			...
		}
	'''
	psycopg2.open_cursor()
	psycopg2.execute_qry("SELECT id, name, taxon_id FROM taxon_name;")
	res = psycopg2.cursor.fetchall()
	psycopg2.close_connection()
	for row in res:
		taxon_name_id = row[0]
		taxon_name = row[1]
		taxon_id = row[2]
		
		# Will do some QA/QC as well
		
		# Add taxon to all_taxon_names
		if taxon_name in taxon_names:
			if taxon_names[taxon_name][0] == taxon_name_id:
				print('WARNING:\n'
					f' assert taxon_names[taxon_name][0] == taxon_name_id\n'
					f' taxon_name: {taxon_name}\n'
					f' taxon_names[taxon_name][0]: {taxon_names[taxon_name][0]}\n'
					f' taxon_name_id: {taxon_name_id}'
				)
		else:
			taxon_names[taxon_name] = [taxon_name_id, taxon_id]

	# Get bacfile row names
	bacfile_row_names = {}
	'''
		bacfile_row_names = {
			'ROW_NAME': ID,
			...
		}
	'''
	psycopg2.open_cursor()
	psycopg2.execute_qry("SELECT * FROM bacfile_row_name;")
	res = psycopg2.cursor.fetchall()
	psycopg2.close_connection()
	for row in res:
		id = row[0]
		row_name = row[1]
		bacfile_row_names[row_name] = id

	# Get system control samples from database
	db_samples = {}
	'''
		system_control_samples = {
			ID: [
				'NAME', 0: Name of the sample, may be duplicates
				'DATE' 	1: Date of the sample
			],
			...
		}
	'''
	psycopg2.open_cursor()
	psycopg2.execute_qry("SELECT * FROM system_control_sample;")
	res = psycopg2.cursor.fetchall()
	psycopg2.close_connection()
	for row in res:
		id = row[0]
		name = row[1]
		date = row[2]
		db_samples[id] = [name, date]

	# Get system control samples from bacfile
	# Alias system_control_sample_qty
	samples = {}
	'''
		sc_sample_names = {
			'SAMPLE_NAME': [
				COL_IDX,      0: Column index in bacfile
				DB_SAMPLE_ID, 1: Sample ID in database
				DATE          2: Date of the sample
			],
			...
		}
	'''
	# Get sample col indexes
	with open(in_fpath) as f:
		line = f.readline()
		ar_line = line.rstrip('\n').split('\t')
		for col_idx in range(1, len(ar_line), 3):
			sample_name = ar_line[col_idx]
			# print(sample_name)
			samples[sample_name] = [col_idx]
			# Get DB sample ID
			for db_sample_id, items in db_samples.items():
				db_sample_name = items[0]
				# print(db_sample_name)
				db_sample_date = items[1]
				if db_sample_name == sample_name:
					samples[sample_name].extend([db_sample_id, db_sample_date])
	del db_samples

	# Create meta files and reads files
	for sample_name, items in samples.items():
		
		# Open bacfile
		with open(in_fpath) as f_in:
			# Skip header row
			f_in.readline()

			col_idx = items[0]
			id = items[1]
			date = items[2]
			year = date.year
			month = date.month
			day = date.day
			this_out_dpath = os.path.join(out_dpath, f'{year}-{month}', f'{day}')
			if not os.path.isdir(this_out_dpath):
				os.makedirs(this_out_dpath)

			# Make meta file
			# Read lines until Taxonomy:
			with open(os.path.join(this_out_dpath, f'{id}-meta.tsv'), 'w') as f_out:
				for line in f_in:
					# Read lines until Taxonomy:
					if line.startswith('Taxonomy:'):
						break
					
					ar_line = line.rstrip('\n').split('\t')
					row_name = ar_line[0]
					# Get bacfile row name id
					row_name_id = bacfile_row_names[row_name]
					# QA/QC row name
					assert row_name in bacfile_row_names
					reads_qty = ar_line[col_idx]
					f_out.write('\t'.join(map(str, [row_name_id, reads_qty])) + '\n')

			# Make species reads file
			# Write the pathogen reads
			with open(os.path.join(this_out_dpath, f'{id}.tsv'), 'w') as f_out:
				for line in f_in:
					ar_line = line.rstrip('\n').split('\t')
					taxon_name = ar_line[0]
					reads_qty = ar_line[col_idx]
					taxon_name_id = taxon_names[taxon_name][0]
					taxon_id = taxon_names[taxon_name][1]
					f_out.write('\t'.join(map(str, [taxon_name_id, taxon_id, reads_qty])) + '\n')


def import_bacfiles_to_db(
	in_fpath: str = '/home/sl/dev/parti-cular/parti_seq/pipelines/short_reads/background/Bacfile_bac_plasmid_M411011A.txt',
	out_dpath: str = '/media/sl/ExtremePro/Parti-Seq/system-controls/bacfiles'
):
	"""_summary_

	Args:
			in_fpath (str): Bacfile with system control samples
			out_fpath_prefix (str): Prefix for separate output files
	"""
	this_log_dpath = os.path.join('logs', 'modules', 'parti-cular')
	if not this_log_dpath: os.makedirs(this_log_dpath)
	this_log_fpath = os.path.join(this_log_dpath, 'import_bacfiles_to_db.log')
	# Reset this function's log
	with open(this_log_fpath, 'w') as f:
		pass

	psycopg2 = Psycopg2(
		conn_config = db.postgresql.databases.parti_cular_local
	)

	# Get NIH taxon IDs
	taxon_names = {}
	'''
		taxon_names = {
			'NAME': NIH_TAXON_ID,
			...
		}
	'''
	psycopg2.open_cursor()
	psycopg2.execute_qry("SELECT * FROM taxon_representative_view;")
	res = psycopg2.cursor.fetchall()
	psycopg2.close_connection()
	for row in res:
		taxon_name = row[0]
		nih_taxon_id = row[1]
		
		# Will do some QA/QC as well
		
		# Add taxon to all_taxon_names
		if taxon_name in taxon_names:
			assert taxon_names[taxon_name] == nih_taxon_id
		else:
			taxon_names[taxon_name] = nih_taxon_id
	
	representative_taxon_names = {}
	'''
		representative_taxon_names = {
			'REPRESENTATIVE_NAME': [
				'NIH_TAXON_ID', | 0 | NIH taxon ID of representative taxon name
				'GROUP_NAME',   | 1 | Taxon Group Name from taxon_names
				'NIH_TAXON_ID'  | 2 | NIH taxon ID of group taxon name
			],
			...
		}
	'''
	for row in res:
		taxon_name = row[0]
		nih_taxon_id = row[1]
		repr_taxon_name = row[2]
		repr_nih_taxon_id = row[3]
		representative_taxon_names[repr_taxon_name] = [
			repr_nih_taxon_id,
			taxon_name,
			nih_taxon_id
		]

	# Get bacfile row names
	bacfile_row_names = {}
	'''
		bacfile_row_names = {
			'ROW_NAME': ID,
			...
		}
	'''
	res = psycopg2.execute_qry("SELECT * FROM bacfile_row_name;")
	for row in res:
		id = row[0]
		row_name = row[1]
		bacfile_row_names[row_name] = id

	# Get system control samples from bacfile
	# Alias system_control_sample_qty
	samples = {}
	'''
		sc_sample_names = {
			'SAMPLE_NAME': {
				'col_idx': COL_IDX,
				'id': ID
			},
			...
		}
	'''
	# Get sample col indexes
	with open(in_fpath) as f:
		line = f.readline()
		ar_line = line.rstrip('\n').split('\t')
		for col_idx in range(1, len(ar_line), 3):
			sample_name = ar_line[col_idx]
			samples[sample_name] = {
				'col_idx': col_idx
			}

	print(f'Found {len(samples)} system control samples.')

	for sample_name in samples:
		res = psycopg2.execute_qry(
			qry = str(
				"INSERT INTO system_control_sample ("
				"	name,"
				"	date,"
				"	method,"
				"	host,"
				"	spike_cell_qty,"
				"	lab_id,"
				"	sample_source_id,"
				"	reagent_fabrication_location_id,"
				"	reagent_lot_id"
				") VALUES (%s, %s, %s, %s, %s, %s, %s, %s, %s)"
			),
			values = (
				sample_name,
				datetime.date(2024, 5, 4),
				'PWF',
				'H',
				'10^4',
				0,
				0,
				0,
				0,
			),
			returning = 'id'
		)
		id = res[0][0]
		samples[sample_name]['id'] = id


	with open(in_fpath) as f_in:
		# Skip header row
		f_in.readline()

		# Create meta files and reads files
		for sample_name in samples:
			col_idx = samples[sample_name]['col_idx']
			id = samples[sample_name]['id']

			# Read lines until Taxonomy:
			with open(os.path.join(out_dpath, f'{id}-meta.tsv'), 'w') as f_out:
				for line in f_in:
					# Read lines until Taxonomy:
					if line.startswith('Taxonomy:'):
						break
				
					ar_line = line.rstrip('\n').split('\t')
					row_name = ar_line[0]
					reads_qty = ar_line[col_idx]
					f_out.write('\t'.join([row_name, str(reads_qty)]) + '\n')

			
			# Write the pathogen reads
			with open(os.path.join(out_dpath, f'{id}.tsv'), 'w') as f_out:
				for line in f_in:
					ar_line = line.rstrip('\n').split('\t')
					species_name = ar_line[0]
					reads_qty = ar_line[col_idx]
					nih_taxon_id = taxon_names[species_name]
					f_out.write('\t'.join(map(str, [nih_taxon_id, reads_qty])) + '\n')
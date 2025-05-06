
import json
import os


def gen_db_table_files(
	config_fpath: str = '/media/sl/ExtremePro/Parti-Seq/system-controls/db-data-staging.json',
	input_dpath: str = '/media/sl/ExtremePro/Parti-Seq/system-controls/db-data-staging',
	output_dpath: str = '/media/sl/ExtremePro/Parti-Seq/system-controls/db-tables-import'
):
	
	# Make output dir if not already
	if not os.path.exists(output_dpath):
		os.makedirs(output_dpath)
	
	config = json.load(open(config_fpath))

	for table in config:

		table_name = table['name']
		fname =f'{table_name.replace("_", "-")}.tsv'
		print(f'Converting {fname}...')
		input_fpath = os.path.join(input_dpath, fname)
		out_fpath = os.path.join(output_dpath, fname)

		# Get IDs for table relationships
		relationships = {}
		'''
			relationships = {
				'fpath': str,
				'map': dict
			}
		'''
		for input_idx, col in enumerate(table['cols']):
			# Get column name
			col_name = list(col.keys())[0]
			if isinstance(col[col_name], list):
				relationships[col_name] = {
					'desc': col[col_name],
					'input_idx': input_idx,
					'fpath': str,
					'map': dict()
				}

				# Get table relationship
				relationship = col[col_name]

				# Get foreign table name and validate
				foreign_table_name = relationship[0].split('.')[0]
				if foreign_table_name != relationship[1].split('.')[0]:
					raise Exception(f"Error: foreign table name {foreign_table_name} is inconsistent.")
				relationships[col_name]['fpath'] = os.path.join(
					output_dpath,
					f'{foreign_table_name.replace("_", "-")}.tsv'
				)

				# Get column name and its index (it's offset by 1 because of ID column)
				# Example: 'taxon_name.name'
				this_col_name = relationship[0].split('.')[1]
				# Default col idx is 0 for table ID column
				this_col_idx = 0
				# Example: 'taxon_name.id'
				that_col_name = relationship[1].split('.')[1]
				# Default col idx is 0 for table ID column
				that_col_idx = 0
				for config_table in config:
					if config_table['name'] == foreign_table_name:
						for idx, config_col in enumerate(config_table['cols'], 1):
							config_col_name = list(config_col.keys())[0]
							if config_col_name == this_col_name:
								this_col_idx = idx
							elif config_col_name == that_col_name:
								that_col_idx = idx
					

				# Map foreign table values to this table relationship
				with open(relationships[col_name]['fpath']) as f:
					for line in f:
						ar_line = line.rstrip('\n').split('\t')
						this_value = ar_line[this_col_idx]
						that_value = ar_line[that_col_idx]
						relationships[col_name]['map'][this_value] = that_value



		with open(input_fpath) as f_in, \
				 open(out_fpath, 'w') as f_out:
			
			next_id = 0

			for line in f_in:
				ar_line = line.rstrip('\n').split('\t')

				# Replace related values with IDs
				for relationship, items in relationships.items():
					input_idx = items['input_idx']
					swap_value = ar_line[input_idx]
					
					try:
						if swap_value not in {'NULL', None, 'None'}:
							new_value = items['map'][swap_value]
						# Change null values to blank ''
						else:
							new_value = ''
					except Exception as e:
						print(
							"Error:\n"
							f" relationship: {items['desc']}\n"
							f" input_idx: {input_idx}\n"
							f" swap_value: {swap_value}"
						)
						raise
					
					ar_line[input_idx] = new_value

				out_line = [str(next_id)]
				out_line.extend(ar_line)
				out_line = '\t'.join(out_line)
				f_out.write(f'{out_line}\n')
				next_id += 1
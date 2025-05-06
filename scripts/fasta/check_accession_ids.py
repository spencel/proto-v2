
import json
import os

from modules import file_sys


def check_for_duplicates(
	fpath: str
):
	
	accession_ids = set()
	duplicates_qty = 0

	with open(fpath) as f:
		# Skip first line
		next(f)
		for i, line in enumerate(f):
			accession_id = line.split('\t')[2]
			if accession_id in accession_ids:
				print(f'Duplicate at line {i}: {accession_id}')
				duplicates_qty += 1
			else:
				accession_ids.add(accession_id)
	
	print(f'Accession IDs:')
	for i, accession_id in enumerate(accession_ids):
		print(f' {accession_id}')
		if i > 2:
			print(' ...')
			break
	print(f'Accession ID qty: {len(accession_ids)}')
	print(f'Duplicate accession ID qty: {duplicates_qty}')


def compare_lists(
	fpaths: list[str],
	accession_ids_col_idx: list[int],
	is_header_row: list[bool],
	results_fpath: str|None = None
):
	results: dict = dict()

	for this_i in range(0, len(fpaths) - 1):
		results[this_i] = dict()
		this_fpath = fpaths[this_i]
		results[this_i]['fpath'] = this_fpath
		this_col_idx = accession_ids_col_idx[this_i]
		this_is_header = is_header_row[this_i]
		results[this_i]['self_duplicates'] = []
		this_acc_list = set()
		with open(this_fpath) as f:
			# Skip header row
			if this_is_header:
				next(f)
			for line in f:
				acc_id = line.strip('\n').split('\t')[this_col_idx]
				if acc_id in this_acc_list:
					results[this_i]['self_duplicates'].append(acc_id)
				else:
					this_acc_list.add(acc_id)
		results[this_i]['qty'] = len(this_acc_list)

		that_i = this_i + 1
		results[that_i] = dict()
		that_fpath = fpaths[that_i]
		results[that_i]['fpath'] = that_fpath
		that_col_idx = accession_ids_col_idx[that_i]
		that_is_header = is_header_row[that_i]
		results[that_i]['self_duplicates'] = []
		that_acc_list = set()
		with open(that_fpath) as f:
			# Skip header row
			if that_is_header:
				next(f)
			for line in f:
				acc_id = line.strip('\n').split('\t')[that_col_idx]
				if acc_id in that_acc_list:
					results[that_i]['self_duplicates'].append(acc_id)
				else:
					that_acc_list.add(acc_id)
		results[that_i]['qty'] = len(that_acc_list)

		print(f'Checking for missing for fpath: {this_i}')
		results[this_i]['missing'] = []
		for that_acc_id in that_acc_list:
			if that_acc_id not in this_acc_list:
				results[this_i]['missing'].append(that_acc_id)

		print(f'Checking for missing for fpath: {that_i}')
		results[that_i]['missing'] = []
		for this_acc_id in this_acc_list:
			if this_acc_id not in that_acc_list:
				results[that_i]['missing'].append(this_acc_id)

	if not results_fpath:
		results_fpath = os.path.join(
			file_sys.File.get_dpath(fpaths[0]),
			'compare-lists-results.json'
		)
	print(results)
	json.dump(results, open(results_fpath, 'w'), indent=2)

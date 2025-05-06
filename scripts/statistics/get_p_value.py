
import json
import os

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.stats as stats
import seaborn as sns

from modules.math import Combinatorics


DATA_DPATH = 'data/scripts/statistics/get_p_value/'


def get_p_value():

################################################################################
# Setting up

	non_perc_dpath = f'{DATA_DPATH}/non-perc'
	perc_dpath = f'{DATA_DPATH}/perc'
	corr_dpath = f'{DATA_DPATH}/correlation_matrix'
	os.makedirs(non_perc_dpath, exist_ok=True)
	os.makedirs(perc_dpath, exist_ok=True)
	os.makedirs(corr_dpath, exist_ok=True)

	results = {
		'reads': {
			'shapiro_wilk_test': {},
			'correlation_matrix': {

			}
		},
		'perc': {
			'shapiro_wilk_test': {}
		}
	}

# END Setting up
################################################################################
# Preparing data - making dataframes

	zscore_thresholds = [3.5, 3.0, 2.5, 2.0, 1.5, 1.0, 0.5]
	dfs = {
		'reads': {
			'low_halo': {
				'data': None,
				'mod_zscore': {
					'all': None
				},
				'zscore': { # Sample Z-scores, not population
					'all': None 
				}
			},
			'no_low_halo': {
				'data': None,
				'mod_zscore': {
					'all': None
				},
				'zscore': {
					'all': None
				}
			}
		},
		'perc': {
			'low_halo': {
				'data': None,
				'mod_zscore': {
					'all': None
				},
				'zscore': {
					'all': None
				}
			},
			'no_low_halo': {
				'data': None,
				'mod_zscore': {
					'all': None
				},
				'zscore': {
					'all': None
				}
			}
		},
		'correlation_matrix': pd.DataFrame(
			columns = [
				'Col1',
				'Col2',
				'reads or perc',
				'low or no halo',
				'zscore or mod zscore',
				'threshold',
				'p correlation'
			]
		)
	}

	dfs['reads']['low_halo']['data'] = pd.read_csv(
		f'{DATA_DPATH}/data.tsv',
		sep = '\t',
		header = 0
	)

	# Low halo data
	dfs['reads']['low_halo']['data']['Total halo reads'] = dfs['reads']['low_halo']['data']['A halo'] + dfs['reads']['low_halo']['data']['I halo']
	dfs['reads']['low_halo']['data']['Total nonhuman reads'] = dfs['reads']['low_halo']['data']['Classified reads'] + dfs['reads']['low_halo']['data']['Unclassified reads']

	# No low halo data
	# Remove outliers: total halo reads < 1000 
	df = dfs['reads']['low_halo']['data']
	drop_indexes = set()
	for i in range(0, len(df.columns)):
		col_name = df.columns[i]
		# Skip non z-score columns
		if col_name != 'Total halo reads':
			continue
		col_data = df.iloc[:, i]
		for j, reads in enumerate(col_data):
			if np.abs(reads) < 1_000:
				drop_indexes.add(j)
	dfs['reads']['no_low_halo']['data'] = df.drop(drop_indexes)

	# Percent low halo data
	df = dfs['reads']['low_halo']['data']
	df_perc = pd.DataFrame()
	for col_name in df:
		# Exclude raw reads and QC reads
		if col_name in ('Raw reads', 'QC reads'):
			continue
		df_perc[col_name] = df[col_name] / df['QC reads']
	dfs['perc']['low_halo']['data'] = df_perc

	# Percent no low halo data
	df = dfs['reads']['no_low_halo']['data']
	df_perc = pd.DataFrame()
	for col_name in df:
		# Exclude raw reads and QC reads
		if col_name in ('Raw reads', 'QC reads'):
			continue
		df_perc[col_name] = df[col_name] / df['QC reads']
	dfs['perc']['no_low_halo']['data'] = df_perc
	del df_perc

	# Get all modified z-scores
	df = dfs['reads']['low_halo']['data']
	df_mod_zscore = pd.DataFrame()
	for col_name in df:
		col_data = df[col_name]
		median = np.median(col_data)
		# Median absolute deviation (MAD)
		mad = np.median(np.abs(col_data - median))
		if mad == 0:
			df_mod_zscore[col_name] = np.zeros_like(col_data)
		else:
			df_mod_zscore[col_name] = 0.6745 * (col_data - median) / mad
	dfs['reads']['low_halo']['mod_zscore']['all'] = df_mod_zscore
	
	drop_indexes = set()
	for threshold in zscore_thresholds:
		for col_name in df_mod_zscore.columns:
			col_data = df_mod_zscore[col_name]
			for index in df_mod_zscore.index:
				ele = col_data.loc[index]
				if np.abs(ele) > threshold:
					drop_indexes.add(index)
		dfs['reads']['low_halo']['mod_zscore'][threshold] = df_mod_zscore.drop(drop_indexes)

	df = dfs['reads']['no_low_halo']['data']
	df_mod_zscore = pd.DataFrame()
	for col_name in df:
		col_data = df[col_name]
		median = np.median(col_data)
		# Median absolute deviation (MAD)
		mad = np.median(np.abs(col_data - median))
		if mad == 0:
			df_mod_zscore[col_name] = np.zeros_like(col_data)
		else:
			df_mod_zscore[col_name] = 0.6745 * (col_data - median) / mad
	dfs['reads']['no_low_halo']['mod_zscore']['all'] = df_mod_zscore

	drop_indexes = set()
	for threshold in zscore_thresholds:
		for col_name in df_mod_zscore.columns:
			col_data = df_mod_zscore[col_name]
			for index in df_mod_zscore.index:
				ele = col_data.loc[index]
				if np.abs(ele) > threshold:
					drop_indexes.add(index)
		dfs['reads']['no_low_halo']['mod_zscore'][threshold] = df_mod_zscore.drop(drop_indexes)

	df = dfs['perc']['low_halo']['data']
	df_mod_zscore = pd.DataFrame()
	for col_name in df:
		col_data = df[col_name]
		median = np.median(col_data)
		# Median absolute deviation (MAD)
		mad = np.median(np.abs(col_data - median))
		if mad == 0:
			df_mod_zscore[col_name] = np.zeros_like(col_data)
		else:
			df_mod_zscore[col_name] = 0.6745 * (col_data - median) / mad
	dfs['perc']['low_halo']['mod_zscore']['all'] = df_mod_zscore

	drop_indexes = set()
	for threshold in zscore_thresholds:
		for col_name in df_mod_zscore.columns:
			col_data = df_mod_zscore[col_name]
			for index in df_mod_zscore.index:
				ele = col_data.loc[index]
				if np.abs(ele) > threshold:
					drop_indexes.add(index)
		dfs['perc']['low_halo']['mod_zscore'][threshold] = df_mod_zscore.drop(drop_indexes)

	df = dfs['perc']['no_low_halo']['data']
	df_mod_zscore = pd.DataFrame()
	for col_name in df:
		col_data = df[col_name]
		median = np.median(col_data)
		# Median absolute deviation (MAD)
		mad = np.median(np.abs(col_data - median))
		if mad == 0:
			df_mod_zscore[col_name] = np.zeros_like(col_data)
		else:
			df_mod_zscore[col_name] = 0.6745 * (col_data - median) / mad
	dfs['perc']['no_low_halo']['mod_zscore']['all'] = df_mod_zscore

	drop_indexes = set()
	for threshold in zscore_thresholds:
		for col_name in df_mod_zscore.columns:
			col_data = df_mod_zscore[col_name]
			for index in df_mod_zscore.index:
				ele = col_data.loc[index]
				if np.abs(ele) > threshold:
					drop_indexes.add(index)
		dfs['perc']['no_low_halo']['mod_zscore'][threshold] = df_mod_zscore.drop(drop_indexes)

	del df_mod_zscore

	# Get all sample z-scores
	df = dfs['reads']['low_halo']['data']
	df_zscore = pd.DataFrame()
	for col_name in df:
		col_data = df[col_name]
		df_zscore[col_name] = (col_data - col_data.mean()) / col_data.std(ddof=1)
	dfs['reads']['low_halo']['zscore']['all'] = df_zscore

	drop_indexes = set()
	for threshold in zscore_thresholds:
		for col_name in df_zscore.columns:
			col_data = df_zscore[col_name]
			for index in df_zscore.index:
				ele = col_data.loc[index]
				if np.abs(ele) > threshold:
					drop_indexes.add(index)
		dfs['reads']['low_halo']['zscore'][threshold] = df_zscore.drop(drop_indexes)
	
	df = dfs['reads']['no_low_halo']['data']
	df_zscore = pd.DataFrame()
	for col_name in df:
		col_data = df[col_name]
		df_zscore[col_name] = (col_data - col_data.mean()) / col_data.std(ddof=1)
	dfs['reads']['no_low_halo']['zscore']['all'] = df_zscore

	drop_indexes = set()
	for threshold in zscore_thresholds:
		for col_name in df_zscore.columns:
			col_data = df_zscore[col_name]
			for index in df_zscore.index:
				ele = col_data.loc[index]
				if np.abs(ele) > threshold:
					drop_indexes.add(index)
		dfs['reads']['no_low_halo']['zscore'][threshold] = df_zscore.drop(drop_indexes)
	
	df = dfs['perc']['low_halo']['data']
	df_zscore = pd.DataFrame()
	for col_name in df:
		col_data = df[col_name]
		df_zscore[col_name] = (col_data - col_data.mean()) / col_data.std(ddof=1)
	dfs['perc']['low_halo']['zscore']['all'] = df_zscore

	drop_indexes = set()
	for threshold in zscore_thresholds:
		for col_name in df_zscore.columns:
			col_data = df_zscore[col_name]
			for index in df_zscore.index:
				ele = col_data.loc[index]
				if np.abs(ele) > threshold:
					drop_indexes.add(index)
		dfs['perc']['low_halo']['zscore'][threshold] = df_zscore.drop(drop_indexes)
	
	df = dfs['perc']['no_low_halo']['data']
	df_zscore = pd.DataFrame()
	for col_name in df:
		col_data = df[col_name]
		df_zscore[col_name] = (col_data - col_data.mean()) / col_data.std(ddof=1)
	dfs['perc']['no_low_halo']['zscore']['all'] = df_zscore

	drop_indexes = set()
	for threshold in zscore_thresholds:
		for col_name in df_zscore.columns:
			col_data = df_zscore[col_name]
			for index in df_zscore.index:
				ele = col_data.loc[index]
				if np.abs(ele) > threshold:
					drop_indexes.add(index)
		dfs['perc']['no_low_halo']['zscore'][threshold] = df_zscore.drop(drop_indexes)

	del df_zscore

	# Save dataframes for QC
	dfs['reads']['low_halo']['data'].to_csv(
		f'{DATA_DPATH}/reads-low-halo.tsv',
		sep = '\t', index = False, header = True
	)
	dfs['reads']['low_halo']['mod_zscore']['all'].to_csv(
		f'{DATA_DPATH}/reads-low-halo-mod-zscore-all.tsv',
		sep = '\t', index = False, header = True
	)
	for threshold in zscore_thresholds:
		dfs['reads']['low_halo']['mod_zscore'][threshold].to_csv(
			f'{DATA_DPATH}/reads-low-halo-mod-zscore-{threshold}.tsv',
			sep = '\t', index = False, header = True
		)
	dfs['reads']['low_halo']['zscore']['all'].to_csv(
		f'{DATA_DPATH}/reads-low-halo-zscore-all.tsv',
		sep = '\t', index = False, header = True
	)
	for threshold in zscore_thresholds:
		dfs['reads']['low_halo']['zscore'][threshold].to_csv(
			f'{DATA_DPATH}/reads-low-halo-zscore-{threshold}.tsv',
			sep = '\t', index = False, header = True
		)
	dfs['reads']['no_low_halo']['data'].to_csv(
		f'{DATA_DPATH}/reads-no-low-halo.tsv',
		sep = '\t', index = False, header = True
	)
	dfs['reads']['no_low_halo']['mod_zscore']['all'].to_csv(
		f'{DATA_DPATH}/reads-no-low-halo-mod-zscore-all.tsv',
		sep = '\t', index = False, header = True
	)
	for threshold in zscore_thresholds:
		dfs['reads']['no_low_halo']['mod_zscore'][threshold].to_csv(
			f'{DATA_DPATH}/reads-no-low-halo-mod-zscore-{threshold}.tsv',
			sep = '\t', index = False, header = True
		)
	dfs['reads']['no_low_halo']['zscore']['all'].to_csv(
		f'{DATA_DPATH}/reads-no-low-halo-zscore-all.tsv',
		sep = '\t', index = False, header = True
	)
	for threshold in zscore_thresholds:
		dfs['reads']['no_low_halo']['zscore'][threshold].to_csv(
			f'{DATA_DPATH}/reads-no-low-halo-zscore-{threshold}.tsv',
			sep = '\t', index = False, header = True
		)
	dfs['perc']['low_halo']['data'].to_csv(
		f'{DATA_DPATH}/perc-low-halo.tsv',
		sep = '\t', index = False, header = True
	)
	dfs['perc']['low_halo']['mod_zscore']['all'].to_csv(
		f'{DATA_DPATH}/perc-low-halo-mod-zscore-all.tsv',
		sep = '\t', index = False, header = True
	)
	for threshold in zscore_thresholds:
		dfs['perc']['low_halo']['mod_zscore'][threshold].to_csv(
			f'{DATA_DPATH}/perc-low-halo-mod-zscore-{threshold}.tsv',
			sep = '\t', index = False, header = True
		)
	dfs['perc']['low_halo']['zscore']['all'].to_csv(
		f'{DATA_DPATH}/perc-low-halo-zscore-all.tsv',
		sep = '\t', index = False, header = True
	)
	for threshold in zscore_thresholds:
		dfs['perc']['low_halo']['zscore'][threshold].to_csv(
			f'{DATA_DPATH}/perc-low-halo-zscore-{threshold}.tsv',
			sep = '\t', index = False, header = True
		)
	dfs['perc']['no_low_halo']['data'].to_csv(
		f'{DATA_DPATH}/perc-no-low-halo.tsv',
		sep = '\t', index = False, header = True
	)
	dfs['perc']['no_low_halo']['mod_zscore']['all'].to_csv(
		f'{DATA_DPATH}/perc-no-low-halo-mod-zscore-all.tsv',
		sep = '\t', index = False, header = True
	)
	for threshold in zscore_thresholds:
		dfs['perc']['no_low_halo']['mod_zscore'][threshold].to_csv(
			f'{DATA_DPATH}/perc-no-low-halo-mod-zscore-{threshold}.tsv',
			sep = '\t', index = False, header = True
		)
	dfs['perc']['no_low_halo']['zscore']['all'].to_csv(
		f'{DATA_DPATH}/perc-no-low-halo-zscore-all.tsv',
		sep = '\t', index = False, header = True
	)
	for threshold in zscore_thresholds:
		dfs['perc']['no_low_halo']['zscore'][threshold].to_csv(
			f'{DATA_DPATH}/perc-no-low-halo-zscore-{threshold}.tsv',
			sep = '\t', index = False, header = True
		)

# END Preparing data - making dataframes
################################################################################	
# Getting correlation matrices
	
	for key, df_mod_zscore in dfs['reads']['low_halo']['mod_zscore'].items():
		df = dfs['reads']['low_halo']['data']
		df = df.loc[df.index.isin(df_mod_zscore.index)]
		df_corr = df.corr()
		df_corr.to_csv(
			f'{corr_dpath}/reads-low-halo-mod-zscore-{key}-correlation-matrix.txt',
			sep = '\t', index = True, header = True
		)
		for i, col_name in enumerate(df_corr.columns):
			col_data = df_corr[col_name]
			for j in range(i + 1, len(col_data)):
				row_name = df_corr.index[j]
				p_corr = col_data.values[j]
				# Skip correlations with a column to itself
				if (col_name == row_name
					  or np.isnan(p_corr)):
					continue
				if key == 'all':
					key = 'none'
				dfs['correlation_matrix'].loc[len(dfs['correlation_matrix'])] = [
					col_name,
					row_name,
					'reads',
					'low halo',
					'mod zscore',
					key,
					p_corr
				]

	for key, df_zscore in dfs['reads']['low_halo']['zscore'].items():
		df = dfs['reads']['low_halo']['data']
		df = df.loc[df.index.isin(df_zscore.index)]
		df_corr = df.corr()
		df_corr.to_csv(
			f'{corr_dpath}/reads-low-halo-zscore-{key}-correlation-matrix.txt',
			sep = '\t', index = True, header = True
		)
		if len(df.index) < 3: # Must have at least 3 rows for corr matrix
			continue
		for i, col_name in enumerate(df_corr.columns):
			col_data = df_corr[col_name]
			for j in range(i + 1, len(col_data)):
				row_name = df_corr.index[j]
				p_corr = col_data.values[j]
				# Skip correlations with a column to itself
				if (col_name == row_name
					  or np.isnan(p_corr)):
					continue
				if key == 'all':
					key = 'none'
				dfs['correlation_matrix'].loc[len(dfs['correlation_matrix'])] = [
					col_name,
					row_name,
					'reads',
					'low halo',
					'zscore',
					key,
					p_corr
				]

	for key, df_mod_zscore in dfs['reads']['no_low_halo']['mod_zscore'].items():
		df = dfs['reads']['low_halo']['data']
		df = df.loc[df.index.isin(df_mod_zscore.index)]
		df_corr = df.corr()
		df_corr.to_csv(
			f'{corr_dpath}/reads-no-low-halo-mod-zscore-{key}-correlation-matrix.txt',
			sep = '\t', index = True, header = True
		)
		if len(df.index) < 3: # Must have at least 3 rows for corr matrix
			continue
		for i, col_name in enumerate(df_corr.columns):
			col_data = df_corr[col_name]
			for j in range(i + 1, len(col_data)):
				row_name = df_corr.index[j]
				p_corr = col_data.values[j]
				# Skip correlations with a column to itself
				if (col_name == row_name
					  or np.isnan(p_corr)):
					continue
				if key == 'all':
					key = 'none'
				dfs['correlation_matrix'].loc[len(dfs['correlation_matrix'])] = [
					col_name,
					row_name,
					'reads',
					'no low halo',
					'mod zscore',
					key,
					p_corr
				]

	for key, df_zscore in dfs['reads']['no_low_halo']['zscore'].items():
		df = dfs['reads']['low_halo']['data']
		df = df.loc[df.index.isin(df_zscore.index)]
		df_corr = df.corr()
		df_corr.to_csv(
			f'{corr_dpath}/reads-no-low-halo-zscore-{key}-correlation-matrix.txt',
			sep = '\t', index = True, header = True
		)
		if len(df.index) < 3: # Must have at least 3 rows for corr matrix
			continue
		for i, col_name in enumerate(df_corr.columns):
			col_data = df_corr[col_name]
			for j in range(i + 1, len(col_data)):
				row_name = df_corr.index[j]
				p_corr = col_data.values[j]
				# Skip correlations with a column to itself
				if (col_name == row_name
					  or np.isnan(p_corr)):
					continue
				if key == 'all':
					key = 'none'
				dfs['correlation_matrix'].loc[len(dfs['correlation_matrix'])] = [
					col_name,
					row_name,
					'reads',
					'no low halo',
					'zscore',
					key,
					p_corr
				]

	for key, df_mod_zscore in dfs['perc']['low_halo']['mod_zscore'].items():
		df = dfs['perc']['low_halo']['data']
		df = df.loc[df.index.isin(df_mod_zscore.index)]
		df_corr = df.corr()
		df_corr.to_csv(
			f'{corr_dpath}/perc-low-halo-mod-zscore-{key}-correlation-matrix.txt',
			sep = '\t', index = True, header = True
		)
		if len(df.index) < 3: # Must have at least 3 rows for corr matrix
			continue
		for i, col_name in enumerate(df_corr.columns):
			col_data = df_corr[col_name]
			for j in range(i + 1, len(col_data)):
				row_name = df_corr.index[j]
				p_corr = col_data.values[j]
				# Skip correlations with a column to itself
				if (col_name == row_name
					  or np.isnan(p_corr)):
					continue
				if key == 'all':
					key = 'none'
				dfs['correlation_matrix'].loc[len(dfs['correlation_matrix'])] = [
					col_name,
					row_name,
					'perc',
					'low halo',
					'mod zscore',
					key,
					p_corr
				]

	for key, df_zscore in dfs['perc']['low_halo']['zscore'].items():
		df = dfs['perc']['low_halo']['data']
		df = df.loc[df.index.isin(df_zscore.index)]
		df_corr = df.corr()
		df_corr.to_csv(
			f'{corr_dpath}/perc-low-halo-zscore-{key}-correlation-matrix.txt',
			sep = '\t', index = True, header = True
		)
		if len(df.index) < 3: # Must have at least 3 rows for corr matrix
			continue
		for i, col_name in enumerate(df_corr.columns):
			col_data = df_corr[col_name]
			for j in range(i + 1, len(col_data)):
				row_name = df_corr.index[j]
				p_corr = col_data.values[j]
				# Skip correlations with a column to itself
				if (col_name == row_name
					  or np.isnan(p_corr)):
					continue
				if key == 'all':
					key = 'none'
				dfs['correlation_matrix'].loc[len(dfs['correlation_matrix'])] = [
					col_name,
					row_name,
					'perc',
					'low halo',
					'zscore',
					key,
					p_corr
				]

	for key, df_mod_zscore in dfs['perc']['no_low_halo']['mod_zscore'].items():
		df = dfs['perc']['low_halo']['data']
		df = df.loc[df.index.isin(df_mod_zscore.index)]
		df_corr = df.corr()
		df_corr.to_csv(
			f'{corr_dpath}/perc-no-low-halo-mod-zscore-{key}-correlation-matrix.txt',
			sep = '\t', index = True, header = True
		)
		if len(df.index) < 3: # Must have at least 3 rows for corr matrix
			continue
		for i, col_name in enumerate(df_corr.columns):
			col_data = df_corr[col_name]
			for j in range(i + 1, len(col_data)):
				row_name = df_corr.index[j]
				p_corr = col_data.values[j]
				# Skip correlations with a column to itself
				if (col_name == row_name
					  or np.isnan(p_corr)):
					continue
				if key == 'all':
					key = 'none'
				dfs['correlation_matrix'].loc[len(dfs['correlation_matrix'])] = [
					col_name,
					row_name,
					'perc',
					'no low halo',
					'mod zscore',
					key,
					p_corr
				]

	for key, df_zscore in dfs['perc']['no_low_halo']['zscore'].items():
		df = dfs['perc']['low_halo']['data']
		df = df.loc[df.index.isin(df_zscore.index)]
		df_corr = df.corr()
		df_corr.to_csv(
			f'{corr_dpath}/perc-no-low-halo-zscore-{key}-correlation-matrix.txt',
			sep = '\t', index = True, header = True
		)
		if len(df.index) < 3: # Must have at least 3 rows for corr matrix
			continue
		for i, col_name in enumerate(df_corr.columns):
			col_data = df_corr[col_name]
			for j in range(i + 1, len(col_data)):
				row_name = df_corr.index[j]
				p_corr = col_data.values[j]
				# Skip correlations with a column to itself
				if (col_name == row_name
					  or np.isnan(p_corr)):
					continue
				if key == 'all':
					key = 'none'
				dfs['correlation_matrix'].loc[len(dfs['correlation_matrix'])] = [
					col_name,
					row_name,
					'perc',
					'no low halo',
					'zscore',
					key,
					p_corr
				]

	dfs['correlation_matrix'].to_csv(
			f'{corr_dpath}/correlation-matrix.txt',
			sep = '\t', index = False, header = True
		)

# END Getting correlation matrices
################################################################################	


# 	col_qty = len(df.columns)
# 	print(f'col_qty: {col_qty}')
# 	print(f'column headers: {df.columns}')

# 	col_paired_combination_qty = Combinatorics.get_binomial_coefficient(
# 		n = len(df.columns),
# 		k = 2
# 	)
# 	print(f'Paired combination count: {col_paired_combination_qty}')
	

# 	for i in range(0, len(df.columns)):
		
# 		col1 = df.iloc[:, i]
# 		col1_name = df.columns[i]

# 		# Q-Q Plot
# 		stats.probplot(col1, dist='norm', plot=plt)
# 		plt.title(f'Q-Q Plot - {col1_name}')
# 		plt.savefig(f'{non_perc_dpath}/q-q-{col1_name}.png')
# 		plt.close('all')

# 		# Shapiro-Wilk Test
# 		stat, p = stats.shapiro(col1)
# 		results['non_percent']['shapiro_wilk_test'][col1_name] = {'p': p}

# 		for j in range(i + 1, len(df.columns)):

# 			col2 = df.iloc[:, j]
# 			col2_name = df.columns[j]
# 			# Histogram
# 			sns.histplot(col1, kde=True, label=col1_name, color='blue')
# 			sns.histplot(col2, kde=True, label=col2_name, color='red')
# 			plt.legend()
# 			plt.title(f'Histogram - {col1_name} - {col2_name}')
# 			plt.savefig(f'{non_perc_dpath}/histogram-{col1_name}-{col2_name}.png')
# 			plt.close('all')

# # END
# ################################################################################
# # Percent of QC reads

# 	# Prepare data
# 	df_perc = pd.DataFrame()
# 	for col_name in df.columns:
# 		# Exclude raw reads and QC reads
# 		if col_name in ('Raw reads', 'QC reads'):
# 			continue
# 		df_perc[col_name] = df[col_name] / df['QC reads']
# 	df_perc.to_csv(
# 		f'{DATA_DPATH}/reads-perc.tsv',
# 		sep = '\t',
# 		index = False,
# 		header = True
# 	)

# 	for i in range(0, len(df_perc.columns)):
		
# 		col1 = df_perc.iloc[:, i]
# 		col1_name = df_perc.columns[i]

# 		# Q-Q Plot
# 		stats.probplot(col1, dist='norm', plot=plt)
# 		plt.title(f'Q-Q Plot - {col1_name}')
# 		plt.savefig(f'{perc_dpath}/q-q-{col1_name}.png')
# 		plt.close('all')

# 		# Shapiro-Wilk Test
# 		stat, p = stats.shapiro(col1)
# 		results['percent']['shapiro_wilk_test'][col1_name] = {'p': p}

# 		for j in range(i + 1, len(df_perc.columns)):

# 			col2 = df_perc.iloc[:, j]
# 			col2_name = df_perc.columns[j]
# 			# Histogram
# 			sns.histplot(col1, kde=True, label=col1_name, color='blue')
# 			sns.histplot(col2, kde=True, label=col2_name, color='red')
# 			plt.legend()
# 			plt.title(f'Histogram - {col1_name} - {col2_name}')
# 			plt.savefig(f'{perc_dpath}/histogram-{col1_name}-{col2_name}.png')
# 			plt.close('all')

# # END Percent of QC reads
# ################################################################################
# # Get Modified Z-score

# 	for i in range(0, len(df.columns)):
# 		col_name = df.columns[i]
# 		data = df.iloc[:, i]
# 		median = np.median(data)
# 		# Median absolute deviation (MAD)
# 		mad = np.median(np.abs(data - median))
# 		if mad == 0:
# 			mod_z_scores = np.zeros_like(data)
# 		else:
# 			mod_z_scores = 0.6745 * (data - median) / mad

# 		df[f'{col_name} Mod Z-score'] = mod_z_scores
# 	df.to_csv(
# 		f'{DATA_DPATH}/reads-mad-z-scores.tsv',
# 		sep = '\t',
# 		index = False,
# 		header = True
# 	)

# 	# Remove outliers
# 	df_no_outliers = df.copy()
# 	drop_indexes = set()
# 	for i in range(0, len(df_no_outliers.columns)):
# 		col_name = df_no_outliers.columns[i]
# 		# Skip non z-score columns
# 		if 'Mod Z-score' not in col_name:
# 			continue
# 		data = df_no_outliers.iloc[:, i]
# 		for j, ele in enumerate(data):
# 			if np.abs(ele) > mod_z_score_threshold:
# 				drop_indexes.add(j)
# 	df_no_outliers = df_no_outliers.drop(drop_indexes)
# 	df_no_outliers.to_csv(
# 		f'{DATA_DPATH}/reads-mad-z-scores-filtered-{mod_z_score_threshold}.tsv',
# 		sep = '\t',
# 		index = False,
# 		header = True
# 	)

# 	df_no_outliers.corr().to_csv(
# 		f'{DATA_DPATH}/reads-mad-z-scores-filtered-{mod_z_score_threshold}-correlation-matrix.txt',
# 		sep = '\t',
# 		index = True,
# 		header = True
# 	)

# # END Get Modified Z-score
# ################################################################################
# # Get Modified Z-score of percents

# 	for i in range(0, len(df_perc.columns)):
# 		col_name = df_perc.columns[i]
# 		data = df_perc.iloc[:, i]
# 		median = np.median(data)
# 		# Median absolute deviation (MAD)
# 		mad = np.median(np.abs(data - median))
# 		if mad == 0:
# 			mod_z_scores = np.zeros_like(data)
# 		else:
# 			mod_z_scores = 0.6745 * (data - median) / mad

# 		df_perc[f'{col_name} Mod Z-score'] = mod_z_scores
# 	df_perc.to_csv(
# 		f'{DATA_DPATH}/perc-mad-z-scores.tsv',
# 		sep = '\t',
# 		index = False,
# 		header = True
# 	)

# 	perc_no_outliers_dfs = {
# 		3.5: None,
# 		3.0: None,
# 		2.5: None,
# 		2.0: None,
# 	}

# 	# Remove outliers
# 	df_perc_no_outliers = df_perc.copy()
# 	drop_indexes = set()
# 	for i in range(0, len(df_perc_no_outliers.columns)):
# 		col_name = df_perc_no_outliers.columns[i]
# 		# Skip non z-score columns
# 		if 'Mod Z-score' not in col_name:
# 			continue
# 		data = df_perc_no_outliers.iloc[:, i]
# 		for j, ele in enumerate(data):
# 			if np.abs(ele) > mod_z_score_threshold:
# 				drop_indexes.add(j)
# 	df_perc_no_outliers = df_perc_no_outliers.drop(drop_indexes)
# 	df_perc_no_outliers.to_csv(
# 		f'{DATA_DPATH}/perc-mad-z-scores-filtered-{mod_z_score_threshold}.tsv',
# 		sep = '\t',
# 		index = False,
# 		header = True
# 	)

# 	df_perc_no_outliers.corr().to_csv(
# 		f'{DATA_DPATH}/perc-mad-z-scores-filtered-{mod_z_score_threshold}-correlation-matrix.txt',
# 		sep = '\t',
# 		index = True,
# 		header = True
# 	)

# # END Get Modified Z-score of percents
# ################################################################################

	json.dump(
		results,
		open(f'{DATA_DPATH}/results.json', 'w'),
		indent = 2
	)
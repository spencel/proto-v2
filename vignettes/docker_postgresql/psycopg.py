import json
import logging
import os

import psycopg2


logging.basicConfig(level=logging.DEBUG)


class Psycopg2:
	
################################################################################
# Properties

	default_col_args = {
		'id': 'SERIAL PRIMARY KEY',
		'name': 'VARCHAR(255) NOT NULL',
		'date': 'DATE'
	}

################################################################################
# Constructor

	def __init__(self, conn_config):
		"""_summary_

		Args:
				conn_config (_type_): {
					"db_name": str,
					"user": str,
					"password": str,
					"host": str,
					"port": int
				}
		"""
		self.conn_config = conn_config
		self.connection = None
		self.cursor = None
		self.is_auto_commit = False

################################################################################
# Connections & Cursors

	def open_connection(self):
		self.connection = psycopg2.connect(
			dbname=self.conn_config['dbname'],
			user=self.conn_config['user'],
			password=self.conn_config['password'],
			host=self.conn_config['host'],
			port=self.conn_config['port']
		)
		if self.is_auto_commit:
			self.connection.autocommit = True
	
	def set_auto_commit(self, is_auto_commit):
		self.is_auto_commit = is_auto_commit
		self.connection.autocommit = is_auto_commit
	
	def open_cursor(self):

		if not self.connection:
			self.open_connection()
		elif self.connection:
			if self.connection.closed:
				self.open_connection()

		if not self.cursor:
			self.cursor = self.connection.cursor()
		elif self.cursor:
			if self.cursor.closed:
				self.cursor = self.connection.cursor()

	def close_cursor(self):
		if self.cursor:
			self.cursor.close()
			logging.debug('Cursor close.')

	def close_connection(self):
		self.close_cursor()
		self.connection.close()
		logging.debug('Connection closed.')
	

################################################################################
# Query Executions

	def exe_qry(self, qry, values = None, returning: str|None=None):
		"""_summary_

		Args:
				qry (_type_): _description_

		Raises:
				Exception: _description_

		Returns:
				_type_: If applicable, returns all records by default.
		"""
		try:
			
			is_close_conn_after = False
			if not self.cursor:
				self.open_cursor()
				is_close_conn_after = True
			elif self.cursor:
				if self.cursor.closed:
					self.open_cursor()
					is_close_conn_after = True

			if returning:
				# Remove ';' if present
				qry = qry.replace(';', '')
				# Append returning statement
				qry = f"{qry}\n  RETURNING {returning};"
			
			if values:
				self.cursor.execute(
					qry,
					values
				)
			else:
				self.cursor.execute(qry)
			
			if not self.is_auto_commit:
				logging.debug('Transaction committed.')
				self.connection.commit()
			elif self.connection.autocommit:
				logging.debug('Transaction auto-committed.')
			if is_close_conn_after:
				self.close_connection()
		
		except psycopg2.OperationalError as e:
				self.connection.rollback()
				self.close_connection()
				raise Exception(e)


################################################################################
# Checking if a Schema element Exists

	def is_db_exists(self, db_name):
		qry = 'SELECT 1 FROM pg_database WHERE datname = %s'
		values = [db_name]
		self.exe_qry(
			qry = qry,
			values = values
		)


################################################################################
# Adding & Deleting Databases

	def add_db(self, db_name):
		self.exe_qry(f"SELECT datname FROM pg_database WHERE datname = '{db_name}'")
		res = self.cursor.fetchone()
		if not res:
			self.set_auto_commit(True)
			self.exe_qry(f'CREATE DATABASE {db_name}')
			self.set_auto_commit(False)
			logging.debug(f'Created database: {db_name}')
		else:
			logging.debug(f'Database already exists: {db_name}')
	
	def del_db(self, db_name):
		self.set_auto_commit(True)
		self.exe_qry(f'DROP DATABASE IF EXISTS {db_name}')
		self.set_auto_commit(False)


################################################################################
# Importing databases

	def _gen_table(self, meta):
		table_name = meta['name']
		qry_str = f"CREATE TABLE IF NOT EXISTS {table_name} ("
		qry_str = f"{qry_str}\n  id {self.default_col_args['id']}"

		for col in meta['cols']:
			# Add comma to previous column statement
			qry_str = f'{qry_str},'

			col_name = list(col.keys())[0]
			args = col[col_name]
			
			if col_name in self.default_col_args and args == "":
				args = self.default_col_args[col_name]
			
			if col_name[-3:] == '_id' and args == "":
				ref_table = col_name[:-3]
				foreign_key = ref_table
				args = f'BIGINT,'
				args = f'{args}\n  CONSTRAINT fk_{foreign_key}'
				args = f'{args}\n  FOREIGN KEY ({col_name})'
				args = f'{args}\n  REFERENCES {ref_table}(id)'
			
			if args.startswith('FOREIGN|'):
				ar_args = args.split('|')
				ref_table = ar_args[-1].split('.')[0]
				foreign_key: str
				if len(ar_args) > 2:
					foreign_key = ar_args[1]
				else:
					foreign_key = ref_table
				ref_col = ar_args[-1].split('.')[1]
				args = f'BIGINT,'
				args = f'{args}\n  CONSTRAINT fk_{foreign_key}'
				args = f'{args}\n  FOREIGN KEY ({col_name})'
				args = f'{args}\n  REFERENCES {ref_table}({ref_col})'

			
			qry_str = f"{qry_str}\n  {col_name} {args}"
		
		# Close query statement
		qry_str = f'{qry_str}\n);'
		logging.debug(f'qry_str: {qry_str}')
		
		self.exe_qry(qry_str)


	def _gen_materialized_view(self, meta):
		# Example:
		# CREATE MATERIALIZED VIEW IF NOT EXISTS genbank_view AS
		# SELECT 
		# 	t.nih_id AS taxon_nih_id,
		# 	tn.name AS taxon_name
		# FROM taxon t
		# JOIN taxon_name tn ON t.id = tn.taxon_id;
		view_name = meta['name']
		# qry_str = f"CREATE MATERIALIZED VIEW IF NOT EXISTS {view_name} AS\n"
		qry_str = f"CREATE MATERIALIZED VIEW {view_name} AS\n"
		qry_str = f"{qry_str}SELECT"

		tables = []

		for col in meta['cols']:
			col_name = list(col)[0]
			relation = col[col_name]
			qry_str = f"{qry_str}\n  {relation} AS {col_name},"
			table_name = relation.split('.')[0]
			if table_name not in tables:
				tables.append(table_name)
		
		# Remove last comma too
		qry_str = f"{qry_str[:-1]}\nFROM {tables[0]}"

		for join in meta['joins']:
			from_relation = list(join.keys())[0]
			from_table = from_relation.split('.')[0]
			to_relation = join[from_relation]
			if '.' not in from_relation:
				from_relation = f'{from_relation}.id'
			if '.' not in to_relation:
				to_relation = f'{to_relation}.id'
			qry_str = f"{qry_str}\nJOIN {from_table} ON {from_relation} = {to_relation}"

		qry_str = f"{qry_str};"
		logging.debug(f'qry_str: {qry_str}')
		
		self.exe_qry(qry_str)


	def import_tables(self, schema_fpath, tables_dpath):
		"""Assumes table name is filename and files have headers containing column
		names. Foreign keys are indicated by column name, eg, 
		reagent_fabrication_facility_id relates to reagent_fabrication_facility
		table.

		Args:
				tables_dpath (_type_): Must only contain files containing table data.
		"""

		# tables = {}
		# '''
		# 	tables = {
		# 		'TABLE_NAME': {
		# 			'fpath': str,
		# 			'cols': [
		# 				{'COL_NAME': 'args'},
		# 				...
		# 			]
		# 		},
		# 		...
		# 	}
		# '''
		
		schema_json = json.load(open(schema_fpath))
		# for table in schema_json:
		# 	fpath = os.path.join(
		# 		tables_dpath,
		# 		table['name'].replace('_', '-') + '.tsv'
		# 	)
		# 	tables[table['name']] = {
		# 		'fpath': fpath,
		# 		'cols': table['cols']
		# 	}

		# Create tables & views
		for meta in schema_json:
			
			if 'type' in meta:

				# Handle materialized view
				if meta['type'] == 'materialized_view':
					# Handle materialized views after data is loaded into database
					continue

			else:

				# Handle table
				self._gen_table(meta)
		
		# Load data into tables
		for meta in schema_json:

			if 'type' in meta:
				# Skip materialized view
				if meta['type'] == 'materialized_view':
					continue
				
			if 'name' in meta:
				table_name = meta['name']
				print(f'importing {table_name}')
				table_fpath = os.path.join(
					tables_dpath,
					table_name.replace('_', '-') + '.tsv'
				)
				with open(table_fpath) as f:
					logging.debug(f'Loading {table_name}...')
					self.cursor.copy_from(f, table_name, sep='\t', null='')
					self.connection.commit()
					
					# Update the next value of the primary key, copy doesn't do this
					self.cursor.execute(f"SELECT MAX(id) FROM {table_name}")
					res = self.cursor.fetchall()
					curr_id = res[0][0]
					next_id = curr_id + 1
					self.cursor.execute(f"ALTER SEQUENCE {table_name}_id_seq RESTART WITH {next_id}")
					self.connection.commit()
		
		# Now create materialized views
		for meta in schema_json:
			if 'type' in meta:
				if meta['type'] == 'materialized_view' \
				and 'fpath' not in meta:
					self._gen_materialized_view(meta)
		
		# Raw SQL materialized views
		for meta in schema_json:
			if 'type' in meta:
				# Skip materialized view
				if meta['type'] != 'materialized_view':
					continue

				if 'fpath' in meta:
					fpath = os.path.join(
						os.path.dirname(schema_fpath),
						meta['fpath']
					)

					lines: str
					with open(fpath) as f:
						lines = ''.join(f.readlines())
					ar_sql = lines.split(';')
					for sql in ar_sql:
						# Skip empty sql
						if sql == '':
							continue
						self.exe_qry(sql)

################################################################################
# Generating database and loading data

	@staticmethod
	def gen_db(
		new_db_name: str, 
		schema_fpath: str,
		tables_dpath: str,
		# db_config_name: str = None,
		db_conn_config: dict = None,
		override_existing_db: bool = False
	):
		"""_summary_

		Args:
				new_db_name (str): _description_
				schema_dpath (str): _description_
				tables_dpath (str): _description_
				db_config_name (str, optional): _description_. Defaults to None.
				db_conn_config (dict, optional): Use this if db_config_name is not
					provided. Must be something other than the new database that will be
					deleted if it already exists. Defaults to None.
		"""
		# if db_config_name:
		# 	db_config = db_configs.postgresql.databases[db_config_name]
		# 	db_conn_config = {
		# 		'dbname': db_config.dbname,
		# 		'user': db_config.user,
		# 		'password': db_config.password,
		# 		'host': db_config.host,
		# 		'port': db_config.port,
		# 	}
		
		# Delete and create database
		psycopg2 = Psycopg2(db_conn_config)
		psycopg2.open_cursor()

		if not override_existing_db:
			if psycopg2.is_db_exists(new_db_name):
				return Exception(f"{new_db_name} already exists.")

		psycopg2.del_db(new_db_name)
		psycopg2.add_db(new_db_name)

		psycopg2.close_connection()

		db_conn_config = {
			'dbname': new_db_name,
			'user': db_conn_config['user'],
			'password': db_conn_config['password'],
			'host': db_conn_config['host'],
			'port': db_conn_config['port']
		}
		
		# Connect to newly created database
		psycopg2 = Psycopg2(db_conn_config)
		psycopg2.open_cursor()

		psycopg2.import_tables(schema_fpath, tables_dpath)

		psycopg2.close_connection()
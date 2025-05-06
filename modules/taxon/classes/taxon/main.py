
import json
import logging
import os
import uuid

from django.db import models


import modules as m


# Configs
bio_config = json.load(open(os.path.join("data", "bio.json"), 'r'))


class Taxon(models.Model):

	class Meta:
		db_table = "taxonomy"

	id = models.UUIDField(
		primary_key=True,
		default=uuid.uuid4,
	)
	name = models.CharField(max_length=255)
	rank = models.CharField(max_length=255)
	common_name = models.CharField(max_length=255)
	parent_id = models.UUIDField()

	@classmethod
	def dump_taxon_tree(cls,
		tree_filepath = None,
		indent_size = 4,
	):

		tree_file = open(tree_filepath, 'w')

		# First select all roots, ie, they have no parent ID
		root_taxons = cls.objects.filter(parent_id=None)

		# Recursive function
		def dump_taxon_tree_recurse(children, level=0):

			# Increment level
			this_level = level + 1

			# Loop through children
			for child in children:
				indent_str = ("|" + " "*(indent_size-1)) * this_level
				tree_file.write(f"{indent_str}{child.name}\n")

				sub_children = cls.objects.filter(parent_id=child.id)
				if sub_children:
					dump_taxon_tree_recurse(sub_children, this_level)

		# Initialize recursion
		for root_taxon in root_taxons:
			# logging.debug(f"{scope_pref}: root_taxon: {root_taxon}")

			tree_file.write(f"{root_taxon.name}\n")

			# Get children
			children = cls.objects.filter(parent_id=root_taxon.id)

			if children:
				dump_taxon_tree_recurse(children)

		tree_file.close()


	@staticmethod
	def export_taxon_tree(
		taxon_relations_fpath,
		taxon_names_fpath,
		export_type: str = 'json',
		indent_size: int = 4,
		out_fpath: str = ''
	):
		# Can't use taxon names to build tree because of homonyms

		# taxon_names_fpath
		# 'TAXON_NAME'  'TYPE'  'TAXON_ID'
		# ...
		#  TYPE: standard, unknown, basionym, etc.
		# We want standard

		# Get taxon names
		# A taxon name can be a homonym
		taxon_names = dict()
		# {
		#   'TAXON_ID': 'TAXON_NAME',
		#   ...
		# }
		with open(taxon_names_fpath) as f:
			for line in f:
				ar_line = line.rstrip('\n').split('\t')
				type = ar_line[1]
				if type == 'standard':
					taxon_name = ar_line[0]
					taxon_id = ar_line[2]
					taxon_names[taxon_id] = taxon_name
				
		# Get roots, they have no parent ID
		root_taxons = []
		# Get taxon data
		taxon_relations = dict()
		# {
		#   'TAXON_ID': 'PARENT_TAXON_ID'
		#   ...
		# }
		with open(taxon_relations_fpath) as f:
			for line in f:
				ar_line = line.rstrip('\n').split('\t')
				taxon_id = ar_line[0]
				parent_id = ar_line[1]
				taxon_relations[taxon_id] = parent_id
				if parent_id == 'None':
					taxon_name = taxon_names[taxon_id]
					root_taxons.append(
						[taxon_id, taxon_name]
					)

		# Sort by 2nd col (the name col)
		root_taxons = sorted(root_taxons, key=lambda x: x[1])
				
		f_tree = open(out_fpath, 'w')

		# Recursive function
		def recurse_export_taxon_tree(children, level=0):
			 
			# Increment level
			this_level = level + 1

			# Loop through children
			for i, child in enumerate(children):
				child_id = child[0]
				child_name = child[1]

				indent_str = (" " + " "*(indent_size-2)) * this_level

				f_tree.write(f'{indent_str} {child_name}\n')

				# Get subchildren
				sub_children = []
				for sub_child_id, parent_id in taxon_relations.items():
					if child_id == parent_id:
						sub_child_name = taxon_names[sub_child_id]
						sub_children.append(
							[sub_child_id, sub_child_name]
						)
				
				# Sort by 2nd col (the name col)
				sub_children = sorted(sub_children, key=lambda x: x[1])
				print(f'{child_name}')
				print(f'{sub_children}')
				

				if len(sub_children) > 0:
					recurse_export_taxon_tree(sub_children, this_level)

		

		for root_taxon in root_taxons:
			root_taxon_id = root_taxon[0]
			root_taxon_name = root_taxon[1]

			f_tree.write(f'{root_taxon_name}\n')

			# Get children
			children = []
			for taxon_id, parent_id in taxon_relations.items():
				if root_taxon_id == parent_id:
					taxon_name = taxon_names[taxon_id]
					children.append(
						[taxon_id, taxon_name]
					)

			# Sort by 2nd col (the name col)
			children = sorted(children, key=lambda x: x[1])
			# print(f'{root_taxon_name}')
			# print(f'{children}')

			if len(children) > 0:
				recurse_export_taxon_tree(children)

		f_tree.close()



	# This function gets the taxonomy ID or creates the lineage
	@classmethod
	def resolve_lineage(cls,
		name = None,
		common_name = None,
		lineage = None
	):
		
		# Flesh out lineage in case it doesn't exist
		prev_taxon = None
		for lineage_name in lineage:

			# logging.debug(f"{scope_pref}: lineage_name: {lineage_name}")

			parent_id = None
			if prev_taxon:
				parent_id = prev_taxon.id
			# logging.debug(f"{scope_pref}: parent_id: {parent_id}")


			this_taxon, created = cls.objects.get_or_create(
				parent_id = parent_id,
				name = lineage_name
			)

			prev_taxon = this_taxon

		parent_id = prev_taxon.id
		# logging.debug(f"{scope_pref}: name: {name}")
		# logging.debug(f"{scope_pref}: common_name: {common_name}")
		# logging.debug(f"{scope_pref}: parent_id: {parent_id}")
		# Add this taxon
		this_taxon, created = cls.objects.get_or_create(
			parent_id = parent_id,
			name = name,
			defaults = {
				"common_name": common_name
			}
		)

		return this_taxon
		
		
	@classmethod
	def load_data_from_file(cls,
		filepath = None
	):

		with open(filepath, 'r') as f:
			# Skip header row
			for line in f:
				headers = line
				break
			# Data
			col_idx = {
				"name": 0,
				"common_name": 1,
				"lineage": 2
			}
			for line in f:
				# Convert line to list
				line = line.rstrip("\n").split("\t")
				name = line[col_idx["name"]]
				common_name = line[col_idx["common_name"]]
				lineage = line[col_idx["lineage"]].split("|")

				cls.resolve_lineage(
					name = name,
					common_name = common_name,
					lineage = lineage
				)

				
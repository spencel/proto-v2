
import json
import os


def main(
	taxon_id_relations_fpath: str = '/media/sl/ExtremePro/Parti-Seq/database-data/db-data-staging/taxon-ancestor.tsv',
		# TAXON_ID	PARENT_TAXON_ID
	taxon_names_fpath: str = '/media/sl/ExtremePro/Parti-Seq/database-data/db-data-staging/taxon-name.tsv',
		# TAXON_NAME	???	TAXON_ID
	taxon_representative_names_fpath: str = '/media/sl/ExtremePro/Parti-Seq/database-data/db-data-staging/taxon-representative.tsv',
		# TAXON_REPRESENTATIVE_NAME	TAXON_NAME
	taxon_ranks_fpath: str = '/media/sl/ExtremePro/Parti-Seq/database-data/db-data-staging/taxon.tsv',
		# TAXON_ID	RANK_NAME
	exclude_taxon_names: list = ['Metazoa'], # Homo sapiens
	exclude_ranks: list = ['no rank']
):
	
	# {
  #   text: "Parent 1",
	# 	icon: "glyphicon glyphicon-stop",
	# 	selectedIcon: "glyphicon glyphicon-stop",
	# 	color: "#000000",
	# 	backColor: "#FFFFFF",
	# 	href: "#node-1",
	# 	selectable: true,
	# 	state: {
	# 		checked: true,
	# 		disabled: true,
	# 		expanded: true,
	# 		selected: true
	# 	},
	# 	tags: ['available'],
  #   nodes: [
  #     {
  #       text: "Child 1",
  #       nodes: [
  #         {
  #           text: "Grandchild 1"
  #         },
	#					...
  #       ]
  #     },
  #     {
  #       text: "Child 2"
  #     },
	#			...
  #   ]
  # },
  # ...

	out_dpath = os.path.join(
		'data', 'scripts', 'bio', 'gen_taxon_tree_json_for_bootstrap_treeview'
	)
	if not os.path.exists(out_dpath):
		os.makedirs(out_dpath)


################################################################################
# Making Taxon ID Tree

	root_taxon_ids = []
	taxon_id_relations = dict()
	# {
	#   'TAXON_ID': 'PARENT_TAXON_ID'
	#   ...
	# }
	with open(taxon_id_relations_fpath) as f:
		for line in f:
			ar_line = line.rstrip('\n').split('\t')
			taxon_id = ar_line[0]
			parent_taxon_id = ar_line[1]
			taxon_id_relations[taxon_id] = parent_taxon_id
			if parent_taxon_id == 'None':
					root_taxon_ids.append(taxon_id)

	id_tree = dict()

	# Recursive function
	def recurse_taxon_id_tree(subtree, child_ids):

		# Loop through children
		for i, child_id in enumerate(child_ids):

			subtree[child_id] = dict()

			# Get subchildren
			subchild_ids = []
			for subchild_id, parent_id in taxon_id_relations.items():
				if child_id == parent_id:
					subchild_ids.append(subchild_id)			

			if len(subchild_ids) > 0:
				recurse_taxon_id_tree(subtree[child_id], subchild_ids)

	for root_taxon_id in root_taxon_ids:

		id_tree[root_taxon_id] = dict()

		# Get children
		child_ids = []
		for child_id, parent_id in taxon_id_relations.items():
			if root_taxon_id == parent_id:
				child_ids.append(child_id)

		if len(child_ids) > 0:
			recurse_taxon_id_tree(id_tree[root_taxon_id], child_ids)
		else:
			id_tree[root_taxon_id] = None
	

	with open(os.path.join(out_dpath, 'id-tree.json'), 'w') as f:
		f.write(json.dumps(id_tree, indent=2))


################################################################################
# Making Taxon Tree & Sorting
	
	taxon_id_x_names = dict()
	# {
	#   'TAXON_ID': [
	#			'TAXON_NAME',
	#			...
	#		]
	#   ...
	# }
	with open(taxon_names_fpath) as f:
		for line in f:
			ar_line = line.rstrip('\n').split('\t')
			# TAXON_NAME	???	TAXON_ID
			taxon_name = ar_line[0]
			taxon_id = ar_line[2]
			if taxon_id not in taxon_id_x_names:
				taxon_id_x_names[taxon_id] = [taxon_name]
			else:
				taxon_id_x_names[taxon_id].append(taxon_name)
	
	with open(os.path.join(out_dpath, 'taxon-id-x-names.json'), 'w') as f:
		f.write(json.dumps(taxon_id_x_names, indent=2))

	taxon_name_x_id = dict()
	for taxon_id, taxon_names in taxon_id_x_names.items():
		for taxon_name in taxon_names:
			taxon_name_x_id[taxon_name] = taxon_id

	# # Add taxon representative nams to the taxon_id_x_names dict
	# with open(taxon_representative_names_fpath) as f:
	# 	for line in f:
	# 		ar_line = line.rstrip('\n').split('\t')
	# 		# TAXON_REPRESENTATIVE_NAME	TAXON_NAME
	# 		taxon_representative_name = ar_line[0]
	# 		taxon_name = ar_line[1]
	# 		taxon_id = taxon_name_x_id[taxon_representative_name]
	# 		if taxon_name not in taxon_id_x_names[taxon_id]:
	# 			taxon_id_x_names[taxon_id].append(taxon_name)
	
	# with open(os.path.join(out_dpath, 'taxon-id-x-names.json'), 'w') as f:
	# 	f.write(json.dumps(taxon_id_x_names, indent=2))

	taxon_ranks = dict()
	# {
	#   'TAXON_ID': 'TAXON_RANK'
	#   ...
	# }
	with open(taxon_ranks_fpath) as f:
		for line in f:
			ar_line = line.rstrip('\n').split('\t')
			taxon_id = ar_line[0]
			taxon_rank = ar_line[1]
			taxon_ranks[taxon_id] = taxon_rank
	
	taxon_tree = []

	def recurse_taxon_tree(sub_id_tree, sub_taxon_tree):
		for taxon_id in sub_id_tree:
			taxon_names = taxon_id_x_names[taxon_id]

			for taxon_name in taxon_names:
				taxon_rank = taxon_ranks[taxon_id]
				sub_taxon_tree.append({
					'taxon_name': taxon_name,
					'taxon_id': taxon_id,
					'taxon_rank': taxon_rank
				})
		
				if sub_id_tree[taxon_id]:
					sub_taxon_tree[-1]['children'] = []
					sub_taxon_tree[-1]['children'] = recurse_taxon_tree(
						sub_id_tree[taxon_id],
						sub_taxon_tree[-1]['children']
					)
		
		sub_taxon_tree = sorted(sub_taxon_tree, key=lambda x: x['taxon_name'])
		return sub_taxon_tree



	for taxon_id in id_tree:
		taxon_names = taxon_id_x_names[taxon_id]

		for taxon_name in taxon_names:
			taxon_rank = taxon_ranks[taxon_id]
			taxon_tree.append({
				'taxon_name': taxon_name,
				'taxon_id': taxon_id,
				'taxon_rank': taxon_rank
			})
		
			if id_tree[taxon_id]:
				taxon_tree[-1]['children'] = []
				taxon_tree[-1]['children'] = recurse_taxon_tree(
					id_tree[taxon_id],
					taxon_tree[-1]['children']
				)

	taxon_tree = sorted(taxon_tree, key=lambda x: x['taxon_name'], reverse=True)
	
	with open(os.path.join(out_dpath, 'taxon-tree.json'), 'w') as f:
		f.write(json.dumps(taxon_tree, indent=2))
	

################################################################################
# Removing Excluded Taxons

	# Only species in Metazoa is Homo sapiens


################################################################################
# Making Taxon Tree for Bootstrap

	bootstrap_tree = []

	def recurse_bootstrap_tree(children, nodes):
		for child in children:
			taxon_rank = f" {child['taxon_rank']}" if child['taxon_rank'] not in exclude_ranks else ''
			link = f"<link href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={child['taxon_id']}\">{child['taxon_id']}</link>"
			text = f"{child['taxon_name']} (ID:{link}){taxon_rank}"
			nodes.append({
				'text': text
			})
			if 'children' in child:
				nodes[-1]['nodes'] = []
				nodes[-1]['nodes'] = recurse_bootstrap_tree(
					child['children'],
					nodes[-1]['nodes']
				)
		return nodes

	for taxon in taxon_tree:
		taxon_rank = f" {taxon['taxon_rank']}" if taxon['taxon_rank'] not in exclude_ranks else ''
		link = f"<link href=\"https://www.ncbi.nlm.nih.gov/Taxonomy/Browser/wwwtax.cgi?id={taxon['taxon_id']}\">{taxon['taxon_id']}</link>"
		text = f"{taxon['taxon_name']} (ID:{link}){taxon_rank}"
		bootstrap_tree.append({
			'text': text
		})
		if 'children' in taxon:
			bootstrap_tree[-1]['nodes'] = []
			bootstrap_tree[-1]['nodes'] = recurse_bootstrap_tree(
				taxon['children'],
				bootstrap_tree[-1]['nodes']
			)
	
	with open(os.path.join(out_dpath, 'bootstrap-tree.json'), 'w') as f:
		f.write(json.dumps(bootstrap_tree, indent=2))
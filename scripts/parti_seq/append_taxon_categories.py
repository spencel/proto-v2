
import json

from modules.fasta import Fasta
from modules.file_sys import File


# Adds new organism names to the Tax_cate.txt file
# Also exports a list of all organisms and their categories in the current genome database
def append_taxon_categories(
	fasta_fpaths: list[str] = [
		'/media/sl/ExtremePro/Parti-Seq/long-reads-db/Step1_ref.fa',
		'/media/sl/ExtremePro/Parti-Seq/long-reads-db/Step2_ref.fa',
		'/media/sl/ExtremePro/Parti-Seq/long-reads-db/Step3_ref.fa'
	],
	nih_genbank_accession_data_fpath: str = '/media/sl/ExtremePro/Parti-Seq/database-data/db-data-staging/nih-genbank-accession.tsv',
	taxon_categories_fpath: str = '/media/sl/ExtremePro/Parti-Seq/long-reads-db/Tax_cate.txt',
	taxon_name_data_fpath: str = '/media/sl/ExtremePro/Parti-Seq/database-data/db-data-staging/taxon-name.tsv',
	taxon_ancestor_data_fpath: str = '/media/sl/ExtremePro/Parti-Seq/database-data/db-data-staging/taxon-ancestor.tsv'
):
	
	accession_ids_fpaths = []
	for fasta_fpath in fasta_fpaths:
		this_fasta = Fasta(
			fpath = fasta_fpath
		)
		out_fpath = File.get_fpath_without_extension(
        fasta_fpath
      ) + '-accession-ids.txt'
		accession_ids_fpaths.append(out_fpath)
		# this_fasta.export_accession_ids(
		# 	out_fpath = out_fpath
		# )
		del this_fasta
	
	# ~-accession-ids.txt
	#		0 accession_id
	accession_ids = set()
	for accession_ids_fpath in accession_ids_fpaths:
		with open(accession_ids_fpath) as f:
			for line in f:
				accession_id = line.rstrip('\n')
				if accession_id not in ['Imtechella_halotolerans', 'Allobacillus_halotolerans']:
					accession_ids.add(accession_id)
				
	# nih-genbank-accession.tsv
	# 	0	accession_id (without version)
	#		1	accession_version
	#		3	taxon_id
	taxon_ids = dict()
	# {
	#		ACCESSION_ID: TAXON_ID,
	#		...
	#	}
	with open(nih_genbank_accession_data_fpath) as f:
		for line in f:
			ar_line = line.rstrip('\n').split('\t')
			accession_id = f'{ar_line[0]}.{ar_line[1]}'
			taxon_id = ar_line[3]
			taxon_ids[accession_id] = taxon_id

	#	taxon-name.tsv
	#		0	taxon_name
	#		2	taxon_id
	taxon_names = dict()
	# {
	#		TAXON_ID: TAXON_NAME,
	#		...
	#	}
	with open(taxon_name_data_fpath) as f:
		for line in f:
			ar_line = line.rstrip('\n').split('\t')
			taxon_name = ar_line[0]
			taxon_id = ar_line[2]
			taxon_names[taxon_id] = taxon_name

	# taxon-ancestors.tsv
	#		0	taxon_id
	#		1 parent_taxon_id
	parent_taxon_ids = dict()
	# {
	#		TAXON_ID: PARENT_TAXON_ID,
	#		...
	#	}
	with open(taxon_ancestor_data_fpath) as f:
		for line in f:
			ar_line = line.rstrip('\n').split('\t')
			taxon_id = ar_line[0]
			parent_taxon_id = ar_line[1]
			parent_taxon_ids[taxon_id] = parent_taxon_id
	# print(json.dumps(parent_taxon_ids, indent=2))
	
	used_taxon_ids = []
	taxon_categories = dict()
	# {
	#		ORGANISM_NAME: CATEGORY_NAME,
	#		...
	#	}
	for accession_id in accession_ids:
		taxon_id = taxon_ids[accession_id]

		if taxon_id in used_taxon_ids:
			continue
		else:
			used_taxon_ids.append(taxon_id)

		taxon_name = taxon_names[taxon_id]
		while True:
			next_taxon_name = taxon_names[taxon_id]
			
			match next_taxon_name.lower():
				case 'bacteria' \
				| 'archaea':
					taxon_categories[taxon_name] = 'Bacteria'
					break
				case 'fungi':
					taxon_categories[taxon_name] = 'Fungi'
					break
				case 'viruses':
					taxon_categories[taxon_name] = 'Virus'
					break
				case 'cellular organisms':
					taxon_categories[taxon_name] = 'parasites'
					break

			parent_taxon_id = parent_taxon_ids[taxon_id]
			taxon_id = parent_taxon_id
	
	with open('current-organisms.tsv', 'w') as f:
		for org_name, category_name in taxon_categories.items():
			f.write('\t'.join([org_name, category_name]) + '\n')
	
	#	Tax_cate.txt
	#		0	taxon_name
	#		1	category (bacteria, fungi, virus, parasites)
	existing_taxon_categories = dict()
	# {
	#		TAXON_NAME: TAXON_CATEGORY,
	#		...
	#	}
	with open(taxon_categories_fpath) as f:
		for line in f:
			ar_line = line.rstrip('\n').split('\t')
			taxon_name = ar_line[0]
			taxon_category = ar_line[1]
			taxon_categories[taxon_name] = taxon_category
	
	for taxon_name, taxon_category in taxon_categories.items():
		if taxon_name not in existing_taxon_categories:
			existing_taxon_categories[taxon_name] = taxon_category
	
	with open(taxon_categories_fpath, 'w') as f:
		for taxon_name, taxon_category in existing_taxon_categories.items():
			f.write('\t'.join([taxon_name, taxon_category]) + '\n')

import json

from modules.fasta import Fasta, DEFLINE_SELECTION
from modules.file_sys import File


def count_genomes(
	fasta_fpaths: list[str] = [
		'/media/sl/ExtremePro/Parti-Seq/long-reads-db/Step1_ref.fa',
		'/media/sl/ExtremePro/Parti-Seq/long-reads-db/Step2_ref.fa',
		'/media/sl/ExtremePro/Parti-Seq/long-reads-db/Step3_ref.fa'
	],
	nih_genbank_accession_data_fpath: str = '/media/sl/ExtremePro/Parti-Seq/database-data/db-data-staging/nih-genbank-accession.tsv',
	taxon_name_data_fpath: str = '/media/sl/ExtremePro/Parti-Seq/database-data/db-data-staging/taxon-name.tsv',
	taxon_categories_fpath: str = '/media/sl/ExtremePro/Parti-Seq/long-reads-db/Tax_cate.txt'
):
	
	accession_ids_fpaths = []
	for fasta_fpath in fasta_fpaths:
		fasta = Fasta(
			fpath = fasta_fpath
		)
		out_fpath = File.get_fpath_without_extension(
			fasta_fpath
		) + '-accession-ids.txt'
		fasta.export_accession_ids(
			out_fpath,
			selections = [DEFLINE_SELECTION.NON_MITOCHONDRIA]
		)
		accession_ids_fpaths.append(out_fpath)

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

	#	Tax_cate.txt
	#		0	taxon_name
	#		1	category (bacteria, fungi, virus, parasites)
	taxon_categories = dict()
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
	
	used_taxon_ids = []
	counts = {
		'bacteria': 0,
		'fungi': 0,
		'virus': 0,
		'parasites': 0
	}
	for accession_id in accession_ids:
		taxon_id = taxon_ids[accession_id]
		taxon_name = taxon_names[taxon_id]
		if taxon_name in taxon_categories \
		and taxon_id not in used_taxon_ids:
			taxon_category = taxon_categories[taxon_name].lower()
			used_taxon_ids.append(taxon_id)
			counts[taxon_category] += 1
		elif taxon_name not in taxon_categories:
			print(taxon_name)
	
	print(json.dumps(counts, indent=2))
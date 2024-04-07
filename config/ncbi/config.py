
import json
import os

from ..util import *


ncbi = {

  # https://www.ncbi.nlm.nih.gov/refseq/
  'refseq': {
    # https://www.ncbi.nlm.nih.gov/books/NBK21091/table/ch18.T.refseq_accession_numbers_and_mole/?report=objectonly/
    'accession_prefixes': {}
  },
  
  # https://www.ncbi.nlm.nih.gov/genbank/
  'genbank': {
    # https://www.ncbi.nlm.nih.gov/genbank/acc_prefix/
    'accession_prefixes': {}
  }
}


# Env
env = json.load(open(
  os.path.join('env', 'ncbi.json'), 'r'
))
for key in env:
  ncbi[key] = env[key]

# RefSeq
with open(
  os.path.join('config', 'ncbi', 'refseq-accession-prefixes.tsv'), 'r'
) as f:
  # Skip first line
  for line in f:
    break
  # Continue on 2nd line
  for line in f:
    ar_line = line.rstrip('\n').split('\t')
    accession_pref = ar_line[0]
    molecule_type = ar_line[1]
    comment = ar_line[2]
    ncbi['refseq']['accession_prefixes'][accession_pref] = {
      'molecule_type': molecule_type,
      'commnent': comment
    }

# Genbank
with open(
  os.path.join('config', 'ncbi', 'genbank-accession-prefixes.tsv'), 'r'
) as f:
  # Skip first line
  for line in f:
    break
  # Continue on 2nd line
  for line in f:
    ar_line = line.rstrip('\n').split('\t')
    accession_pref = ar_line[0]
    molecule_type = ar_line[1]
    comment = ar_line[2]
    ncbi['genbank']['accession_prefixes'][accession_pref] = {
      'molecule_type': molecule_type,
      'commnent': comment
    }

ncbi = create_dotdict(ncbi, is_mutable=False)

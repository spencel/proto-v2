
import commands
import modules as m
import test
import scripts
from modules.test import Class_A


def main_test():
  m.test.test()
  m.test.Class_A.test()
  Class_A.test()
  m.test.Class_A.Class_C.test()
  m.test.Class_A.Class_D.test()
  m.test.Class_B.test()
  m.test.helpers.test()
  m.test.helper_scripts.test()
  m.test.scripts.test()
  
# def test():
#   m.taxon.Taxon.export_taxon_tree(
#     '/media/sl/ExtremePro/Parti-Seq/system-controls/db-data-staging/taxon-ancestor.tsv',
#     '/media/sl/ExtremePro/Parti-Seq/system-controls/db-data-staging/taxon-name.tsv',
#     out_fpath = 'temp-taxon.tree'
#   )

if __name__ == "__main__":
  # Example command:
  #   python -m fire commands/fastq/main.
  pass

import os
import re

import modules as m

from ..fasta import Fasta


class Fastas():
  
  def __init__(self,
    fastas: list[Fasta]|list[str]
  ):
    """_summary_

    Args:
        fastas (list[m.fasta.Fasta] | list[str]): A list of Fasta instances, fasta filepaths, or directories.
    """
    self.fastas = []

    for fasta in fastas:
      if isinstance(fasta, str):
        path = fasta
        if os.path.isfile(path):
          fpath = path
          self.fastas.append(Fasta(fpath=fpath))
        elif os.path.isdir(path):
          dpath = path
          for name in os.listdir(dpath):
            if not re.search('.f(a|asta)$', name):
              continue
            fpath = os.path.join(dpath, name)
            self.fastas.append(Fasta(fpath=fpath))
      else:
        self.fastas.append(fasta)

  
  def export(self,
    seq_hash_type: str|None = m.hash.SHA256
  ) -> list[dict]:
    """_summary_

    Args:
        is_deflines (bool, optional): _description_. Defaults to True.
        seq_hash (str | None, optional): _description_. Defaults to None.

    Returns:
        dict: _description_
    """
    out_dpath = None
    out_metrics = []
    
    for fasta in self.fastas:
      out_metrics.append(
        fasta.export(
          out_dpath = out_dpath,
          seq_hash_type = seq_hash_type
        )
      )
    
    return out_metrics
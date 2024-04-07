
import config as c


class RefSeq():
  """_summary_
  https://www.ncbi.nlm.nih.gov/refseq/
  
  Returns:
      _type_: _description_
  """

  meta = c.ncbi.refseq

  @classmethod
  def get_meta(cls):
    print(cls.meta)
    return cls.meta

import modules as m


class Text():
  
  file_extensions = {
    '\t': '.tsv',
    ',': '.csv'
  }
  
  @classmethod
  def vertical_list_to_table(cls,
    fpath: str,
    lines_per_row: int,
    out_fpath: str|None = None,
    delim: str = '\t'
  ) -> str:
    """_summary_

    Args:
        fpath (str): _description_
        lines_per_row (int): How many lines in the text file that should be in each row.
        out_fpath (str | None, optional): _description_. Defaults to None.
        delim (str, optional): _description_. Defaults to '\t'.

    Returns:
        str: out_fpath
    """
    file_extension = cls.file_extensions[delim]

    if not out_fpath:
      out_fpath = m.File.get_fpath_without_extension(
        fpath
      ) + '-table' + file_extension

    with open(fpath, 'r') as f_in, \
      open(out_fpath, 'w') as f_out:


      out_line = []
      for line in f_in:
        out_line.append(line[:-1])
        break

      for i, line in enumerate(f_in, start=2):
        out_line.append(line[:-1])

        if i % lines_per_row == 0:
          f_out.write(delim.join(out_line) + '\n')
          out_line = []
    
    return out_fpath
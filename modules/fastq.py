
import gzip
import random
import sys

import config
import modules as m


class Fastq():

  def __init__(
    self,
    fpath
  ):
    self.fpath = fpath


  @staticmethod
  def generate_mock_from_file(
    input_fpath,
    output_fpath,
    base_qty = None,
    file_size = None
  ) -> 'Fastq':
    
    f_out = open(output_fpath, 'w')

    base_swap = {
      'A': 'T',
      'T': 'G',
      'G': 'C',
      'C': 'A',
      'N': 'A'
    }

    if base_qty:
      
      base_written_qty = 0
      while base_written_qty < base_qty:

        # Handle gz files
        if input_fpath.endswith(".gz"):
          with gzip.open(input_fpath, "rb") as f_in:
            i = 1
            for line in f_in:
              sys.stdout.write("\r")
              line = line.decode('utf-8')
              # Sequence line is every 4 lines sequence line starts at line 2
              if (i - 2) % 4 == 0:
                i_base = random.randint(0, len(line) - 2)
                old_base = line[i_base]
                new_base = base_swap[old_base]
                line = line[:i_base] + new_base + line[i_base + 1:]
                base_written_qty += len(line) - 1
              f_out.write(line)
              if (base_written_qty) % 1_000_000 == 0:
                sys.stdout.write(f" {base_written_qty / 1_000_000} Mbp")
              i += 1
              if base_written_qty >= base_qty:
                # Write last 2 lines
                f_out.write(next(f_in).decode('utf-8'))
                f_out.write(next(f_in).decode('utf-8'))
                break
        
        else:
          with open(input_fpath, "r") as f_in:
            i = 1
            for line in f_in:
              # Sequence line is every 4 lines sequence line starts at line 2
              if (i - 2) % 4 == 0:
                i_base = random.randint(0, len(line) - 2)
                old_base = line[i_base]
                new_base = base_swap[old_base]
                line = line[:i_base] + new_base + line[i_base + 1:]
                base_written_qty += len(line) - 1
              f_out.write(line)
              if (base_written_qty) % 1_000_000 == 0:
                sys.stdout.write(f" {base_written_qty / 1_000_000} Mbp")
              i += 1
              if base_written_qty >= base_qty:
                # Write last 2 lines
                f_out.write(next(f_in))
                f_out.write(next(f_in))
                break

      print(f" {base_written_qty / 1_000_000} Mbp")

    f_out.close()

    return Fastq(output_fpath)
  

  def generate_mock_from_self(
    self,
    output_fpath,
    base_qty = None,
    file_size = None
  ) -> 'Fastq':
    return self.generate_mock_from_file(
      self.fpath,
      output_fpath,
      base_qty,
      file_size
    )
  
  
  @staticmethod
  def copy(fpath, out_fpath, copy_seq_qty=None):
    # Each read or sequence is 4 lines
    copy_line_qty = 4 * copy_seq_qty
    fname = m.File.get_fname(fpath)

    f_out = open(out_fpath, 'w')
    seq_qty = 0
    print(f"Copying {fname}...")
    # Handle gz files
    if fpath.endswith(".gz"):
      with gzip.open(fpath, "rb") as f:
        i = 1
        for line in f:
          f_out.write(line.decode('utf-8'))
          sys.stdout.write("\r")
          # Sequence line is every 4 lines sequence line starts at line 2
          if (i - 2) % 4 == 0:
            seq_qty += 1
          i += 1
          if (i + 1) % 1_000_000 == 0:
            sys.stdout.write(f" seq_qty: {seq_qty}")
          if i > copy_line_qty:
            break
    
    else:
      with open(fpath, "r") as f:
        i = 1
        for line in f:
          f_out.write(line)
          sys.stdout.write("\r")
          # Sequence line is every 4 lines sequence line starts at line 2
          if (i - 2) % 4 == 0:
            seq_qty += 1
          i += 1
          if (i + 1) % 1_000_000 == 0:
            sys.stdout.write(f" seq_qty: {seq_qty}")
          if i > copy_line_qty:
            break
    print(f" seq_qty: {seq_qty}")

  
  # @staticmethod
  # def analyze(fpath):

  

  ##############################################################################
  # Get Sequence and Base Quantity

  # Example fastq.gz file:
  # b'@A01050:505:HYL7YDSXY:2:1101:12897:1000 1:N:0:GTATTGAC+GGAATTCC\n'
  # b'TTTCAGGTGAGTTGGGGCCTCTAGTTAATGTACCTGGACTCTTTGCCTGTAACTTCTCTGTTTATT\n'
  # b'+\n'
  # b'FFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFFF:FF:FFFFFFFFFF\n'
  # b'@A01050:505:HYL7YDSXY:2:1101:19117:1000 1:N:0:GTATTGAC+GGAATTCC\n'
  # b'AAACAGTGCTTACAACTTCTGCCAACTCATGAAAGCAGAGCCCTGCTGGCAGCCTAATGAAATGCA\n'
  # b'+\n'

  @staticmethod
  def sm_get_base_and_seq_qty(fpath):
    base_qty = 0
    seq_qty = 0
    fname = m.File.get_fname(fpath)

    print(f"Reading {fname}...")
    # Handle gz files
    if fpath.endswith(".gz"):
      with gzip.open(fpath, "rb") as f:
        i = 1
        for line in f:
          sys.stdout.write("\r")
          # Sequence line is every 4 lines sequence line starts at line 2
          if (i - 2) % 4 == 0:
            base_qty += len(line) - 1
            seq_qty += 1
          i += 1
          if (i + 1) % 1_000_000 == 0:
            sys.stdout.write(f" seq_qty: {seq_qty}")
    
    else:
      with open(fpath, "r") as f:
        i = 1
        for line in f:
          sys.stdout.write("\r")
          # Sequence line is every 4 lines sequence line starts at line 2
          if (i - 2) % 4 == 0:
            base_qty += len(line) - 1
            seq_qty += 1
          i += 1
          if (i + 1) % 1_000_000 == 0:
            sys.stdout.write(f" seq_qty: {seq_qty}")
    print(f" seq_qty: {seq_qty}")
    
    return base_qty, seq_qty
  
  def get_base_and_seq_qty(self):
    self.sm_get_base_and_seq_qty(self.fpath)
          
  
  ##############################################################################
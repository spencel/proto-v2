
import gzip
import json

from devtools import debug as devtools_debug

import modules as m


class Debugger():

  devtools_debug = devtools_debug

  def __init__(
    self,
    debug_is_on = False,
    default_print_file_last_line_number = 20,
    log_fpath = None,
  ):

    if type(debug_is_on) is bool:
      self.debug_is_on = debug_is_on
  
    if type(debug_is_on) is int:
      if debug_is_on == 1:
        self.debug_is_on = True
      else:
        self.debug_is_on = False
    
    if type(debug_is_on) is  str:
      if debug_is_on.lower() in ["1", "true"]:
        self.debug_is_on = True
      else:
        self.debug_is_on = False
    
    if self.debug_is_on:
      print("Debug mode is on.")
    
    self.default_print_file_last_line_number = default_print_file_last_line_number

    self.log_fpath = log_fpath

    # Reset file
    if self.log_fpath:
      with open(self.log_fpath, 'w') as f:
        pass

  
  # Will either print to CLI or log file
  def output(self, message):
    if self.log_fpath:
      with open(self.log_fpath, 'a') as f:
          f.write(message + "\n")
    else:
      print(message)


  # Print raw message
  def print(self, message):
    if self.debug_is_on:
      self.output(message)
      

  # Simple print: "var: value"
  def print_var(self, var_name, value):
    if self.debug_is_on:
      self.output(f"{var_name}: {value}")
  
  
  # Prints file lines to see what it looks like
  def print_file_lines(
    self,
    fpath,
    last_line_number = None
  ):
    if not self.debug_is_on:
      return
    
    if not last_line_number:
      last_line_number = self.default_print_file_last_line_number

    self.output(
      f"Printing first {last_line_number} lines of {fpath}"
    )

    def read_lines(f, mode="txt"):
      i = 1
      for line in f:
          match mode:
            case ".gz":
              self.output(line.decode("utf-8").strip("\n"))
            case _:
              self.output(line.strip("\n"))
          i += 1
          if i > last_line_number:
            break

    # Handle gz files
    if fpath.endswith(".gz"):
      with gzip.open(fpath, "rb") as f:
        read_lines(f, ".gz")
    
    # Treat any other file as a text file
    else:
      with open(fpath, 'r') as f:
        read_lines(f, ".txt")
    
    self.output(
      f"\nDone printing first {last_line_number} lines of {fpath}\n"
    )

  @staticmethod
  def write_to_file(
    data,
    fpath = 'debug-dump.txt'
  ):
    
    with open(fpath, 'w') as f:
      f.write(m.Json.json_to_pretty_str(data))
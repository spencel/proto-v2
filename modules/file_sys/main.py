
import hashlib
import os

from .classes import *
import modules as m


def remove_extension(filepath: str) -> str:
  return os.path.splitext(filepath)[0]


def split_path(path: str) -> list[str]:
  path_delim = os.path.sep
  return path.split(path_delim)


def path_list_to_dict(paths: list[str]) -> dict:

  file_tree_dict = dict()

  def recurse(file_tree_dict_branch, path_parts):
    part = path_parts[0]
    if part not in file_tree_dict_branch:
      file_tree_dict_branch[part] = dict()
      if len(path_parts) > 1:
        file_tree_dict_branch[part] = recurse(
          file_tree_dict_branch[part], path_parts[1:]
        )
    elif len(path_parts) > 1:
      file_tree_dict_branch[part] = recurse(
        file_tree_dict_branch[part], path_parts[1:]
      )
    
    return file_tree_dict_branch

  for path in paths:
    path_parts = m.file_sys.split_path(path)
    part = path_parts[0]
    if part not in file_tree_dict:
      file_tree_dict[part] = dict()
      if len(path_parts) > 1:
        file_tree_dict[part] = recurse(file_tree_dict[part], path_parts[1:])
    elif len(path_parts) > 1:
      file_tree_dict[part] = recurse(file_tree_dict[part], path_parts[1:])
  
  return file_tree_dict


global_paths_txt = ''
def path_dict_to_txt(file_tree_dict: dict) -> str:
  global global_paths_txt

  def recurse(file_tree_dict_branch, level):
    global global_paths_txt
    for key in file_tree_dict_branch:
      global_paths_txt += '\t'*level + key + "\n"
      recurse(file_tree_dict_branch[key], level+1)
  
  global_paths_txt = ''
  recurse(file_tree_dict, 0)
  
  return global_paths_txt

# Compare files 2 directories via checksums
def compare_dirs(
  this_dpath: str,
  that_dpath: str
) -> dict:
  
  results = {
    this_dpath: {
      'files': {},
      # 'FILENAME': {
      #   'path': 'FILEPATH',
      #   'checksum': 'CHECKSUM
      # }
      # ...
      'missing': []
      # 'FILENAME', ...
    },
    that_dpath: {
      'files': {},
      'missing': []
    },
    # Checksum matches even filenames differ
    'checksum_matches': [],
    # {
    #   'this_fname': 'FILENAME',
    #   'that_fname': 'FILENAME',
    #   'checksum': 'CHECKSUM'
    # }
    'checksum_mismatches': []
    # 'FILENAME', ...
  }

  def get_checksum(fpath):
    print(f' hashing {fpath}')
    hasher = hashlib.sha256()
    with open(fpath, 'rb') as f:
      for chunk in iter(lambda: f.read(4096), b''):
        hasher.update(chunk)
    return hasher.hexdigest()
  
  def get_data(dpath):
    data = {} # fname: fpath, checksum
    for name in os.listdir(dpath):
      path = os.path.join(dpath, name)
      if os.path.isfile(path):
        checksum = get_checksum(path)
        data[name] = {
          'path': path,
          'checksum': checksum
        }
    return data

  these_fnames = get_data(this_dpath)
  results[this_dpath]['files'] = these_fnames
  
  those_fnames = get_data(that_dpath)
  results[that_dpath]['files'] = those_fnames
  
  # Get filenames and missing filenames
  for this_fname in these_fnames:
    if this_fname not in those_fnames:
      results[that_dpath]['missing'].append(this_fname)
  for that_fname in those_fnames:
    if that_fname not in these_fnames:
      results[this_dpath]['missing'].append(that_fname)
  
  # Get checksum mismatches by file name
  for this_fname, data in these_fnames.items():
    this_checksum = data['checksum']
    if this_fname in those_fnames:
      that_checksum = those_fnames[this_fname]['checksum']
    else:
      continue
    if this_checksum != that_checksum:
      results['checksum_mismatches'].append(this_fname)
  
  # Get all checksum matches even filenames differ
  for this_fname, data in these_fnames.items():
    this_checksum = data['checksum']
    for that_fname, data in those_fnames.items():
      that_checksum = data['checksum']
      if this_checksum == that_checksum:
        results['checksum_matches'].append({
          'this_fname': this_fname,
          'that_fname': that_fname,
          'checksum': this_checksum
        })

  return results


def get_filename(filepath, with_extension=True):
  filename = os.path.basename(filepath)
  if with_extension:
    return filename
  else:
    return os.path.splitext(filename)[0]


def make_dirs(dpath):
  if not os.path.exists(dpath):
    os.makedirs(dpath)
    return 1
  else:
    return 0
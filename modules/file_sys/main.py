
import os

import modules as m


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
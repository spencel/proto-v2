
from datetime import datetime
import json
import traceback

import modules as m


# Static
def load_as_dot_notation(
  json_path,
  is_mutable = True
):
  return m.Dict.create_dotdict(json.load(open(json_path, 'r')), is_mutable)


# Open a json file and convert its json a dictionary/tuple
def load_file(fpath):
  try:

    json_data = None
    with open(fpath, 'r') as f:
      json_data = json.load(f)
    return json_data

  # Catch any exception
  except Exception as e:
    print('\nUnexpected exception:')
    print(traceback.format_exc())
    raise e


# Used mostly for debugging
def to_pretty_str(arg_json):
  try:
    return json.dumps(arg_json, indent=2)
  except:
    # logging.info("not json serializable")
    # Line of code below doesn't work yet?
    # printStr = '%s: [not json serializable] %s' % (printStr, json.dumps(make_jsonable(value), indent=2))
    return arg_json


def to_table(
  json_data: dict,
  col_names: dict = {'key_col_name': None},
  print_header_row: bool = True
) -> list[list]:
  
  # Set key column name to key if not given
  key_col_name = col_names['key_col_name'] if col_names['key_col_name'] else 'key'

  col_names_map = {
    key_col_name: 0
  }
  col_idx_map = {
    0: key_col_name
  }
  col_qty: int

  # Get column names from first item in dictionary
  for key, items in json_data.items():

    for col_idx, subkey in enumerate(items, start=1):
      col_name = subkey
      # Get column index
      col_names_map[col_name] = col_idx
      col_idx_map[col_idx] = col_name

    break

  col_qty = len(col_names_map)

  table_lines = []

  # Add header row
  if print_header_row:
    new_line = [None] * col_qty
    for col_idx, col_name in col_idx_map.items():
      new_line[col_idx] = col_name
    table_lines.append(new_line)

  # Generate table
  for i_line, (key, items) in enumerate(json_data.items()):
    new_line = [None] * col_qty
    new_line[0] = key

    for col_name, value in items.items():
      col_idx = col_names_map[col_name]
      new_line[col_idx] = value
    
    table_lines.append(new_line)

  return table_lines
   


@staticmethod
def _custom_encoder(obj):
  if isinstance(obj, datetime):
    return obj.isoformat()
  

# Save to file
def save_to_file(json_data, fpath, indent=2):
  with open(fpath, 'w') as f:
    json.dump(
      json_data,
      f,
      indent = indent,
      default = _custom_encoder
      )
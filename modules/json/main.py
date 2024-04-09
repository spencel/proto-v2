
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
  col_names: dict = {'key_col_name': None}
) -> list[list]:
  
  key_col_name = col_names['key_col_name'] if col_names['key_col_name'] else 'key'

  col_idxs = {
    key_col_name: 0
  }
  col_qty: int

  table_lines = [[key_col_name]]
  
  # Get define column names from first item in dictionary
  for line_idx, (key, subdict) in enumerate(json_data.items(), start=1):
    table_lines.append([key])

    for col_idx, (subkey, value) in enumerate(subdict.items(), start=1):
      col_idxs[subkey] = col_idx
      col_name = subkey
      table_lines[0].append(col_name)
      table_lines[line_idx].append(value)

    break

  col_qty = len(col_idxs)
  
  # Continue generating the rest of the lines
  for line_idx, (key, subdict) in enumerate(json_data.items(), start=2):
    print(f'col_qty: {col_qty}')
    new_line = [None] * col_qty
    new_line[0] = key

    print(f'new_line: {new_line}')
    table_lines.append(new_line)
    print(f'line_idx: {line_idx}')
    print(table_lines)

    for subkey, value in subdict.items():
      col_name = subkey
      col_idx = col_idxs[col_name]
      table_lines[line_idx][col_idx] = (value)

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
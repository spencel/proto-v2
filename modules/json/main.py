
from datetime import datetime
import json
import traceback

import modules as m


# Static
def load_json_as_dot_notation(
  json_path,
  is_mutable = True
):
  return m.Dict.create_dotdict(json.load(open(json_path, 'r')), is_mutable)


# Open a json file and convert its json a dictionary/tuple
def load_json_file(fpath):
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
def json_to_pretty_str(arg_json):
  try:
      return json.dumps(arg_json, indent=2)
  except:
      # logging.info("not json serializable")
      # Line of code below doesn't work yet?
      # printStr = '%s: [not json serializable] %s' % (printStr, json.dumps(make_jsonable(value), indent=2))
      return arg_json


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
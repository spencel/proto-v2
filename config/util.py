
import json


class _dotdict_mutable(dict):
  __getattr__ = dict.get
  __setattr__ = dict.__setitem__
  __delattr__ = dict.__delitem__


class _dotdict_constant(dict):
  __getattr__ = dict.get

def create_dotdict(data, is_mutable=True):
    if is_mutable:
      data = _dotdict_mutable(data)
    else:
      data = _dotdict_constant(data)
    for key in data.keys():
      if isinstance(data[key], dict):
        data[key] = create_dotdict(data[key], is_mutable)
    return data

def load_json_as_dot_notation(json_path, is_mutable=True):
    return create_dotdict(json.load(open(json_path, 'r')), is_mutable)
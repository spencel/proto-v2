

class _dotdict_mutable(dict):
  __getattr__ = dict.get
  __setattr__ = dict.__setitem__
  __delattr__ = dict.__delitem__


class _dotdict_constant(dict):
  __getattr__ = dict.get


def create_dotdict(
  data: dict|None = None,
  is_mutable: bool =True
) -> _dotdict_mutable|_dotdict_constant:
  """_summary_

  Args:
      data (dict | None, optional): _description_. Defaults to None.
      is_mutable (bool, optional): _description_. Defaults to True.

  Returns:
      _dotdict_mutable|_dotdict_constant: _description_
  """
  if not data:
    data = dict()
  
  if is_mutable:
    data = _dotdict_mutable(data)
  else:
    data = _dotdict_constant(data)
  for key in data.keys():
    if isinstance(data[key], dict):
      data[key] = create_dotdict(data[key], is_mutable)
  return data
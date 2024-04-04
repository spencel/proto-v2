




class _dotdict_mutable(dict):
  __getattr__ = dict.get
  __setattr__ = dict.__setitem__
  __delattr__ = dict.__delitem__


class _dotdict_constant(dict):
  __getattr__ = dict.get


class Dict():


  @classmethod
  def create_dotdict(cls, data, is_mutable=True):
    if is_mutable:
      data = _dotdict_mutable(data)
    else:
      data = _dotdict_constant(data)
    for key in data.keys():
      if isinstance(data[key], dict):
        data[key] = cls.create_dotdict(data[key], is_mutable)
    return data
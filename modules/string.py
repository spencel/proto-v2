
import config
import modules as m


class String():
    

  @staticmethod
  def try_str_to_number(string) -> int|float|str:
    # Check if it's an integer
    try:
      string_as_int = int(string.replace(',', ''))
      return string_as_int
    except:
      pass
    
    # Check if it's a float
    # This must come after the int() check, because a string representing aa
    # integer can also be cast as a float.
    try:
      string_as_float = float(string.replace(',', ''))
      return string_as_float
    except:
      pass
    
    # String isn't a number
    return string
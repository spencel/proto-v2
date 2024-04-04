
conversions = {

  # Digital Computing
  # Units of information
  'bit': 1,
  'byte': 8, # bits
  'kilobyte': 8 * 1_000,
  'kibibyte': 8 * 1_0024,
  'megabyte': 8 * 1_000_000,
  'mebibyte': 8 * 1_048_576,
  'gigabyte': 8 * 1_000_000_000,
  'gibibyte': 8 * 1_073_741_824

}

def convert(qty, from_unit, to_unit):
  # Eg 3000 bytes to kilobytes = 3 kilobytes
  # result = 3000 / (8 * 1_000) * 8  = 3
  # result = qty / conversions['kilobyte'] * conversions['byte'] = 3
  # result = qty / conversions[to_unit] * conversions[from_unit] = 3
  result = qty / conversions[to_unit] * conversions[from_unit]
  return result
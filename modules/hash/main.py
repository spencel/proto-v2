
import hashlib


SHA256 = 'sha256'
HASHE_TYPES = (SHA256)


def get_hash(
  input: str,
  hash_type: str = SHA256
) -> str:
  """_summary_

  Args:
      input (str, required): _description_.
      hash_type (str, optional): _description_. Defaults to 'sha256'.

  Returns:
      str: hex digest
  """
  hash_str: str
  
  if hash_type == SHA256:
    hash = hashlib.sha256()
    hash.update(input.encode())
    hash_str = hash.hexdigest()
  
  return hash_str
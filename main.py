
# import logging

import commands
import modules as m
from modules.test import Class_A


# logging.basicConfig(
#   filename = 'logs/main.log',
#   filemode = 'w',
#   level =logging.DEBUG
# )


def test():
  m.test.test()
  m.test.Class_A.test()
  Class_A.test()
  m.test.Class_A.Class_C.test()
  m.test.Class_A.Class_D.test()
  m.test.Class_B.test()
  m.test.helpers.test()
  m.test.helper_scripts.test()
  m.test.scripts.test()
  


if __name__ == "__main__":
  # python -m fire commands/fastq/main.
  pass
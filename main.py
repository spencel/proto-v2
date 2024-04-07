
import commands
import modules as m
from modules.test import Class_A

def test():
  m.test.test()
  m.test.Class_A.test()
  Class_A.test() # from .main import *
  m.test.Class_A.Class_C.test()
  m.test.Class_A.Class_D.test()
  m.test.Class_B.test()

if __name__ == "__main__":
  # python -m fire commands/fastq/main.
  pass
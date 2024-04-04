
import commands
import modules

def test():
  modules.test.test()
  modules.test.Class_A.test()
  modules.test.Class_A.Class_C.test()
  modules.test.Class_A.Class_D.test()
  modules.test.Class_B.test()

if __name__ == "__main__":
  # python -m fire commands/fastq/main.
  pass

import modules as m


class Log():

    DEFAULT_LOG_FPATH = "./default.log"

    @classmethod
    def reset(cls, fpath=None):
        if not fpath:
            fpath = cls.DEFAULT_LOG_FPATH
        with open(fpath, 'w') as f:
            pass


    @classmethod
    def write(cls, text, fpath=None):
        if not fpath:
            fpath = cls.DEFAULT_LOG_FPATH
        with open(fpath, 'a') as f:
            f.write(text+"\n")
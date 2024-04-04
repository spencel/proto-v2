
import os
import re
import shutil


class Dir():


    @staticmethod
    def make(dpath):
        os.mkdir(dpath)

    @staticmethod
    # make_dir_if_not_exist
    def make_if_not_exist(dpath):
        """
        Make a directory if it doesn't exist and required parent directories
        Returns:
        True  - if the directory didn't exist
        False - if the directory did exist
        :param dpath:
        :return:
        """
        if not os.path.isdir(dpath):
            os.makedirs(dpath)
            return True
        else:
            return False


    @staticmethod
    def delete(dpath):
        """
        This will delete a directory and its contents if it exists.
        :param dpath:
        :return:
        """
        if os.path.isdir(dpath):
            shutil.rmtree(dpath)


    @classmethod
    def replace(cls, dpath):
        cls.delete(dpath)
        cls.make(dpath)

    @staticmethod
    def get_fpaths(dpath, regular_expression=".", recursive=False) -> list[str]:
        """
        Returns a list of filepaths in the specified directory.
        :param dpath:
        :return:
        """
        fpaths = []

        if recursive:
            for root, dnames, fnames in os.walk(dpath):
                for fname in fnames:
                    if re.search(regular_expression, fname):
                        fpaths.append(os.path.join(
                            root, fname
                        ))

        else:
            for item in os.listdir(dpath):
                path = os.path.join(dpath, item)
                if os.path.isfile(path) \
                        and re.search(regular_expression, path):
                    fpaths.append(path)
                # elif os.path.isdir(path) and recursive:
                #     fpaths = recurse(fpaths, path)

        return fpaths
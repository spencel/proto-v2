
import re


# Regular Expressions
class RegEx():


    @classmethod
    def clean(cls, raw_string):
        """
        If a raw string is used as a match pattern, escapes necessary characters.
        :param raw_string:
        :return:
        """
        return raw_string.replace(".", "\.")


    # Todo: move this to a different class
    @classmethod
    def get_who_label_from_gisaid_variant(cls, covv_variant):
        who_label = re.sub("^[A-Z]* ([A-Za-z0-9/]*) .*", "\\1", covv_variant)
        if not who_label:
            return None
        return who_label
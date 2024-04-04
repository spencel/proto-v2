
import datetime
import logging
import os

# Configs



class Time():

    @staticmethod
    def get_one_year_ago():
        return datetime.date.today() - datetime.timedelta(days=365)
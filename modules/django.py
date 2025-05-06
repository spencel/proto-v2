
import django
from django.conf import settings

from config import db as db_configs

class Django():

  @staticmethod
  def init():
    if settings.configured:
      pass
    else:
      # logging.debug(f"{scope_pref}: Django setup")
      db_config = db_configs.postgresql.databases.proto
      settings.configure(
        INSTALLED_APPS=[
          'modules',
        ],
        DATABASES={
          'default': {
            "ENGINE": db_config.django_engine,
            "NAME": db_config.database,
            "USER": db_config.user,
            "PASSWORD": db_config.password,
            "HOST": db_config.host,
            "PORT": db_config.port,
          }
        }
      )
      django.setup()

import django
from django.conf import settings

import config as c

class Django():

  @staticmethod
  def init():
    if settings.configured:
      pass
    else:
      # logging.debug(f"{scope_pref}: Django setup")
      settings.configure(
        INSTALLED_APPS=[
          'modules',
        ],
        DATABASES={
          'default': {
            "ENGINE": c.db.iron_scratch.django_engine,
            "NAME": c.db.iron_scratch.database,
            "USER": c.db.iron_scratch.user,
            "PASSWORD": c.db.iron_scratch.password,
            "HOST": c.db.iron_scratch.host,
            "PORT": c.db.iron_scratch.port,
          }
        }
      )
      django.setup()
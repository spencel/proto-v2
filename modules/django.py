
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
            "ENGINE": c.db.proto.django_engine,
            "NAME": c.db.proto.database,
            "USER": c.db.proto.user,
            "PASSWORD": c.db.proto.password,
            "HOST": c.db.proto.host,
            "PORT": c.db.proto.port,
          }
        }
      )
      django.setup()
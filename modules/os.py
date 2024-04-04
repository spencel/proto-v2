
import os
import re
import shutil
import uuid

from django.db import models


# Configs



class Os(models.Model):

    class Meta:
        db_table = "file_system"

    id = models.UUIDField(
        primary_key=True,
        default=uuid.uuid4,
    )
    filepath = models.TextField()




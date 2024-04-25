
import os
from pathlib import Path

from .util import *

from .ncbi import *

# Config
root_dpath = os.getcwd()
this_dpath = Path(__file__).parent


# Classes/Functions
# Temporary fix, could do some of these configs in their respective modules or
# create modules where they would belong.




# Get paths
paths = load_json_as_dot_notation(
  os.path.join(this_dpath, "paths.json"),
  is_mutable=True
)
paths.genomes = os.path.join(
  paths.dir_names.data,
  paths.dir_names.genome_data
)
paths.gisaid_data = os.path.join(
  paths.dir_names.data,
  paths.dir_names.gisaid_data
)
paths.blastdb = os.path.join(
  paths.dir_names.data,
  paths.dir_names.blast_databases
)

# Get environments
env_dpath = os.path.join(
  root_dpath,
  paths.dir_names.environment
)
db = load_json_as_dot_notation(
  os.path.join(env_dpath, "db.json"),
  is_mutable=False
)

# Set aliases
pg_schema_dpath = os.path.join(
  root_dpath,
  paths.dir_names.postgres_database
)

# Get BLASTN json
blast = load_json_as_dot_notation(
  os.path.join(this_dpath, "blast.json"),
  is_mutable=False
)

# Get AWS
aws = load_json_as_dot_notation(
  os.path.join(
    this_dpath,
    'aws',
    'config.json'
  ),
  is_mutable = False
)
aws.data_dpath = os.path.join(
   paths.dir_names.data,
   paths.dir_names.modules,
   aws.module_name
)

# Get Basepair API
basepair_api = json.load(open(os.path.join(
   this_dpath,
   'basepair_api',
   'sl-basepair.json'
), 'r'))


# Setup logging
logging = load_json_as_dot_notation((os.path.join(
   this_dpath,
   'logging',
   'config.json'
)))
logging.default_log_fpath = os.path.join(
	logging.log_dname,
	logging.default_log_fname
)
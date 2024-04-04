
import json
import os

import basepair as basepair

import config as cf
import modules as m


class Basepair():
  
  data_dpath = os.path.join(
    'data', 'modules', 'basepair_api'
  )
  analyses_fpath = os.path.join(
    data_dpath, 'basepair-analyses.json'
  )
  projects_fpath = os.path.join(
    data_dpath, 'basepair-projects.json'
  )
  samples_fpath = os.path.join(
    data_dpath, 'basepair-samples.json'
  )
  metrics_fpath = os.path.join(
    data_dpath, 'metrics.json'
  )

  # Constructor
  def __init__(self,
    api_config = cf.basepair_api,
    is_connect: bool = False
  ):
    self.api_config = api_config
    self.api = None
    if is_connect:
      self.api = basepair.connect(api_config)
  

  def connect(self, 
    api_config = cf.basepair_api
  ):
    self.api = basepair.connect(api_config)


  @staticmethod
  def _get_user_id_from_path(path: str) -> str:
    # Example path: '/api/v2/users/5371'
    user_id = path.split('/')[-1]
    return user_id


  @staticmethod
  def sm_redo_analysis(
    analysis_id,
    instance_type = None,
    basepair_config = None,
    basepair_api = None
  ):
    if not basepair_config:
      basepair_config = cf.basepair_api
    if not basepair_api:
      basepair_api = basepair.connect(basepair_config)
    
    basepair_api.restart_analysis(
      uid = analysis_id,
      instance_type = instance_type
    )


  def get_metrics(self):

    metrics = dict()
    '''
    {
      'user_id': {
        'username': str
        'project_qty': int
        'analysis_qty': int
        'sample_qty': int
        'sample_storage_GB': float (GB)
      }
    }
    '''
    def init_data(username: str = '') -> dict:
      return {
        'username': username,
        'project_qty': 0,
        'analysis_qty': 0,
        'sample_qty': 0,
        'sample_storage_GB': 0
      }

    # projects = self.api.get_projects()
    # with open(self.projects_fpath, 'w') as f:
    #   f.write(json.dumps(projects, indent=2))
    projects =  json.load(open(self.projects_fpath))

    # analyses = self.api.get_analyses()
    # with open(self.analyses_fpath, 'w') as f:
    #   f.write(json.dumps(analyses, indent=2))
    analyses =  json.load(open(self.analyses_fpath))

    # samples = self.api.get_samples()
    # with open(self.samples_fpath, 'w') as f:
    #   f.write(json.dumps(samples, indent=2))
    samples =  json.load(open(self.samples_fpath))

    for project in projects:
      user_id = self._get_user_id_from_path(project['owner'])
      if user_id not in metrics:
        username = project['owner__username']
        metrics[user_id] = init_data(username)
      metrics[user_id]['project_qty'] += 1
  
    for analysis in analyses:
      user_id = self._get_user_id_from_path(analysis['owner'])
      if user_id not in metrics:
        username = analysis['owner__username']
        metrics[user_id] = init_data(username)
      metrics[user_id]['analysis_qty'] += 1
    
    for sample in samples:
      user_id = self._get_user_id_from_path(sample['owner'])
      if user_id not in metrics:
        username = sample['owner__username']
        metrics[user_id] = init_data(username)
      metrics[user_id]['sample_qty'] += 1
      filesize_GB = m.units.convert(
        sample['meta']['filesize'],
        'byte',
        'gigabyte'
      )
      metrics[user_id]['sample_storage_GB'] += filesize_GB
    
    with open(self.metrics_fpath, 'w') as f:
      f.write(json.dumps(metrics, indent=2))

    return metrics
    

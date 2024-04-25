
import boto3

from modules import logging
from modules import dict

import config as cf


class Ec2():
  
  log = logging.Log(__qualname__)

  def __init__(self,
    aws_profile_name: str = cf.aws.mbsl,
  ):
    self.aws_profile_name = aws_profile_name
    self.log.write(None, 'aws_profile', aws_profile_name)
    self.session = boto3.Session(
      profile_name = aws_profile_name
    )
    self.client = self.session.client('ec2')

  @staticmethod
  def get_launch_spec(
    max_inst_qty: int = 1,
    min_inst_qty: int = 1,
    launch_template_id: str|None = None,
    launch_template_name: str|None = None,
    launch_template_version: str|None = None
  ):
    
    launch_spec = {}

    # LaunchTemplate
    launch_spec['launch_template'] = {}
    if launch_template_id:
      launch_spec['launch_template']['LaunchTemplateId'] = launch_template_id
    if launch_template_name:
      launch_spec['launch_template']['LaunchTemplateName'] = launch_template_name
    if launch_template_version:
      launch_spec['launch_template']['Version'] = launch_template_version

    # Instance Quantity
    # MaxCount
    launch_spec['max_inst_qty'] = max_inst_qty
    # MinCount
    launch_spec['min_inst_qty'] = min_inst_qty

    launch_spec = dict.create_dotdict(launch_spec)
    
    return launch_spec

    
    
  def request_spot_instance(self,
    launch_template_id: str = cf.aws.ec2.default_launch_template_id
  ):
    
    launch_spec = self.get_launch_spec(
      launch_template_id = launch_template_id,
      launch_template_version = '1'
    )
    self.log.write(launch_spec, 'launch_spec')
      
    res = self.client.run_instances(
      LaunchTemplate = launch_spec.launch_template,
      MaxCount = launch_spec.max_inst_qty,
      MinCount = launch_spec.min_inst_qty
    )
    self.log(res, 'res')
    spot_inst_req_id = res['SpotInstanceRequests'][0]['SpotInstanceRequestId']

    self.log.write(spot_inst_req_id, 'spot_inst_req_id')
  

  def describe_instances(self):
    # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/opsworks/client/describe_instances.html#describe-instances
    res = self.client.describe_instances()
    self.log.write('describe_instances', 'res', res)
  
  def describe_regions(self):
    res = self.client.describe_regions()
    self.log.write('describe_regions', 'res', res)
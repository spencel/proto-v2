
import time

import boto3

from modules import logging
from modules import dict
from modules.aws.classes.S3 import S3

import config as cf


# https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/ec2.html
class Ec2():
  
  log = logging.Log(__qualname__)

  def __init__(self,
    aws_profile_name: str = cf.aws.sdk.profile_names.micronbrane_spencer,
  ):
    self.aws_profile_name = aws_profile_name
    self.log.debug(None, 'aws_profile', aws_profile_name)
    self.session = boto3.Session(
      profile_name = aws_profile_name
    )
    self.client = self.session.client('ec2')

  
  def make_network_interfaces_obj(self):
    pass


  def get_launch_spec(self,
    # IAM Instance Profile
    instance_profile_iam_role_arn: str|None = None,
    instance_profile_iam_role_id: str|None = None,
    # Instance Type
    instance_type: str|None = None,
    # Instance Market Options
    market_type: str|None = None,
    # Tag Specifications
    instance_name: str|None = None,
    # Instance Quantity
    max_inst_qty: int = 1,
    min_inst_qty: int = 1,
    # Launch Template
    launch_template_id: str|None = None,
    launch_template_name: str|None = None,
    launch_template_version: str|None = None,
    # User Data
    user_data: str|None = None
  ):
    
    launch_spec = {}

    if instance_profile_iam_role_arn or instance_profile_iam_role_id:
      launch_spec['iam_instance_profile'] = {}
      if instance_profile_iam_role_arn:
        launch_spec['iam_instance_profile']['Arn'] = instance_profile_iam_role_arn
      if instance_profile_iam_role_id:
        launch_spec['iam_instance_profile']['Id'] = instance_profile_iam_role_id

    # Instance Type
    if instance_type:
      launch_spec['instance_type'] = instance_type

    # Launch Template
    if launch_template_id or launch_template_name or launch_template_version:
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

    # Instance Market Options
    if market_type:
      launch_spec['instance_market_options'] = {
        'MarketType': market_type
      }

    # Tag Specifications
    # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/ec2/client/request_spot_instances.html#:~:text=request%20was%20created.-,TagSpecifications,-(list)%20%E2%80%93
    # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/ec2/client/create_tags.html#create-tags
    # For spot instances, must use client.create_tags() after request_spot_instances()
    # The instance Name cannot be set from the Launch Template
    # https://ap-southeast-1.console.aws.amazon.com/ec2/home?region=ap-southeast-1#LaunchTemplates:
    if instance_name:
      launch_spec['tag_specs'] = []
      launch_spec['tag_specs'].append({
        'ResourceType': 'spot-instances-request',
        'Tags': [{
          'Key': 'Name',
          'Value': instance_name
        }]
      })

    # User Data
    if user_data:
      launch_spec['user_data'] = user_data

    launch_spec = dict.create_dotdict(launch_spec)
    
    return launch_spec
  
  # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/ec2/client/terminate_instances.html
  def terminate_instances(self,
    instance_ids = list
  ):
    res = self.client.terminate_instances(
      InstanceIds = instance_ids
    )
    self.log.debug('terminate_instances', 'res', res)
    return res

  # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/opsworks/client/describe_instances.html#describe-instances
  def describe_instances(self,
    instance_ids: list|None = None
  ):
    res = self.client.describe_instances(
      InstanceIds = instance_ids
    )
    self.log.debug('describe_instances', 'res', res)
    return res
  

  def describe_regions(self):
    res = self.client.describe_regions()
    self.log.debug('describe_regions', 'res', res)
    return res
  

  # Don't use client.request_spot_instances(), it's deprecated
  # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/ec2/client/run_instances.html#
  # Launching a spot instance: https://ec2spotworkshops.com/launching_ec2_spot_instances/runinstances_api.html
  # Spot Instance Advisor: https://aws.amazon.com/ec2/spot/instance-advisor/
  def request_spot_instances_v1(self,
    # IAM Instance Profile
    instance_profile_iam_role_arn: str = cf.aws.ec2.default_instance_profile_iam_role_arn,
    # Instance Type
    instance_type: str = cf.aws.ec2.default_instance_type,
    # Launch Template
    launch_template_id: str = cf.aws.ec2.default_launch_template_id,
    # Tag Specifications
    instance_name: str = cf.aws.ec2.default_instance_name
  ):
    
    # Transferring files from S3 to EC2
    # https://docs.aws.amazon.com/AWSEC2/latest/UserGuide/AmazonS3.html
    user_data = str(
      'curl "https://awscli.amazonaws.com/awscli-exe-linux-x86_64.zip" -o "awscliv2.zip"\n'
      'unzip awscliv2.zip\n'
      'sudo ./aws/install\n'
      'aws s3 cp s3://micronbrane.test/Tax_cate_230511.txt /home/Tax_cate_230511.txt\n'
      'aws s3 cp /home/Tax_cate_230511.txt s3://micronbrane.test/test/Tax_cate_230511.txt\n'
      # 'sudo shutdown --poweroff now -h\n'
    )
    
    launch_spec = self.get_launch_spec(
      # IAM Instance Profile
      instance_profile_iam_role_arn = instance_profile_iam_role_arn,
      # Instance Type
      instance_type = instance_type,
      # Instance Market Options
      market_type = 'spot',
      # Launch Template
      launch_template_id = launch_template_id,
      # User Data
      user_data = user_data
    )
    self.log.debug('request_spot_instances', 'launch_spec', launch_spec)
    
    run_instances_res = self.client.run_instances(
      IamInstanceProfile = launch_spec.iam_instance_profile,
      InstanceMarketOptions = launch_spec.instance_market_options,
      InstanceType = launch_spec.instance_type,
      LaunchTemplate = launch_spec.launch_template,
      MaxCount = launch_spec.max_inst_qty,
      MinCount = launch_spec.min_inst_qty,
      UserData = launch_spec.user_data
    )
    self.log.debug('request_spot_instances', 'run_instances_res', run_instances_res)

    instance_id = run_instances_res['Instances'][0]['InstanceId']

    create_tags_res = self.client.create_tags(
      Resources = [instance_id],
      Tags = [{
        'Key': 'Name',
        'Value': instance_name
      }]
    )
    self.log.debug('request_spot_instances', 'create_tags_res', create_tags_res)

    # Wait for spot instance to launch
    is_waiting = True
    while is_waiting:
      # Wait 5 seconds
      time.sleep(5)
      self.log.info('request_spot_instances', value='Instance is pending.')

      # Get instance state
      # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/ec2/client/run_instances.html#:~:text=instance%2C%20if%20applicable.-,State%20(dict)%20%E2%80%93,-The%20current%20state
      describe_instances_res = self.describe_instances(
        instance_ids = [instance_id]
      )
      state_name: str
      for instance in describe_instances_res['Reservations'][0]['Instances']:
        if instance['InstanceId'] == instance_id:
          state_name = instance['State']['Name']
      
      # Check if the instance is running
      if state_name == cf.aws.ec2.instance.states.running.name:
        is_waiting = False
        self.log.info('request_spot_instances', value='Instance is running.')
      
      # If the instance is neither pending or running, terminate it
      elif state_name != cf.aws.ec2.instance.states.pending.name:
        # Terminate instance
        self.log.info('request_spot_instances', value='Something went wrong, terminating instance.')
        terminate_instances_res = self.terminate_instances(
          instance_ids = [instance_id]
        )
        is_waiting = False

    

    # # Terminate instance
    # self.log.info('request_spot_instances', value='Terminating instance...')
    # terminate_instances_res = self.terminate_instances(
    #   instance_ids = [instance_id]
    # )
  

  # https://boto3.amazonaws.com/v1/documentation/api/latest/reference/services/ec2/client/run_instances.html#
  def request_spot_instances_v2(self):

    run_instances_res = self.client.run_instances(
      # ImageId = 'ami-0e97ea97a2f374e3d', # Default Amazon Linux 2023
      LaunchTemplate = {'LaunchTemplateId': 'lt-0ad2f9d5c24f55c68'},
      InstanceType = 't2.micro',
      # InstanceType = 'd3.4xlarge',
      IamInstanceProfile = {'Arn': 'arn:aws:iam::255446977961:instance-profile/test.cubic_station'},
      MaxCount = 1,
      MinCount = 1,
      SecurityGroupIds = ['sg-054a5c1b2428488ca'],
      InstanceMarketOptions = {
        'MarketType': 'spot'
      }
    )
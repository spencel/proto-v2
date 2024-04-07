
# Standard library
import json
import os
import subprocess
import sys

# 3rd Party
import boto3

# Project specific
import config as cf
import modules as m

# Class specific
from .methods import *


# Constants

# Get Bucket Objects as File Tree
TYPE_TREE = 'tree'
TYPE_LIST = 'list'


# Class S3
class S3():


  _service_name = cf.aws.sdk.clients.s3
  _data_dpath = os.path.join(
    cf.aws.data_dpath,
    _service_name
  )


  # Constructor
  def __init__(self,
    is_connect: bool = False,
    profile_name: str = cf.aws.sdk.profile_names.default
  ):
    self.client = None
    if is_connect:
      self.session = boto3.session.Session(
        profile_name = profile_name
      )
      self.client = self.session.client(
        service_name = self._service_name
      )

  
  # Get Buckets
  def get_buckets(self,
    is_save: bool = False,
    save_data_fpath: str|None = None
  ):
    """
    Example command: python -m fire main modules aws Aws S3 --is-connect True --profile-name mbsl  get_buckets --is-save True
    """
    res = self.client.list_buckets()
    if is_save:
      if not save_data_fpath:
        save_data_fpath = os.path.join(
          self._data_dpath,
          'bucket-list.json'
        )
      m.json.save_to_file(res, save_data_fpath)

    return res
  
  
  # Get Bucket Objects
  def get_bucket_objects(self,
    bucket_name: str,
    is_save: bool = False,
    save_data_fpath: str|None = None
  ):
    responses = []
    responses.append(self.client.list_objects_v2(
      Bucket = bucket_name
    ))
    i = 0
    while 'NextContinuationToken' in responses[i]:
      ContinuationToken = responses[i]['NextContinuationToken']
      responses.append(self.client.list_objects_v2(
        Bucket = bucket_name,
        ContinuationToken = ContinuationToken
      ))
      i = len(responses) - 1
    
    if is_save:
      if not save_data_fpath:
        save_data_fpath = os.path.join(
          self._data_dpath,
          f'{bucket_name}.bucket-objects.json'
        )
      m.json.save_to_file(responses, save_data_fpath)

    return responses


  # Get Bucket Objects as File Tree
  def get_bucket_file_tree(self,
    bucket_name: str,
    type: str = TYPE_TREE,
    is_save: bool = False,
    save_data_fpath: str|None = None
  ):
    
    responses = self.get_bucket_objects(
      bucket_name
    )

    filepaths: str|list
    fname = ''

    if type == TYPE_LIST:
      fname = f'{bucket_name}.file-list.txt'
      filepaths = ''
      for res in responses:
        for contents in res['Contents']:
          filepaths += contents['Key'] + '\n'
      
    
    elif type == TYPE_TREE:
      fname = f'{bucket_name}.file-tree.txt'
      filepaths = []
      for res in responses:
        for contents in res['Contents']:
          filepaths.append(contents['Key'])
    
    def write_file_structure(file_paths, save_data_fpath):
      with open(save_data_fpath, 'w') as f:
          for path in file_paths:
              components = path.split('/')
              current_indent = 0
              for component in components:
                  f.write('  ' * current_indent + component + '\n')
                  current_indent += 1
    
    # def print_file_tree(directory, file_output):
    #   with open(file_output, "w") as f:
    #     for root, dirs, files in os.walk(directory):
    #       level = root.replace(directory, '').count(os.sep)
    #       indent = ' ' * 4 * (level)
    #       f.write('{}{}/\n'.format(indent, os.path.basename(root)))
    #       sub_indent = ' ' * 4 * (level + 1)
    #       for file in files:
    #         f.write('{}{}\n'.format(sub_indent, file))

    if is_save:
      if not save_data_fpath:
        save_data_fpath = os.path.join(
          self._data_dpath,
          fname
        )

      if type == TYPE_LIST:
        with open(save_data_fpath, 'w') as f:
          f.write(filepaths)

      elif type == TYPE_TREE:
          
        file_tree_dict = m.file_sys.path_list_to_dict(
          filepaths
        )
          
        m.json.save_to_file(
          file_tree_dict,
          os.path.join(
            self._data_dpath,
            f'{bucket_name}.file-tree.json'
          )
        )

        paths_txt = m.file_sys.path_dict_to_txt(file_tree_dict)
        with open(save_data_fpath, 'w') as f:
          f.write(paths_txt)


    return file_tree_dict


  # Upload File in Parts
  upload_multipart_file = upload_multipart_file


  # Get Multipart Uploads
  @staticmethod
  def get_multipart_uploads(
    profile_name: str,
    bucket_name: str
  ):
    """get_multipart_uploads.

    :profile_name str:
    :bucket_name str: 

    :return: Response JSON converted to dict.
    :rtype: dict
    """

    proc = subprocess.run(
      'aws s3api list-multipart-uploads'
      f' --profile {profile_name}'
      f' --bucket {bucket_name}',
      shell = True,
      stdout = subprocess.PIPE
    )

    return json.loads(proc.stdout)
  

  @staticmethod
  def del_multipart_upload(
    profile_name: str,
    bucket_name: str,
    key: str,
    upload_id: str
  ):
    """del_multipart_upload.
    https://awscli.amazonaws.com/v2/documentation/api/latest/reference/s3api/abort-multipart-upload.html

    :profile_name str:
    :bucket_name str: 

    :return: Response JSON converted to dict.
    :rtype: dict
    """

    proc = subprocess.run(
      'aws s3api abort-multipart-upload'
      f' --profile {profile_name}'
      f' --bucket {bucket_name}'
      f' --key {key}'
      f' --upload-id {upload_id}',
      shell = True,
      stdout = subprocess.PIPE
    )

    return proc.stdout


  @staticmethod
  def del_all_multipart_uploads(
    profile_name: str,
    bucket_name: str
  ):
    """del_all_multipart_uploads.

    :profile_name str:
    :bucket_name str: 

    :return: Quantity of multipart uploads deleted.
    :rtype: int
    """

    multipart_uploads = S3.get_multipart_uploads(
      profile_name,
      bucket_name
    )
    if 'Uploads' not in multipart_uploads:
      print('There are no multipart uploads.')
      return 0
    
    multipart_upload_qty = len(multipart_uploads['Uploads'])
    print(f'There are {multipart_upload_qty} multipart uploads.')


    uploads_deleted_qty = 0
    for upload in multipart_uploads['Uploads']:

      key = upload['Key']
      upload_id = upload['UploadId']

      res = S3.del_multipart_upload(
        profile_name,
        bucket_name,
        key,
        upload_id
      )

      if res == b'':
        uploads_deleted_qty += 1
        sys.stdout.write(f'Deleted {uploads_deleted_qty} multipart uploads.\r')
    
    print(f'Deleted {uploads_deleted_qty} multipart uploads.')


    return uploads_deleted_qty

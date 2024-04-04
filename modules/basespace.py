
# App registration: https://developer.basespace.illumina.com/apps/16649633/details

import json

import requests

import config
import modules as m


CLIENT_ID = 'ab0fcf015af5479699bd654cebf81c24'
CLIENT_SECRET = 'c1690364ef514771bed6764858a6f2cc'
ACCESS_TOKEN = 'c7321deb1e7b40e8981a7c49065912c9'
BASE_URL = 'https://api.basespace.illumina.com/v2'


class BaseSpace():
    

  def __init__(self,
    clients = None
  ):
    pass
    
  @staticmethod
  def get_datasets(
    id = None
  ):
    '''
    Gets all datasets or the one with the specified ID.
    '''
    str_url = str(
      f'{BASE_URL}'
       '/datasets'
    )
    if id:
      str_url += f'/{id}'
    
    str_url += f'?access_token={ACCESS_TOKEN}'
    
    json_content = json.loads(
      requests.get(
        str_url
      ).content.decode('utf-8')
    )
    
    m.Debugger.write_to_file(json_content)
    return json_content


  @staticmethod
  def update_dataset(
    id = None,
    name = None
  ):
    '''
    data: {
      'name': <name>
    }
    '''
    
    json_content = json.loads(
      requests.post(
        url = str(
          f'{BASE_URL}'
          '/datasets'
          f'/{id}'
          f'?access_token={ACCESS_TOKEN}'
        ),
        json = {
          'Name': name
        }
      ).content.decode('utf-8')
    )
    
    m.Debugger.write_to_file(json_content)
    return json_content

  
  @staticmethod
  def get_biosamples(
    biosamplename = None,
    offset = None
  ):
    str_url = str(
      f'{BASE_URL}'
        '/biosamples'
      f'?access_token={ACCESS_TOKEN}'
    )
    if biosamplename:
      str_url += f'&biosamplename={biosamplename}'
    if offset:
      str_url += f'&offset={offset}'

    json_content = json.loads(
      requests.get(
        str_url
      ).content.decode('utf-8')
    )

    m.Debugger.write_to_file(json_content)
    return json_content
  
  
  @staticmethod
  def get_users():

    str_url = str(
      f'{BASE_URL}'
        '/users/current'
      f'?access_token={ACCESS_TOKEN}'
    )

    json_content = json.loads(
      requests.get(
        str_url
      ).content.decode('utf-8')
    )

    m.Debugger.write_to_file(json_content)
    return json_content
#!/usr/bin/python

import os

import boto3

import modules as m


debug = m.Debugger.devtools_debug

cmd = m.Command(__file__)
data_dpath = cmd.data_dpath


def test(

):
  
  # Create an EC2 client
  # ec2_client = boto3.client('ec2')
  aws = m.Aws(
    clients = ['ec2']
  )

  launch_specs = {
    
  }
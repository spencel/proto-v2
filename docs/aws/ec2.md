
Example ec2 describe_instances() response:
```
{
  "Reservations": [
    {
      "Groups": [],
      "Instances": [
        {
          "AmiLaunchIndex": 0,
          "ImageId": "ami-021af1602e4e3a5fe",
          "InstanceId": "i-01e41437371249107",
          "InstanceType": "i3.4xlarge",
          "KeyName": "worker",
          "LaunchTime": "2024-04-25T02:57:56+00:00",
          "Monitoring": {
            "State": "disabled"
          },
          "Placement": {
            "AvailabilityZone": "ap-southeast-1c",
            "GroupName": "",
            "Tenancy": "default"
          },
          "PrivateDnsName": "",
          "ProductCodes": [],
          "PublicDnsName": "",
          "State": {
            "Code": 48,
            "Name": "terminated"
          },
          "StateTransitionReason": "User initiated (2024-04-25 03:01:59 GMT)",
          "Architecture": "x86_64",
          "BlockDeviceMappings": [],
          "ClientToken": "d7112f266756fb5c832d32d5ea899cd128184acf6c921c099db916880ee1",
          "EbsOptimized": false,
          "EnaSupport": true,
          "Hypervisor": "xen",
          "InstanceLifecycle": "spot",
          "NetworkInterfaces": [],
          "RootDeviceName": "/dev/xvda",
          "RootDeviceType": "ebs",
          "SecurityGroups": [],
          "SpotInstanceRequestId": "sir-kcmes6wq",
          "StateReason": {
            "Code": "Client.UserInitiatedShutdown",
            "Message": "Client.UserInitiatedShutdown: User initiated shutdown"
          },
          "Tags": [
            {
              "Key": "Owner",
              "Value": "spencer.lank@micronbrane.com"
            },
            {
              "Key": "Pipeline",
              "Value": "SR-PaRTI-Seq"
            },
            {
              "Key": "Type",
              "Value": "Worker"
            },
            {
              "Key": "Host",
              "Value": "micronbrane.basepairtech.com"
            },
            {
              "Key": "Name",
              "Value": "prod-878"
            },
            {
              "Key": "Env",
              "Value": "prod"
            }
          ],
          "VirtualizationType": "hvm",
          "CpuOptions": {
            "CoreCount": 8,
            "ThreadsPerCore": 2
          },
          "CapacityReservationSpecification": {
            "CapacityReservationPreference": "open"
          },
          "HibernationOptions": {
            "Configured": false
          },
          "MetadataOptions": {
            "State": "pending",
            "HttpTokens": "optional",
            "HttpPutResponseHopLimit": 1,
            "HttpEndpoint": "enabled",
            "HttpProtocolIpv6": "disabled",
            "InstanceMetadataTags": "disabled"
          },
          "EnclaveOptions": {
            "Enabled": false
          },
          "PlatformDetails": "Linux/UNIX",
          "UsageOperation": "RunInstances",
          "UsageOperationUpdateTime": "2024-04-25T02:57:56+00:00",
          "MaintenanceOptions": {
            "AutoRecovery": "default"
          },
          "CurrentInstanceBootMode": "legacy-bios"
        }
      ],
      "OwnerId": "255446977961",
      "RequesterId": "060752065858",
      "ReservationId": "r-0bcf40bf5678af3d3"
    }
  ],
  "ResponseMetadata": {
    "RequestId": "d1ef48af-957c-4b14-aea3-79a106404b23",
    "HTTPStatusCode": 200,
    "HTTPHeaders": {
      "x-amzn-requestid": "d1ef48af-957c-4b14-aea3-79a106404b23",
      "cache-control": "no-cache, no-store",
      "strict-transport-security": "max-age=31536000; includeSubDomains",
      "vary": "accept-encoding",
      "content-type": "text/xml;charset=UTF-8",
      "content-length": "5298",
      "date": "Thu, 25 Apr 2024 03:11:39 GMT",
      "server": "AmazonEC2"
    },
    "RetryAttempts": 0
  }
}
```
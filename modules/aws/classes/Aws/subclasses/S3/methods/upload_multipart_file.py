#!/bin/python

import argparse
import json
import logging
import os
import subprocess
import time
import threading


"""
Note:

A multipart upload can be canceled, aborted, or deleted using:
https://awscli.amazonaws.com/v2/documentation/api/latest/reference/s3api/abort-multipart-upload.html
Example command:
	aws s3api abort-multipart-upload 
	--bucket blast-db-xz
	--key multipart/01
	--upload-id qxzgwfmH_RZeDMBazr0T5sB6HB0wlt7qZYBgdZekTx1nYWRXZeGIV0ffafWHfupZmwAITVaNsbntab30.IlkyLJBwYNZMCSaOP1G4UmG5_fKH5YxGEsH.zFvnjUxsbT8sO4gg710MTh13slM17j9Ng--
	--profile mbsl

Set your profile region to the same as the S3 bucket.
	

Step 1: Create a bucket or select on an existing one.

Step 2: Enable accelerated transfer if necessary.
https://docs.aws.amazon.com/AmazonS3/latest/userguide/transfer-acceleration.html?icmpid=docs_amazons3_console

Step 3: Split the file into parts
https://repost.aws/knowledge-center/s3-multipart-upload-cli
Template command:
	split -b <part-size> -d -a <part-number-length> <input-file> <output-prefix>
Command used:
	split -b 1G -d -a 3 blast_db.tar.xz blast_db.tar.xz

	It will generate files:
	blast_db.tar.xz001
	blast_db.tar.xz002
	...

Step 4: Initilize the mulitpart upload.
https://docs.aws.amazon.com/cli/latest/reference/s3api/create-multipart-upload.html
Used this command to initialize the multipart upload:
	aws s3api create-multipart-upload
	--bucket partiseq.db.nanopore
	--key 'Step3_ref_231125.fa'
	--profile mbsl
	--no-verify-ssl
Response:
{
    "ServerSideEncryption": "AES256",
    "Bucket": "blast-db-xz",
    "Key": "blast_db.tar.xz_",
    "UploadId": "V3YpIYSuS971lBjH9h0u0r1D2byp4G1nMDl5zSO1iwuiGC3R8gDGtvbd2thuuOjFh3MI2D3uMMgtwCzZsh2BxR9c2KZa7crNUEdD0rVF_0TdEGf5THkG_x2FfewxxlL0Fhnv8hYA3vWgYQ2cTTR9jw--"
}

Step 5: Start the multipart upload.
Need to use the upload-part command for each part
https://docs.aws.amazon.com/cli/latest/reference/s3api/upload-part.html
File parts must be in the same directory as the script.
Example command:
	python aws_s3_multipart_uploader.py

Step 6: Get the list of uploaded part numbers and ETags for step 7.
https://awscli.amazonaws.com/v2/documentation/api/latest/reference/s3api/list-parts.html
	aws s3api list-parts
	--bucket blast-db-xz
	--key 'blast_db.tar.xz_'
	--upload-id "yapv.K_bT.6ZQif8HWAux4OCYwkdigV5H9ZSTKKNJ0aIsfQ.BCk31wXOE9VHX8azymM7A9118FOlgDB3hxH0hvav6lrrUf7m0_Gvr2DbO0JOvpQqsS1BEvOYKX.CF2k8o0D9WN.4tuget.PNYYZ1eg--"
	--profile mbsl

Step 7: Complete the multipart upload.
https://awscli.amazonaws.com/v2/documentation/api/latest/reference/s3api/complete-multipart-upload.html
	aws s3api complete-multipart-upload
	--multipart-upload '{"Parts": [{"PartNumber": 1,"ETag": "052ca863fee541adc3550319f2796116"},{"PartNumber": 2,"ETag": "50dd6cea7390f45b2a2df917e949d7e3"},{"PartNumber": 3,"ETag": "af536ce051f59d1401a7d10e3cf61bdd"},{"PartNumber": 4,"ETag": "2325d1b3d5e415b3226f6584226df888"},{"PartNumber": 5,"ETag": "d7871f378f33d621e5bab8f29ddbca39"},{"PartNumber": 6,"ETag": "c783ad7f8b73fc6d1f69e96f7994b144"},{"PartNumber": 7,"ETag": "993b19343bbf6bcfd841cf412632e403"},{"PartNumber": 8,"ETag": "3cafba9985d0517b258303c1f9df6d00"},{"PartNumber": 9,"ETag": "4c0d6b7efadfd31afda56051185211df"},{"PartNumber": 10,"ETag": "c494e39564cdcfedc6f8f2b7b821f4a2"},{"PartNumber": 11,"ETag": "e0e82dd9dc393c79de6638865a784b61"},{"PartNumber": 12,"ETag": "e8ff1239892584a1b4d0b1ee7bb4f693"},{"PartNumber": 13,"ETag": "e0dbc932eafb16b242d3f3914c3853e8"},{"PartNumber": 14,"ETag": "593b95976ffeb12a091349d2f6b5b140"},{"PartNumber": 15,"ETag": "9320d05807622a5b70a81aad90493b78"}]}'
	--bucket blast-database
	--key 'Step1_ref_230828.fa.bwt'
	--upload-id "V3YpIYSuS971lBjH9h0u0r1D2byp4G1nMDl5zSO1iwuiGC3R8gDGtvbd2thuuOjFh3MI2D3uMMgtwCzZsh2BxR9c2KZa7crNUEdD0rVF_0TdEGf5THkG_x2FfewxxlL0Fhnv8hYA3vWgYQ2cTTR9jw--"
	--profile mbsl
"""


# Config
PROGRESS_FILE_EXT = '.progress'
STATUS_NOT_UPLOADED = 'not uploaded'
STATUS_DONE = 'done'


# Functions

# Deletes file parts and progress file
def delete_file_parts(fpath):
	files_deleted_qty = 0
	dpath = os.path.dirname(fpath)
	for name in os.listdir(dpath):
		path = os.path.join(dpath, name)
		if path.startswith(fpath) and path != fpath:
			part_fpath = path
			os.remove(part_fpath)
			files_deleted_qty += 1
	
	return files_deleted_qty


def make_file_parts(
		fpath,
		part_size
	):
	logging.debug(f'fpath: {fpath}')
	fname = os.path.basename(fpath)
	logging.debug(f'fname: {fname}')

	subprocess.run(
		 'split'
		f' -b {part_size}' # Size of a part in bytes
		 ' -d'     # Use numeric suffixes starting at 0, not alphabetic
		# ' -a 3'   # Generate suffixes of length N (default 2)
		f' {fpath}'  # Input file
		f' {fpath}', # Output file part prefix
		shell = True
	)

	return True


def create_multi_part_upload(profile, bucket, key):
	print(f'profile: {profile}')
	print(f'bucket: {bucket}')
	print(f'key: {key}')
	proc = subprocess.run(
		 'aws s3api create-multipart-upload'
		f' --profile "{profile}"'
		f' --bucket "{bucket}"'
		f' --key "{key}"'
		 ' --no-verify-ssl',
		shell = True,
		capture_output = True,
		text = True
	)
	return json.loads(proc.stdout)


def create_progress_file(progress_fpath, aws_profile, aws_s3_bucket, aws_s3_file_key):
	print(f'aws_profile: {aws_profile}')
	fpath = progress_fpath.replace(PROGRESS_FILE_EXT, '')
	logging.debug(f'fpath: {fpath}')
	dpath = os.path.dirname(fpath)
	logging.debug(f'dpath: {dpath}')

	progress_json = {
		'parts': dict()
	}

	for name in os.listdir(dpath):
		path = os.path.join(dpath, name)
		logging.debug(f'path: {path}')
		if path.startswith(fpath) and path != fpath:
			part_fpath = path
			index = part_fpath.replace(fpath, '')
			logging.debug(f'index: {index}')
			progress_json['parts'][index] = {
				'fpath': part_fpath,
				'status': STATUS_NOT_UPLOADED
			}
	
	upload_meta = create_multi_part_upload(aws_profile, aws_s3_bucket, aws_s3_file_key)

	for key, value in upload_meta.items():
		progress_json[key] = value

	
	with open(progress_fpath, 'w') as f:
		f.write(json.dumps(progress_json, indent=2))
	
	return progress_json


def get_progess(fpath, aws_profile, aws_s3_bucket, aws_s3_file_key):
	logging.debug(f'fpath: {fpath}')

	progress_fpath = fpath + PROGRESS_FILE_EXT
	logging.debug(f'progress_fpath: {progress_fpath}')

	if not os.path.isfile(progress_fpath):
		print("Initializing upload...")
		create_progress_file(progress_fpath, aws_profile, aws_s3_bucket, aws_s3_file_key)
	
	else:
		print("Resuming upload...")
	
	return progress_fpath


def read_process_output(process):
	
	# Function to continuously read output of the subprocess
	while True:
			output = process.stdout.readline()
			if output == b'' and process.poll() is not None:
					break
			if output:
					print(output.decode().strip())

def run_subprocess(command):
	
	# Function to run the subprocess and display real-time status
	process = subprocess.Popen(
		command,
		stdout = subprocess.PIPE,
		stderr = subprocess.STDOUT,
		text = True,
		shell = True
	)

	output_thread = threading.Thread(
		target = read_process_output,
		args = (process,),
		daemon = True
	)

	output_thread.start()
	process.wait()
	output_thread.join()


def upload_parts(progress_fpath, aws_profile):

	upload_times = [] # seconds
	avg_upload_time = 0 # seconds

	progress_json = json.load(open(progress_fpath))

	bucket = progress_json['Bucket']
	key = progress_json['Key']
	upload_id = progress_json['UploadId']

	for part_number in progress_json['parts']:
		this_file = progress_json['parts'][part_number]
		s3_part_number = int(part_number) + 1
		part_fpath = this_file["fpath"]
		status = this_file["status"]
		if status == STATUS_NOT_UPLOADED:
			print(f"Starting upload of {part_fpath}...")

			start_time = time.time()
			process = subprocess.run(
				 "aws s3api upload-part"
				f" --bucket {bucket}"
				f" --key '{key}'"
				f" --part-number {s3_part_number}"
				f" --body {part_fpath}"
				f" --upload-id \"{upload_id}\""
				f" --profile {aws_profile}",
				shell = True,
				stdout = subprocess.PIPE,
				stderr = subprocess.PIPE,
				text = True,
				check=True
			)
			elapsed_time = time.time() - start_time
			print(f"Upload time: {elapsed_time:.2f} seconds.")

			process_output = process.stdout
			process_error = process.stderr
			print("Server response:")
			print("stdout:")
			print(process_output + "\n")
			print("stderr:")
			print(process_error + "\n")

			progress_json['parts'][part_number]["status"] = STATUS_DONE
			with open(progress_fpath, 'w') as f:
				json.dump(progress_json, f, indent=2)
			
			upload_times.append(elapsed_time)
			avg_upload_time = sum(upload_times)/len(upload_times)
	
	return True


def get_uploaded_part_list(progress_fpath, aws_profile):

	progress_json = json.load(open(progress_fpath))

	proc = subprocess.run(
		 'aws s3api list-parts'
		f' --profile {aws_profile}'
		f' --bucket {progress_json["Bucket"]}'
		f' --key {progress_json["Key"]}'
		f' --upload-id {progress_json["UploadId"]}',
		shell = True,
		capture_output = True,
		text = True
	)

	part_list = json.loads(proc.stdout)
	logging.debug(part_list)

	return part_list


def complete_upload(progress_fpath, part_list, aws_profile):

	progress_json = json.load(open(progress_fpath))
	part_list = {
		'Parts': part_list['Parts']
	}

	for i, part in enumerate(part_list['Parts']):
			del part_list['Parts'][i]['LastModified']
			del part_list['Parts'][i]['Size']
			part_list['Parts'][i]['ETag'] = part_list['Parts'][i]['ETag'].replace('"', '')
	logging.debug(json.dumps(part_list, indent=2))

	part_list = json.dumps(part_list)

	key = progress_json['Key']
	upload_id = progress_json['UploadId']
	
	cmd_str = str(
		 'aws s3api complete-multipart-upload'
		f' --profile {aws_profile}'
		f' --multipart-upload \'{part_list}\''
		f' --bucket {progress_json["Bucket"]}'
		f' --key "{key}"'
		f' --upload-id "{upload_id}"'
	)
	logging.debug(cmd_str)
	proc = subprocess.run(
		cmd_str,
		shell = True,
		capture_output = True,
		text = True
	)

	res = None
	try:
		res = json.loads(proc.stdout)
	except ValueError as e:
		logging.warning(e)
		res = proc.stdout

	logging.debug(res)

	return res

@staticmethod
def upload_multipart_file(
	aws_profile: str,
	fpath: str,
	size: str,
	bucket_name: str
):
	print(f'aws_profile: {aws_profile}')
	progress_fpath = fpath + '.progress'
	aws_s3_file_key = os.path.basename(fpath)
	make_file_parts(fpath, size)
	get_progess(fpath, aws_profile, bucket_name, aws_s3_file_key)
	upload_parts(progress_fpath, aws_profile)
	part_list = get_uploaded_part_list(progress_fpath, aws_profile)
	complete_upload(progress_fpath, part_list, aws_profile)
	delete_file_parts(fpath)


if __name__ == "__main__":

	parser = argparse.ArgumentParser(
		epilog = str(
			'The SIZE argument is an integer and optional unit (example: 10K is 10*1024).\n'
			'Units are K,M,G,T,P,E,Z,Y (powers of 1024) or KB,MB,... (powers of 1000).\n'
			'Binary prefixes can be used, too: KiB=K, MiB=M, and so on.'
		),
		formatter_class = argparse.RawTextHelpFormatter
	)

	# File to upload to AWS S3 Bucket
	parser.add_argument(
		'--file',
		'-f',
		type = str,
		help = 'Filepath to upload to AWS S3 Bucket'
	)

	# Size of each file part
	parser.add_argument(
		'--size',
		'-s',
		type = str,
		help = 'Put SIZE bytes per output file'
	)
	
	# Bucket name
	parser.add_argument(
		'--bucket',
		'-b',
		type = str,
		help = 'Bucket name'
	)

	parser.parse_args()
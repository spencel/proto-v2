
import datetime
import json
import logging

import config as cf


def _serialize(obj):
	if isinstance(obj, datetime.datetime):
		return obj.isoformat()


class Log():

	default_log_fpath = cf.logging.default_log_fpath
	# Reset log file
	with open(default_log_fpath, 'w') as f:
		pass

	def __init__(self,
		module: str,
		reset_log: bool = True
	):
		# Create logger
		logger = logging.getLogger(__name__)
		# Set logging level
		logger.setLevel(logging.DEBUG)
		logger.propagate = False
		# Create log file handler
		file_handler = logging.FileHandler(self.default_log_fpath)
		# Set logging level of log file handler
		file_handler.setLevel(logging.DEBUG)
		# Create log formatter
		formatter = logging.Formatter('%(levelname)s|%(asctime)s|%(message)s')
		file_handler.setFormatter(formatter)
		logger.addHandler(file_handler)
		self.module = module
		self.logger = logger


	def write(self,
		level: str|None = None,
		name: str|None = None,
		identifier: str|None = None,
		value: any = None,
		is_print_json: bool = True
	):
		
		out_line_gen = f'{self.module}'

		if name:
			out_line_gen += f'.{name}:'
		else:
			out_line_gen += f':'

		if identifier:
			out_line_gen += f'{identifier}:'
		
		if isinstance(value, dict) and is_print_json:
			str_json = json.dumps(
				value,
				indent = 2,
				default = _serialize
			)
			out_line_gen += f'\n{str_json}'
		
		else:
			out_line_gen += f' {value}'
		
		match level:
			case 'debug': self.logger.debug(out_line_gen)
			case 'info': self.logger.info(out_line_gen)
			case 'error': self.logger.error(out_line_gen)


	def debug(self, *args):
		# print(f'kwargs: {kwargs}')
		self.write(
			'debug', # level
			*args
			# **kwargs
		)
	def info(self, *args, **kwargs):
		self.write(
			'info', # level
			*args, **kwargs
		)
	def error(self, *args, **kwargs):
		self.write(
			'error', # level
			*args, **kwargs
		)
	
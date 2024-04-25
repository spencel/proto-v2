
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
		logger = logging.getLogger(__name__)
		logger.setLevel(logging.DEBUG)
		file_handler = logging.FileHandler(self.default_log_fpath)
		file_handler.setLevel(logging.DEBUG)
		formatter = logging.Formatter('%(levelname)s|%(asctime)s|%(message)s')
		file_handler.setFormatter(formatter)
		logger.addHandler(file_handler)
		self.module = module
		self.logger = logger


	def write(self,
		name: str|None = None,
		identifier: str|None = None,
		value: any = None,
		is_print_json: bool = True
	):
		print(f'value: {value}')		
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
			out_line_gen += str(value)
		
		self.logger.debug(out_line_gen)
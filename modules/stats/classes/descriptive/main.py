
import json


class Descriptive():
	"""_summary_
		N - population
		∑ - [population] sum 
		μ (mu) - [population] mean
		Mo - mode[s] - there may be more than 1
		fMo - the frequency of the mode, there is only 1
	"""
	SAMPLE_VARIANCE_TYPE = 'SAMPLE'
	POPULATION_VARIANCE_TYPE = 'POPULATION'
	SAMPLE_STANDARD_DEVIATION_TYPE = SAMPLE_VARIANCE_TYPE
	POPULATION_STANDARD_DEVIATION_TYPE = POPULATION_VARIANCE_TYPE
	SAMPLE_STANDARD_ERROR_TYPE = SAMPLE_STANDARD_DEVIATION_TYPE
	POPULATION_STANDARD_ERROR_TYPE = POPULATION_STANDARD_DEVIATION_TYPE

	
	def __init__(self, data):
		"""_summary_

		Args:
				data (number): a list of numbers
		"""
		self.data = data

	
	def population(self):
		# Alias population
		N = len(self.data)
		return N


	@classmethod
	def get_sum(cls, data):
		# Alias ∑ - [population] sum 
		sum = 0
		for value in data:
			sum += value
		return sum

	def sum(self):
		return self.get_sum(self.data)
	


	def mean(self):
		# μ = Sum / Population
		sum = self.sum()
		N = self.population()
		# Alias μ (mu) - [population] mean
		mean = sum / N
		return mean
	

	def median(self) -> int|float:
		"""If the number of observations (n) is odd, the median is the middle value.
		  You can find it by counting (n + 1) / 2 positions from the beginning of 
			the ordered list.
			If the number of observations is even, the median is the average of the 
			two middle values. You can find them by locating n / 2 and (n / 2) + 1 
			positions from the beginning of the ordered list.
		"""
		# Sort the data in either ascending or descending order
		self.data.sort()
		# Get the population
		N = self.population()
		# Get the index of the middle value
		
		if N % 2 != 0:
			idx = int(N / 2)
			median = self.data[idx]

		else:
			idx_1 = N // 2 - 1
			first_val = self.data[idx_1]
			idx_2 = idx_1 + 1
			second_val = self.data[idx_2]
			median = (first_val + second_val) / 2
	
		return median


	def modes(self):
		"""The values that occur most frequently in a dataset.
		"""
		frequencies = {}
		for value in self.data:
			if value not in frequencies:
				frequencies[value] = 1
			else:
				frequencies[value] += 1

		# There is more than 1 highest frequency, then there are multiple modes
		# Alias Mo - mode
		modes = []
		prev_highest_frequency = 0
		for value, frequency in frequencies.items():
			curr_frequency = frequency
			if curr_frequency > prev_highest_frequency:
				modes = [value]
				prev_highest_frequency = curr_frequency
			elif curr_frequency == prev_highest_frequency:
				modes.append(value)

		# The frequency of the mode, there is only 1
		# Alias fMo
		frequency = prev_highest_frequency

		return modes, frequency
	

	def variance(self, type: str):
		# Get the mean
		mean = self.mean()
		# Get the squared differences
		# https://www.notion.so/squared-differences-the-52b72e0ee85b4d829ced3e8da081fc58
		# (x_i - mean)^2
		squared_diffs = []
		for value in self.data:
			squared_diff = (value - mean)**2
			squared_diffs.append(squared_diff)
		sum_of_squared_diffs = self.get_sum(squared_diffs)
		population = self.population()
		variance: float
		if type == self.SAMPLE_VARIANCE_TYPE:
			variance = sum_of_squared_diffs / (population - 1)
		elif type == self.POPULATION_VARIANCE_TYPE:
			variance = sum_of_squared_diffs / population
		return variance


	def standard_deviation(self, type: str):
		variance: float
		if type == self.SAMPLE_STANDARD_DEVIATION_TYPE:
			variance = self.variance(self.SAMPLE_VARIANCE_TYPE)
		elif type == self.POPULATION_STANDARD_DEVIATION_TYPE:
			variance = self.variance(self.POPULATION_VARIANCE_TYPE)
		
		# Get standard variation by taking the square root of the variance
		standard_deviation = variance ** 0.5

		return standard_deviation
	
	
	def standard_error(self, type: str):
		standard_deviation: float
		if type == self.SAMPLE_STANDARD_ERROR_TYPE:
			standard_deviation = self.standard_deviation(self.SAMPLE_STANDARD_DEVIATION_TYPE)
		elif type == self.POPULATION_STANDARD_ERROR_TYPE:
			standard_deviation = self.standard_deviation(self.POPULATION_STANDARD_DEVIATION_TYPE)
		
		population = self.population()
		standard_error = standard_deviation / (population ** 0.5)

		return standard_error
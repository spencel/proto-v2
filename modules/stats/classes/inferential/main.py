
import scipy
import scipy.stats

from modules.stats.classes.descriptive import Descriptive


class Inferential():

	ONE_SIDED_TEST_TYPE = 'ONE_SIDED_TEST'
	TWO_SIDED_TEST_TYPE = 'TWO_SIDED_TEST'


	def __init__(self, data):
		self.data = data

		
	def test_statistic(self, null_hypothesis):
		# t = ("sample mean" - "hypothesized mean") /
		#		("sample standard deviation" / "sample size"^(1/2))
		hypothesized_mean = null_hypothesis
		descriptive_stats = Descriptive(self.data)
		# The sample mean is the mean of this population subset (not n - 1)
		sample_mean = descriptive_stats.mean()
		sample_standard_deviation = descriptive_stats.standard_deviation(
			type = Descriptive.SAMPLE_STANDARD_DEVIATION_TYPE
		)
		# The sample population is the population of this subset (not n - 1)
		sample_population = descriptive_stats.population()
		test_statistic = (sample_mean  - hypothesized_mean) / \
			(sample_standard_deviation / (sample_population ** 0.5))
		
		return test_statistic
	
	def p_value(self, sides, test_statistic=None):
		"""https://www.statology.org/how-to-calculate-a-p-value-from-a-t-test-by-hand/

		Args:
				test_statistic (_type_, optional): _description_. Defaults to None.
		"""
		p_value: float
		if sides == self.ONE_SIDED_TEST_TYPE:
			if test_statistic > 0:
				p_value = 1 - scipy.stats.norm.cdf(test_statistic)
			else:
				p_value = scipy.stats.norm.cdf(test_statistic)
		elif sides == self.TWO_SIDED_TEST_TYPE:
			p_value = 2 * scipy.stats.norm.cdf(-abs(test_statistic))
		
		return p_value
	
	def t_test(self, hypothesized_mean):
		descriptive_stats = Descriptive(self.data)
		sample_population = descriptive_stats.population()
		sample_mean = descriptive_stats.mean()
		sample_standard_deviation = descriptive_stats.standard_deviation(
			type = Descriptive.SAMPLE_STANDARD_DEVIATION_TYPE
		)
		test_statistic = self.test_statistic()
		# Alias degrees of freedom
		df = sample_population - 1

		return test_statistic, df
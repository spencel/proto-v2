
import math


class Combinatorics():

	# n ≥ k ≥ 0
	@staticmethod
	def get_binomial_coefficient(n, k):
		if k > n:
			return None
		if k < 0:
			return None
		return int(math.factorial(n) \
			/ (math.factorial(k) * math.factorial(n - k)))

	@classmethod
	def combination_formula(cls, n, k):
		return cls.get_binomial_coefficient(n, k)

import json

from modules import stats as s


def stats():

	print('Descriptive Stats')

	# Test median                              v
	descriptive_stats_a = s.Descriptive([0, 1, 2, 3, 4])
	median_a = descriptive_stats_a.median()
	print(f' median_a: {median_a}')
	assert median_a == 2
  
	#                                          v  v
	descriptive_stats_b = s.Descriptive([0, 1, 2, 3, 4, 5])
	median_b = descriptive_stats_b.median()
	print(f' median_b: {median_b}')
	assert median_b == 2.5

	
	data = [0, 0, 0, 0, 1, 1, 2, 3, 7, 8, 9, 10, 10, 10, 4, 4, 5, 6, 6, 6, 6]
	descriptive_stats = s.Descriptive(data)

	population = descriptive_stats.population()
	print(f' population: {population}')

	sum = descriptive_stats.sum()
	print(f' sum: {sum}')

	# Test mean
	mean = descriptive_stats.mean()
	print(f' mean: {mean}')

	# Test mode
	modes, frequency = descriptive_stats.modes()
	print(f' modes: {json.dumps(modes, indent=2)}')
	print(f'  frequency: {frequency}')

	# Test sample variance
	sample_variance = descriptive_stats.variance(type='SAMPLE')
	print(f' sample_variance: {sample_variance}')
	assert sample_variance == 12.633333333333336
	# Test population variance
	population_variance = descriptive_stats.variance(type='POPULATION')
	print(f' population_variance: {population_variance}')
	assert population_variance == 12.031746031746033

	# Test sample standard deviation
	sample_standard_deviation = descriptive_stats.standard_deviation(type='SAMPLE')
	print(f' sample_standard_deviation: {sample_standard_deviation}')
	assert sample_standard_deviation == 3.554340070017687
	# Test population standard deviation
	population_standard_deviation = descriptive_stats.standard_deviation(type='POPULATION')
	print(f' population_standard_deviation: {population_standard_deviation}')
	assert population_standard_deviation == 3.468680733614155

	# Test sample standard error
	sample_standard_error = descriptive_stats.standard_error(type='SAMPLE')
	print(f' sample_standard_error: {sample_standard_error}')
	assert sample_standard_error == 0.7756205912605091
	# Test population standard error
	population_standard_error = descriptive_stats.standard_error(type='POPULATION')
	print(f' population_standard_error: {population_standard_error}')
	assert population_standard_error == 0.7569281915915153


	
	inferential_stats = s.Inferential(data)

	# Test test statistic
	print('test-statistic: hypothesized_mean|test_statistic')
	hypothesized_means = list(range(0, len(data)))
	test_statistics = []
	for hypothesized_mean in hypothesized_means:
		test_statistics.append(inferential_stats.test_statistic(
			null_hypothesis=hypothesized_mean
		))
	for i in range(len(hypothesized_means)):
		print(f' {hypothesized_means[i]}|{test_statistics[i]}')

	# Test one sided p value
	print('test-statistic: hypothesized_mean|test_statistic|one_sided_p_value')
	one_sided_p_values = []
	for test_statistic in test_statistics:
		one_sided_p_values.append(inferential_stats.p_value(
			sides = 'ONE_SIDED_TEST',
			test_statistic = test_statistic
		))
	for i in range(len(hypothesized_means)):
		print(f' {hypothesized_means[i]}|{test_statistics[i]}|{one_sided_p_values[i]}')
	# Test two sided p value
	print('test-statistic: hypothesized_mean|test_statistic|one_sided_p_value|two_sided_p_value')
	two_sided_p_values = []
	for test_statistic in test_statistics:
		two_sided_p_values.append(inferential_stats.p_value(
			sides = 'TWO_SIDED_TEST',
			test_statistic = test_statistic
		))
	for i in range(len(hypothesized_means)):
		print(f' {hypothesized_means[i]}|{test_statistics[i]}|{one_sided_p_values[i]}|{two_sided_p_values[i]}')
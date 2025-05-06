
from typing import Literal

import scipy


# t-test
T_TEST_INDEPENDENT = 'independent'
def t_test(
	series1, 
	series2,
	type: Literal[T_TEST_INDEPENDENT]
):
	match type:
		case T_TEST_INDEPENDENT:
			t_stat, p_value = scipy.ttest_ind(series1, series2)
	return t_stat, p_value
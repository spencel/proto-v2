
import datetime
import json
import math
import os
import pathlib
import shutil
import traceback


# Config



def get_scope_pref(parent, scopeName):
    return '%s.%s' % (parent, scopeName)

# File system paths #########################################################################################################

# Returns True if directory didn't exist, so was created
def create_dpath_if_not_exist(dpath):
    path = pathlib.Path(dpath)
    if path.is_dir():
        return False
    else:
        path.mkdir(parents=True, exist_ok=True)
        return True

# Removes a directory and its contents
def delete_directory(dpath):
    shutil.rmtree(dpath)

def split_path(path):
    path = os.path.normpath(path)
    return path.split(os.sep)


# Time ######################################################################################################################

def get_previous_date(previous_day_qty):
    return datetime.date.today() - datetime.timedelta(days=previous_day_qty - 1)


# Performance Testing #######################################################################################################

class StopWatch():

    def __init__(self):
        self.time_marks = list()

    # Start also resets
    def start(self, print_start=False):
        start_time = datetime.datetime.now()
        self.time_marks = [start_time]
        if print_start:
            print(f"Start time: {str(start_time)}")

    # Split
    def split(self, print_split=False):
        stop_time = datetime.datetime.now()
        self.time_marks.append(stop_time)
        if print_split:
            prev_time = self.time_marks[-2]
            elapsed_time = stop_time - prev_time
            print(f"Split time: {str(elapsed_time)}")

    # Reset simply starts the stopwatch over again
    def reset(self):
        self.start()

# String date to datetime object
def str_to_datetime(str_date, delimeter):

    arr_date = str_date.split("-")
    year = 1
    if len(arr_date) > 0:
        year = int(arr_date[0])
    month = 1
    if len(arr_date) > 1:
        month = int(arr_date[1])
    day = 1
    if len(arr_date) > 2:
        day = int(arr_date[2])

    return datetime.datetime(year, month, day)



# Print an object/variable, and its identifier and scope
# Patterns:
## Get the scope prefix
#### scopePref # = '%s.%s' % (__name__, <function-name>.__name__)
## Or
#### scopePref # = get_scope_pref(parent, scopeName)
## Print a variable name and its value
#### DebugLog # (debug_print(scopePref, '<variable-name>', <variable>))
## Print a variable name and its value as JSON
#### DebugLog # (debug_print(scopePref, '<variable-name>', <variable>, True))
## Print message
#### DebugLog # (debug_print(scopePref, message=<str>))
# To remove debug statements, Ctrl + Shift + R to open replace all window and use these regular expressions:
## Replace "(|.+)scopePref =.+\n" with ""
## Replace "(|.+)DebugLog\(debug_print\(.+\n" with ""
def debug_print(scopePref, name=None, value=None, asJson=False, message=None):
    printStr = '%s:' % (scopePref)

    # If it's just a message, print it then exit function
    if message:
        printStr = '%s %s' % (printStr, message)
        print(printStr)
        return printStr

    # If it's not just a message, then assume we are printing an object and its value
    if name:
        printStr = '%s %s' % (printStr, name)
        if asJson:
            str_value = json_to_str(value)
            if str_value:
                printStr = '%s: %s' % (printStr, str_value)
            else:
                # Line of code below doesn't work yet?
                # printStr = '%s: [not json serializable] %s' % (printStr, json.dumps(make_jsonable(value), indent=2))
                printStr = '%s: [not json serializable] %s' % (printStr, value)
        else:
            printStr = '%s: %s' % (printStr, value)

    print(printStr)
    return printStr

def invert_str(str_arg):
    return "".join(reversed(str_arg))

# Get the complement nucleic acid sequence
def get_cna_seq(na_seq):
    new_na_seq = ""
    for base in na_seq:
        if base == "T":
            new_na_seq += 'a'
        elif base == 'a':
            new_na_seq += "T"
        elif base == "G":
            new_na_seq += "C"
        elif base == "C":
            new_na_seq += "G"
    return new_na_seq

# Get the reverse complement nucleic acid sequence
## Eg, of a reverse primer
def get_rev_cna_seq(na_seq):
    na_seq = invert_str(na_seq)
    return get_cna_seq(na_seq)


# Get the average quantity of days per month based on the orbital period of the the Earth
def get_avg_d_per_m():
    EARTH_ORBITAL_PERIOD = 365.256363004 # https://en.wikipedia.org/wiki/Earth
    return EARTH_ORBITAL_PERIOD / 12

# Get Timestamp for Filename Purposes
# Prevents filename collisions
# Can add additional use cases if needed
def get_fname_timestamp():
    try:

        # Format: YYYYmmddHHMMssssss
        # full year, month, day, hour, minute, microseconds
        return datetime.datetime.now().strftime('%Y%m%d%H%M%f')

    # Catch any exception
    except Exception as e:
        print('\nUnexpected exception:')
        print(traceback.format_exc())
        raise e

def get_timestamped_filepath(dirpath, baseFilename):
    try:

        stamp = get_fname_timestamp()
        filepath = "{}{}{}.csv".format(
                dirpath,
                baseFilename,
                stamp
            ).lower().replace(' ', '-')

        return filepath

    # Catch any exception
    except Exception as e:
        print('\nUnexpected exception:')
        print(traceback.format_exc())
        raise e

def get_env_as_json(env):
    return json.dumps({**{}, **env}, indent=2)

def line_to_list(line):
    return line.strip('\n').split('\t')

def list_to_string(data, delimiter=" | ", quotes=None):
    # data can be a list or set
    scope_pref = __name__
    # logging.debug(f"{scope_pref}: data: {data}")
    # logging.debug(f"{scope_pref}: type(data): {type(data)}")
    str_list = ""
    if delimiter:
        if quotes:
            str_list = f"{quotes}{delimiter}{quotes}".join(map(str, data))
        else:
            str_list = f"{delimiter}".join(map(str, data))
    else:
        str_list = "".join(map(str, data))
    # logging.debug(f"{scope_pref}: str_list A: {str_list}")
    # Remove trailing brackets of set or list
    # if isinstance(data, set): didn't work for some reason with postgres fetch
    # logging.debug(f"{scope_pref}: str_list: {str_list}")
    if str_list[0] == "(" and str_list[-1] == ")":
        str_list = str_list.lstrip("(").rstrip(")")
    # elif isinstance(data, list): didn't work for some reason with postgres fetch
    if str_list[0] == "[" and str_list[-1] == "]":
        # logging.debug(f"{scope_pref}: isinstance(data, list): {isinstance(data, list)}")
        str_list = str_list.lstrip("[").rstrip("]")
    if quotes:
        str_list = quotes + str_list + quotes
    # logging.debug(f"{scope_pref}: str_list B: {str_list}")
    # Remove extra comma if delimiter was specified
    # str_list = str_list[:len(str_list) - len(delimiter)]

    return  str_list

def list_members(object, outfilepath=None):

    def recurse_object(object, level=0, outfile=None):

        def handle_item(item, value):
            if item:
                print(" " * level + item)
                if outfile:
                    outfile.write(f"{' '*level}{item}\n")
            if isinstance(value, (dict, list, set)) or hasattr(value, '__class__'):
                next_level = level + 2
                recurse_object(value, next_level, outfile)

        # Class handler
        if hasattr(object, "__class__") and hasattr(object, "__dict__"):
            # logging.debug(f"{scope_pref}: object: {object}")
            for property, value in vars(object).items():
                handle_item(property, value)

        # Dictionary handler
        if isinstance(object, dict):
            for key, value in object.items():
                handle_item(key, value)

        # List/set handler
        if isinstance(object, list) or isinstance(object, set):
            print(" " * level + "array")
            if outfile:
                outfile.write(f"{' ' * level}array\n")
            for value in object:
                handle_item(None, value)

    if outfilepath:
        with open(outfilepath, 'w') as f:
            recurse_object(object, outfile=f)
    else:
        recurse_object(object)

# GISAID Data column names and index
##  1	covv_accession_id
##  2	covv_add_host_info
##  3	covv_sampling_strategy
##  4	covv_add_location
##  5	covv_assembly_method
##  6	covv_clade
##  7	covv_collection_date
##  8	covsurver_prot_mutations
##  9	gc_content
## 10	covv_gender
## 11	covv_host
## 12	is_high_coverage
## 13	is_reference
## 14	is_complete
## 15	covv_lineage
## 16	pangolin_lineages_version
## 17	covv_location
## 18	n_content
## 19	covv_outbreak
## 20	covv_passage
## 21	covv_patient_age
## 22	covv_patient_status
## 23	covsurver_existingmutlist
## 24	sequence
## 25	sequence_length
## 26	covv_seq_technology
## 27	covv_specimen
## 28	covv_subm_date
## 29	covv_type
## 30	covsurver_uniquemutlist
## 31	covv_variant
## 32	covv_virus_name

# Probability & Statistics

# Combinations and permutations
# https://www.mathplanet.com/education/pre-algebra/probability-and-statistics/combinations-and-permutations
# https://www.calculator.net/permutation-and-combination-calculator.html?cnv=10&crv=2&x=96&y=17
def get_permutation_qty(set_length, subset_length=None, all_subset_lengths=False):

    if not subset_length:
        subset_length = set_length
    # print(f"subset_length: {subset_length}")

    start_length = subset_length
    if all_subset_lengths:
        start_length = 1
    # print(f"start_length: {start_length}")

    permutation_qty = 0

    for length in range(start_length, subset_length+1):
        # print(f"length: {length}")

        permutation_qty += \
            math.factorial(set_length) \
            / math.factorial(set_length - length)

    return int(permutation_qty)

# https://www.calculator.net/permutation-and-combination-calculator.html?cnv=10&crv=2&x=96&y=17
def get_permutation_with_replacement_qty(set_length, subset_length=None, all_subset_lengths=False):

    if not subset_length:
        subset_length = set_length
    # print(f"subset_length: {subset_length}")

    start_length = subset_length
    if all_subset_lengths:
        start_length = 1
    # print(f"start_length: {start_length}")

    permutation_qty = 0

    for length in range(start_length, subset_length+1):
        # print(f"length: {length}")
        permutation_qty += math.pow(set_length, length)
        # print(f"permutation_qty: {permutation_qty}")

    return int(permutation_qty)

def get_combination_qty(set_length, subset_length=None, all_subset_lengths=False):

    if not subset_length:
        subset_length = set_length
    # print(f"subset_length: {subset_length}")

    start_length = subset_length
    if all_subset_lengths:
        start_length = 1
    # print(f"start_length: {start_length}")

    combination_qty = 0

    for length in range(start_length, subset_length + 1):
        # print(f"length: {length}")

        combination_qty += \
            math.factorial(set_length) \
            / math.factorial(length) \
            / math.factorial(set_length - length)

    return int(combination_qty)

def get_combination_with_replacement_qty(set_length, subset_length=None, all_subset_lengths=False):

    if not subset_length:
        subset_length = set_length
    # print(f"subset_length: {subset_length}")

    start_length = subset_length
    if all_subset_lengths:
        start_length = 1
    # print(f"start_length: {start_length}")

    combination_qty = 0

    for length in range(start_length, subset_length + 1):
        # print(f"length: {length}")

        combination_qty += \
            math.factorial(set_length + length - 1) \
            / math.factorial(length) \
            / math.factorial(set_length - 1)
        # print(f"combination_qty: {combination_qty}")

    return int(combination_qty)

def print_permutations_and_combinations(set_length, subset_length=None):
    print(f"perms no rep one: {get_permutation_qty(set_length, subset_length)}")
    print(f"perms no rep all: {get_permutation_qty(set_length, subset_length, True)}")
    print(f"perms replac one: {get_permutation_with_replacement_qty(set_length, subset_length)}")
    print(f"perms replac all: {get_permutation_with_replacement_qty(set_length, subset_length, True)}")
    print(f"combs no rep one: {get_combination_qty(set_length, subset_length)}")
    print(f"combs no rep all: {get_combination_qty(set_length, subset_length, True)}")
    print(f"combs replac one: {get_combination_with_replacement_qty(set_length, subset_length)}")
    print(f"combs replac all: {get_combination_with_replacement_qty(set_length, subset_length, True)}")
"""

There is a pattern to follow when listing combinations, which is a helpful visualization.

Example A
Set:           0 1
Set length:    2
Subset length: 2

Permutations:
	qty: 2
	0 1
	1 0

Permutations with replacement:
	qty: 4
	0 0
	0 1
	1 1
	1 0

Combinations:
	qty: 1
	0 1

Combinations with replacement:
	qty: 3
	0 0
	0 1
	1 1

---

Example B
Set: 0 1 2
Set Length: 3
Subset Length: 3

Permutations:
	qty: 6
	0 1 2 | 0 2 1
	1 0 2 | 1 2 0
	2 0 1 | 2 1 0 

Permutations All Subset Lengths:
    qty: +6
    same as above
    qty: +6
	0 1 | 0 2
	1 0 | 1 2
	2 0 | 2 1
	qty: +3
	1 | 2 | 3
	total qty: 15	

Permutations with replacement:
	qty: 27
    0 0 0 | 0 0 1 | 0 1 0 | 0 1 1
	        0 0 2 | 0 2 0 | 0 2 2
	                0 1 2
	                0 2 1
	1 1 1 | 1 1 0 | 1 0 1 | 1 0 0
	        1 1 2 | 1 2 1 | 1 2 2
	                1 0 2
	                1 2 0
    2 2 2 | 2 2 0 | 2 0 2 | 2 0 0
	        2 2 1 | 2 1 2 | 2 1 1
	                2 0 1
	                2 1 0
	                
Permutations with replacement All Subset Lengths:
    qty: +27
    same as above
    qty: +9
    0 0 | 0 1 | 0 2
	1 0 | 1 1 | 1 2
	2 0 | 2 1 | 2 2
    qty: +3
    1 | 2 | 3
    total qty: 39

Combinations:
	qty: 1
	0 1 2

Combinations All Subset Lengths:
	qty: +1
	0 1 2
	qty: +3
	0 1 | 0 2 | 1 2
	qty: +3
	0 | 1 | 2
	total qty: 7

Combinations with replacement:
	qty: 10
	0 0 0 | 0 0 1 | 0 0 2
	 0 1 1 | 0 1 2 | 0 2 2
	1 1 1 | 1 1 2 | 1 2 2
	2 2 2
	
Combinations with replacement All Subset Lengths:
    qty: +10
    same as above
    qty: +6
	0 0 | 0 1 |	0 2
	1 1 | 1 2
	2 2
	qty: +3
	0 | 1 | 2
	total qty: 19

---

Example C
Set: 0 1 2
Set Length: 3
Subset Length: 2

Permutations:
	qty: 6
	0 1 | 0 2
	1 0 | 1 2
	2 0 | 2 1

Permutations All Subset Lengths:
    qty: +6
    same as above
	qty: +3
	1 | 2 | 3
	total qty: 9

Permutations with replacement:
	qty: 9
	0 0 | 0 1 | 0 2
	1 0 | 1 1 | 1 2
	2 0 | 2 1 | 2 2

Combinations:
	qty: 3
	0 1
	0 2
	1 2

Combinations with replacement:
	qty: 6
	0 0 | 0 1 |	0 2
	1 1 | 1 2
	2 2

"""

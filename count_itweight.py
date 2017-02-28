"""
Author: Melanie van den Bosch
Script: Counting the number of iterations and difference count
for each genomeID from a tmpweight file .
"""

from sys import argv
import re


def get_filenames(argv):
    """ Retrieving and reading filenames from cmd line 
    """
    file_input = argv[1]
    file_output = argv[2]
    return file_input, file_output

def make_records(file_input):
    """ Makes records of each genomeID and its iterations
    """
    with open(file_input) as w_file:
        raw_record = []
        for line in w_file:
            if line.startswith("GCA_"):
				
		raw_record.append(line)		
		yield raw_record
		raw_record = []
            else:
                raw_record.append(line)
		#yield raw_record

def tmpweight_parser(file_input):
    """ Parsing the lines of differences count and genomeIDs into a dictionary
    """
    count_dict = {}
    ID_pattern = re.compile(r"(GCA_\d{9})\s*")
    diffcount_pattern = re.compile(r"differences between tables is\s*(\d+)\s*")
    for line in open(file_input):
        ID_match = ID_pattern.match(line)
        if ID_match:
	   ID = ID_match.group(1)
	   count_dict[ID] = []
	diffcount_match = diffcount_pattern.match(line)
	if diffcount_match:
	   diffcount = diffcount_match.group(1)
           count_dict[ID] += [diffcount]
    return count_dict

		
def write_dict_records(countdict):
    records = []
    dict_records = {}
    for key, value in countdict.items():
        print key
        print value
        dict_records["ID"] = key
        dict_records["diffcount"] = value
        #dict_records = {}
    	records.append(dict_records)

    return records


def write_file(file_input, file_output):
    """ Write the count of iterations and difference in iterations to a file
    """
    with open(file_output, 'w') as outp:
        count_dict = count_it_diff(file_input)
	    





if __name__ == "__main__":
    #get filenames
    infile, outfile = get_filenames(argv)
    dictcount = tmpweight_parser(infile)
    print dictcount
    recs = write_dict_records(dictcount)
    print recs
    

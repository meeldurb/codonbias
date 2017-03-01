"""
Author: Melanie van den Bosch
Script: Counting the number of iterations and difference count
for each genomeID from a tmpweight file .
"""

from sys import argv
import re
from operator import itemgetter
import csv

def get_filenames(argv):
    """ Retrieving and reading filenames from cmd line 
    
    Keyword arguments:
	argv -- list of arguments for this script
    Returns:
         file_input -- string, name of tmpweight file as argument
         file_output -- string, name of file that info needs to be written to
    """
    file_input = argv[1]
    file_output = argv[2]
    return file_input, file_output


def tmpweight_parser(file_input):
    """ Parsing the lines of differences count and genomeIDs into a dictionary
        
    Keyword arguments;
        file_input -- string, name of tmpweight file
    Returns:
        count_dict -- dictionary, contains all the genomeIDs and the difference 
        count for each iteration it has been trough for making the iterative weight 
        tables
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

		
def convert_dict_to_records(countdict):
    """ Convert the dictionary to list of dictionaries with keys:
	
    Keyword arguments:
        countdict -- large dictionary from tmpweight file
    Returns:
        records -- list of dictionaries. Contents are ID, diffcount (difference     
                    count for each iteration, iterations (number of iterations)
    """
    records = []
    dict_records = {}
    for key, value in countdict.items():
        dict_records["ID"] = key
        dict_records["diffcount"] = value
	dict_records["iterations"] = len(value)
        records.append(dict_records.copy())

    return records


def write_file(records, file_output):
    """ Write the count of iterations and difference in iterations to a file
    
    Keyword Arguments: 
        records -- list of dicts, containing the keys ID, iterations and diffcount 
                   for each dict
        file_output -- string, file the dicts need to be written to
    """
    with open(file_output, 'wb') as csvfile:
        writer = csv.writer(csvfile, delimiter=';')
        header = ("genomeID", "difference count", "iterations")
	writer.writerow(header)
        for record in records:
	    if record:
                row = (record["ID"], record["iterations"], \
				record["diffcount"])
                writer.writerow(row)
            else:
                "No file was written"
        
	    
def order_by_genomeID(count_records, reverse=True):
    """ Return a list of iteration per genomeID dictionaries ordered on genomeID
	
    Keyword Arguments: 
	count_records -- list of dictionaries, 
        reverse -- boolean, if False orderes on ascending numerical value
    Returns:
        list, sorted values on genomeID key
    """
    return sorted(count_records, key=itemgetter("ID"), reverse=reverse)




if __name__ == "__main__":
    #get filenames
    infile, outfile = get_filenames(argv)
    # parse the tmpweight file into dictionary
    countdict = tmpweight_parser(infile)
    # convert dict to list of dicts
    recs = convert_dict_to_records(countdict)
    # sort list of dicts on genomeID
    sorted_recs = order_by_genomeID(recs, reverse=False)
    #print sorted_recs
    # write recors to csv file
    write_file(sorted_recs, outfile)
    
    
    

"""
Author: Melanie van den Bosch
Script: Converting a csv input file to a fasta file .
"""

def write_to_fasta(file_input, file_output):
    """ Converting csv file to a fasta file
    """
    for line in file_input:
        if line.startswith("pf"):
            column = line.split(',')
            ID = column[0]
            CDS_dom = column[1]
            file_output.write('>%s\n%s\n' %(ID, CDS_dom))
            
            
        


with open('tmp.csv','r') as file_input:
    with open('tmp.fasta', 'w') as file_output:
        write_to_fasta(file_input, file_output)

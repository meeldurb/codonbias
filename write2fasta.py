"""
Author: Melanie van den Bosch
Script: Converting a csv input file of 
sequences from SPARQL query to a fasta file .
"""


def get_filenames(argv):
    """ Retrieving and opening file from cmd line 
    """
    file_input = argv[1]
    file_output = argv[2]

def write_fasta(file_input, file_output):
    """ Write files to fasta
    """

    with open(file_input,'r') as csv_input:
        with open(file_output, 'w') as fasta_output:
            convert_to_fasta(file_input, file_output)



def convert_to_fasta(file_input, file_output):
    """ Converting csv file to a fasta file
    """
    for line in file_input:
        if line.startswith("<"):
            column = line.split(',')
            ID = column[0]
            CDS_gene = column[1]
            file_output.write('>%s\n%s\n' %(ID, CDS_gene))
        if line.startswith("pf"):
            column = line.split(',')
            ID = column[0]
            CDS_dom = column[1]
            file_output.write('>%s\n%s\n' %(ID, CDS_dom))
                    
            
            

                    




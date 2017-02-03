"""
Author: Melanie van den Bosch
Script: Converting a csv input file of gene
sequences from SPARQL query to a fasta file .
"""

def write_to_fasta(file_input, file_output):
    """ Converting csv file to a fasta file
    """
    for line in file_input:
        if line.startswith("<"):
            column = line.split(',')
            ID = column[0]
            CDS_gene = column[1]
            file_output.write('>%s\n%s\n' %(ID, CDS_gene))
            
#def get_file(argv):
 #   """ Retrieving and opening file from cmd line 
  #  """
   # file_name = argv[1]            
        


with open('tmpcai.csv','r') as file_input:
    with open('tmpcai.fasta', 'w') as file_output:
        write_to_fasta(file_input, file_output)

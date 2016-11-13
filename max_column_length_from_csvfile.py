import math

def read_csv(filename):
    matrix = []
    csvfile = open(filename)
    for line in csvfile:
        line = line.split(',')
        print line
        column = line[0]
        column_len = len(column)
        matrix += [column_len]
    max_count = max(matrix)
    csvfile.close()
    return max_count


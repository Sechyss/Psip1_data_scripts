def file_from_list(inputlist, outfile):
    with open(outfile, 'a') as file:
        for i in inputlist:
            file.write(str(i) + '\n')


def list_from_file(inputfile):
    with open(inputfile, 'r') as file:
        lines = file.read().split('\n')
    return list(lines)


def parse_conf_file(file_path):
    with open(file_path, 'r') as f:
        data = f.readlines()
        d = {}
        for line in data:
            content = line.split()
            d[content[0]] = [bond.split('-') for bond in content[1:]]
    return d

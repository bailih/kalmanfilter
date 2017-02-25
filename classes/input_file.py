from ast import literal_eval

class InputFile:

    def __init__(self, filename, delimiter='\t'):
        head = {}

        with open(filename, 'r') as f:
            flines = [line.rstrip() for line in f]  #  use rstrip to remove newlines

        line = flines.pop(0)

        # If the header exists, parse it
        if line == '--HEAD':
            line = flines.pop(0)
            while line != '--BODY':
                try:
                    key, value = line.replace(' ', '').split('=')
                except ValueError:
                    raise ValueError('HEAD not formatted correctly. Are you missing --BODY?')
                head[key] = value
                line = flines.pop(0)

        self.head = head
        self.lines = flines
        self.delimiter = delimiter

    def adjustConfig(self, kalmanConfig, svsfConfig):
        ''' Adjust the config based on the header of the input file. '''
        order = kalmanConfig[0]

        head = self.head
        if 'cell_capacity' in head:
            kalmanConfig[5]['cnom'] = float(head['cell_capacity'])
            svsfConfig[2]['cnom'] = float(head['cell_capacity'])

        if 'c' in head:
            kalmanConfig[6] = literal_eval(head['c'])
            svsfConfig[3] = literal_eval(head['c'])

        try:
            params = literal_eval(head['o{}'.format(order)])
            for key in params:
                kalmanConfig[5][key] = params[key]
                svsfConfig[2][key] = params[key]

        finally:
            for config in kalmanConfig:
                try:
                    print(kalmanConfig[config])
                except:
                    pass
            return kalmanConfig, svsfConfig


    def __getitem__(self, key):
        ''' Returns a pre-split line, or object from header '''
        head = self.head

        if type(key) is int:
            line = self.lines[key].split(self.delimiter)
            if 'dt' in head:
                line.insert(2, head['dt'])
            return line

        elif type(key) is str:
            return head[key]

        else:
            raise TypeError('InputFile key must be either integer or string!')

    def __len__(self):
        return len(self.lines)


def _test():
    f = InputFile('Inputs/c7_9.txt')
    print(f[4])
    print(f['cnom'])
    print(f[0.2])


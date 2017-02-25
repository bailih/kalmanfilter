from sys import stdout


class ProgressIndicator:
    '''
        Static singleton to allow easy progress bar generation and tagging

        Progress bars are useful, but needing to add counting and bar generation to each
        part of code you would like to track is cumbersome and ugly.
    '''

    s_instance = None

    def __init__(self):
        self.bar = ''
        self.tag = ''

    def init_bar(self, length: 'length of the bar', end: 'what we are counting towards'):
        # We set the initial count to -1 so when we update the bar, it gets put to 0
        self.prog_count = -1
        self.length = length
        self.end    = end
        self.update_bar()

    def update_bar(self):
        self.prog_count += 1
        length = self.length
        count = self.prog_count

        blocks = count // (self.end // length)
        self.bar = '[{}{}]'.format('#' * blocks ,  '-' * (length - blocks))

        self._write()

    def set_tag(self, tag):
        self.tag = tag
        self._write()

    def _write(self):
        stdout.write('\r{} {}'.format(self.bar, self.tag))

    @staticmethod
    def get_instance():
        if ProgressIndicator.s_instance is None:
            ProgressIndicator.s_instance = ProgressIndicator()
        return ProgressIndicator.s_instance
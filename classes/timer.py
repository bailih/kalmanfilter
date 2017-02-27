from datetime import datetime
import sys

class Timer:

    s_instance = None

    def __init__(self):
        self.name_dict = {}
        self.t0 = None

    def time(self, func, args, name, iter=1):
        return_value = None
        for i in range(iter):

            t0 = datetime.now()
            return_value = func(*args)
            t1 = datetime.now()

            if name not in self.name_dict:
                self.name_dict[name] = []
            self.name_dict[name].append((t1 - t0).total_seconds())

        return return_value

    def stopwatch(self, name):
        t1 = datetime.now()
        if self.t0 is not None:
            delta_t = (t1 - self.t0).total_seconds()
            # try-catch is faster than checking with if
            try:
                self.name_dict[name].append(delta_t)
            except KeyError:
                self.name_dict[name] = [delta_t]
            self.t0 = None
        else:
            self.t0 = t1



    def print(self):
        name_dict = self.name_dict
        string = '\n'
        for name in name_dict:
            sum = 0
            for time in name_dict[name]:
                sum += time
            string += '{}: {}s\n'.format(name, sum / len(name_dict[name]))
            with open('jim_jam.csv', 'a') as f :
                f.write('{},\n'.format(sum / len(name_dict[name])))
        sys.stdout.write(string)

    @staticmethod
    def get_instance():
        if Timer.s_instance is None:
            Timer.s_instance = Timer()
        return Timer.s_instance

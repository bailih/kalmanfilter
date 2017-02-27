

class SoCTracker():

    def __init__(self, cell_count: float, cell_capacity: float, init_soc: float, balance_current: float):
        self.balance_current = balance_current

        self.cells = []
        self.cell_count = cell_count
        for i in range(int(cell_count)):
            self.cells.append(Cell(cell_capacity, 0))

        # Keep track of average SoC and minimum SoC for the cells to see which is a better estimate of battery SoC
        self.avg_soc = []
        self.min_soc = []

        self.volts   = []

    def update(self, voltage: float, current: float, duration: float, balance_info: tuple):

        balance_cell, balance_dir = balance_info

        self.volts.append(voltage)

        soc_sum = 0
        min_soc = 9999.0

        # Add current to each cell, adjusting with balance_current if the cell is balancing

        for cell in self.cells:

            if balance_cell == None or (self.cells.index(cell) + 1) != int(balance_cell):
                cell.add_current(current, duration)
            else:
                # balance_dir = 1 means discharging, so if balance_dir = 1, subtract 2 B_C to get -1 B_C
                cell.add_current(current + (self.balance_current - (2 * self.balance_current * balance_dir)), duration)

            cell_soc = cell.get_soc()

            soc_sum += cell_soc

            if cell_soc < min_soc:
                min_soc = cell_soc


        self.avg_soc.append(soc_sum / self.cell_count)
        self.min_soc.append(min_soc)


class Cell():

    def __init__(self, capacity: float, soc: float):
        self.capacity = capacity
        self.charge   = soc * capacity

        # Takes in current in milliamps and duration of current in seconds.
        # Adds to charge in cell in amp-seconds
    def add_current(self, current: float, duration: float):
        self.charge += (current / 1000) * duration

    def get_soc(self) -> float:
        return self.charge / self.capacity
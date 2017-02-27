from .kalman import Kalman
from copy import deepcopy


class MulticellKalman:
    '''
        Allows for Kalman filtering of individual cells within a battery pack, taking balancing
        into account.

        Anecdotal evidence suggests that this doesn't affect the estimation in any meaningful way,
        but this class is left here in case it proves useful in the future.
    '''

    def __init__(self, num_cells: int, balance_current: float, config: list):
        self.num_cells = num_cells
        self.balance_current = balance_current
        self.cells = [Kalman(*deepcopy(config)) for _ in range(num_cells)]

    def update(self, dt: float, voltage: float, current: float, balance_info: (int, int)):
        '''Updates each cell, adjusting current if the cell is being balanced.'''
        for i in range(self.num_cells):
            cell_current = current
            if i + 1 == balance_info[0]:
                cell_current += self.balance_current - 2 * self.balance_current * balance_info[1]
            self.cells[i].read_inputs(cell_current, voltage)
            self.cells[i].update(dt)

    def get_simulated_soc(self) -> float:
        '''Return average of all cells' estimated state of charge.'''
        return sum([cell.get_simulated_soc() for cell in self.cells]) / self.num_cells

    def get_simulated_voltage(self) -> float:
        '''Return average of all cells' estimated voltage.'''
        return sum([cell.get_simulated_voltage() for cell in self.cells]) / self.num_cells

class progressBar():
	def __init__(self, title, length=40):
		self.BAR_LENGTH = length
		self.title = title
		print('{}\t['.format(self.title) + ' ' * self.BAR_LENGTH + ']', end='')

	def update(self, val):
		# round from 0 to self.BAR_LENGTH
		bars = round(val * self.BAR_LENGTH)
		print('\r{}\t['.format(self.title) + '#' * bars + ' ' * (self.BAR_LENGTH - bars) + ']\t{0:.2f}%'.format(
			val * 100), end='')

	def close(self):
		print('')

class PulseParams():
	def __init__(self, sampleRate, rest, cov, cRate, data):

		self.restTime = rest
		self.cov = cov
		self.cRate = cRate
		self.inputData = data

		self.params = {
			'File': '',
			'Make': '',
			'Model': '',
			'COV': '',
			'RC': 0,
			'ID': '',
			'SoH': 0,
			'OCV': 0,
			'Vloss': 0,
			'IRecR1': 0,
			'OCVMinAmp': 0,
			'Current': 0,
			'HalfDischarge': 0,
			'VlossArea': 0,
			'ChargeArea': 0,
			'TotalArea': 0
		}
		self.setBatParams()

		self.rawData = []
		self.time = []
		self.voltage = []
		self.current = []
		self.sampleRate = sampleRate  # seconds
		self.extractData()

		self.initialRestData = []
		self.chargePulseData = []
		self.recoveryData = []
		self.seperatePulseSections()

		self.initOCV = self.voltage[0]
		self.criticalVoltage1 = self.chargePulseData[1][-1]
		self.criticalVoltage2 = self.recoveryData[1][1]
		self.criticalVoltage3 = self.recoveryData[1][int(self.restTime / self.sampleRate)]

	def extractData(self):
		for line in self.inputData:
			if len(self.time) == 0:
				self.time.append(0)
			else:
				self.time.append(self.time[-1] + self.sampleRate)
			self.voltage.append(float(line.split('\t')[0]))
			self.current.append(float(line.split('\t')[1]))

		self.rawData = [self.time, self.voltage, self.current]
		self.sampleRate = self.time[1]

	def setBatParams(self):

		self.params['File'] = 'IDK'
		self.params['SoH'] = 'IDK'
		self.params['ID'] = 'IDK'
		self.params['CoV'] = self.cov

	def seperatePulseSections(self):

		positiveCurrentsIndex = [index for index in range(len(self.current)) if self.current[index] > 0.045]

		endInitialRestPoint = int(positiveCurrentsIndex[0])
		endChargePulsePoint = int(60 / self.sampleRate - 300)  # this will need to be more exact
		step = 0.02
		while abs(self.voltage[endChargePulsePoint] - self.voltage[endChargePulsePoint + 1]) < step:
			endChargePulsePoint += 1
		endChargePulsePoint += 1
		endRecoveryPoint = int(endChargePulsePoint + 30 / self.sampleRate)

		self.initialRestData = [self.time[0:endInitialRestPoint],
								self.voltage[0:endInitialRestPoint],
								self.current[0:endInitialRestPoint]]

		self.chargePulseData = [self.time[endInitialRestPoint:endChargePulsePoint],
								self.voltage[endInitialRestPoint:endChargePulsePoint],
								self.current[endInitialRestPoint:endChargePulsePoint]]

		self.recoveryData = [self.time[endChargePulsePoint:endRecoveryPoint],
							 self.voltage[endChargePulsePoint:endRecoveryPoint],
							 self.current[endChargePulsePoint:endRecoveryPoint]]

	def calculateParams(self):

		self.params['Vloss'] = (self.initOCV - self.criticalVoltage3) * 1000
		self.params['IRecR1'] = (self.criticalVoltage1 - self.criticalVoltage2) * 1000
		self.params['OCVMinAmp'] = (self.initOCV - self.criticalVoltage1) * 1000

		return self.params

	def getFullParams(self):
		return self.params

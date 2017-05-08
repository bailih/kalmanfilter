from NN.rt_nonfinite import *
from NN.bsxfun import *
from math import exp


def classify(X):
	av = [0] * 6
	ixstart = 0
	Xp1 = [0] * 6

	b = [-104.0, 18.9999999999992, -407.999999999999]
	b_b = [0.0181818181818181, 0.0089686098654708, 0.0051948051948052]

	x = [0] * 6
	denominator = 0
	ix = 0
	dv0 = [0] * 6
	c_b = [0] * 6
	a = [1.0752503836391998, -0.86985261212243781,
		 0.0772419454957475, 0.5008906959897762, 0.85408498580274461,
		 1.9326654584254981]

	b_a = [-0.66706991414925121, 0.95080557570418722,
		   0.62163446263342426, 0.63445653420490788, -2.5373942800698028,
		   -4.6183499425093935, -4.3621357450446663, -0.43985387652905145,
		   -3.3014313273015379, -2.005960341272119, 2.9694398632763086,
		   -2.4597280392134646, -1.0309826845538694, 0.95250993761197045,
		   1.0259662333784862, -1.152803474248177, 4.5007910718871935,
		   -1.4008346663398641]

	c_a = [-0.43778760579936565, 1.568023028565501,
		   0.029337439033484181, 0.45795230450593971, 0.32692759236823915,
		   0.02636193992097963]

	numerator = [0] * 4
	d_a = [0.8688706560174112, 4.1275974361448649,
		   -1.2847733400951113, 0.87956649295336131, 1.047118979641422,
		   0.48244621260171711, -0.0011794635735095816, -0.6508117600670168,
		   -0.30416554274614221, 0.8970519760552188, -0.4473707529180837,
		   0.35421556216573669, -0.8823936620502234, 1.4616904273281079,
		   2.1985670479066148, 0.75493660137711338, 0.328048515311972,
		   -1.1006675468043594, -0.2101967497579349, 2.5144813980371978,
		   0.11741758392572139, 0.10341118258170384, 0.49173901637930667,
		   -0.19311826222909464, 1.571035202943557, -2.706372422463168,
		   -5.1269376604274566, 0.54163875959583518, -2.1279709393742494,
		   3.05323695977829, -0.89237826619679683, -0.68889693012446185,
		   2.4046853107229609, -1.4894305740007563, 0.22335321380030143,
		   -1.1867260394483137]

	e_a = [1.4697433960169914, -0.37985323950330468,
		   0.727707309312373, -1.8103688716904811]

	b_av = [0] * 4
	f_a = [-2.1472992909613917, 0.33743332302633383,
		   1.2135711774244253, 0.58329808338968969, -0.685515256200231,
		   -4.4917046251557489, 4.8024703200979388, 0.37043682566991054,
		   5.4144756625124506, -3.1077246609523859, -2.5259666269352339,
		   0.22768046628586849, 1.1004077440458522, -0.79044123867902294,
		   -0.570264481531221, 0.24930837335093, 2.1792806730521415,
		   -1.9999025958887355, 0.50879361640017706, -0.68190177853336709,
		   -3.2077399561421571, 1.7305773053795823, 1.1014295021689073,
		   0.37418313308958567]

	exitg1 = False
	unnamed_idx_0 = 0

	# max min and max input function
	for ixstart in range(3):
		Xp1[ixstart] = X[ixstart] - b[ixstart]

	for ixstart in range(3):
		av[ixstart] = Xp1[ixstart] * b_b[ixstart]

	for ixstart in range(3):
		Xp1[ixstart] = av[ixstart]

	for ixstart in range(3):
		av[ixstart] = Xp1[ixstart] - 1.0

	for ixstart in range(3):
		Xp1[ixstart] = av[ixstart]

	' layer 1 sigmoid symmetric transfer function '
	for ixstart in range(6):
		denominator = 0.0
		for ix in range(3):
			denominator += b_a[ixstart + 6 * ix] * Xp1[ix]

		x[ixstart] = exp(-2.0 * (a[ixstart] + denominator))

	' layer 2 sigmoid symmetric transfer function '
	for ixstart in range(6):
		c_b[ixstart] = c_a[ixstart]
		dv0[ixstart] = 2.0 / (1.0 + x[ixstart]) - 1.0

	for ixstart in range(6):
		denominator = 0.0
		for ix in range(6):
			denominator += d_a[ixstart + 6 * ix] * dv0[ix]

		x[ixstart] = exp(-2.0 * (c_b[ixstart] + denominator))

	' layer 3 '
	for ixstart in range(4):
		numerator[ixstart] = e_a[ixstart]

	for ix in range(6):
		dv0[ix] = 2.0 / (1.0 + x[ix]) - 1.0

	for ix in range(4):
		denominator = 0.0
		for ixstart in range(6):
			denominator += f_a[ix + (ixstart << 2)] * dv0[ixstart]

		b_av[ix] = numerator[ix] + denominator

	' Competitive soft transer function '
	ixstart = 1
	denominator = b_av[0]
	if rtIsNaN(b_av[0]):
		ix = 2
		exitg1 = False
		while ((not exitg1) and (ix < 5)):
			ixstart = ix
			if (not rtIsNaN(b_av[ix - 1])):
				denominator = b_av[ix - 1]
				exitg1 = True
			else:
				ix += 1

	if (ixstart < 4):
		while (ixstart + 1) < 5:
			if (b_av[ixstart] > denominator):
				denominator = b_av[ixstart]
			ixstart += 1

	numerator = bsxfun(b_av, denominator)
	for ixstart in range(4):
		numerator[ixstart] = exp(numerator[ixstart])

	denominator = numerator[0]
	for ixstart in range(3):
		denominator += numerator[ixstart + 1]

	unnamed_idx_0 = denominator
	for ixstart in range(1):
		if (abs(denominator) < 0.000000000000001):
			unnamed_idx_0 = 1.0

	for ixstart in range(4):
		b_av[ixstart] = numerator[ixstart] / unnamed_idx_0

	for ixstart in range(4):
		numerator[ixstart] = b_av[ixstart]

	' output '
	Y = bsxfun(numerator, -1.0)
	for ixstart in range(4):
		Y[ixstart] = Y[ixstart] / 2.0

	temp = 0
	ret = 0
	for i0 in range(4):
		if Y[i0] > temp:
			ret = i0
			temp = Y[i0]

	return ret

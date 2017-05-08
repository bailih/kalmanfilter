from NN.rtGetNaN import *
from NN.rtGetInf import *

rtNaN = 0
rtNaNF = 0
rtInf = 0
rtInfF = 0
rtMinusInf = 0
rtMinusInfF = 0


def rt_InitInfAndNaN(realsize):
	global rtNaN
	global rtNaNF
	global rtInf
	global rtInfF
	global rtMinusInf
	global rtMinusInfF

	rtNaN = rtGetNaN()
	rtNaNF = rtGetNaNF()
	rtInf = rtGetInf()
	rtInfF = rtGetInfF()
	rtMinusInf = rtGetMinusInf()
	rtMinusInfF = rtGetMinusInfF()


def rt_isInf(val):
	global rtInf
	global rtMinusInf
	return 1 if val == rtInf or val == rtMinusInf else 0


def rt_isInfF(val):
	global rtInfF
	global rtMinusInfF
	return 1 if val == rtInfF or val == rtMinusInfF else 0


def rtIsNaN(val):
	return 0


def rtIsNaNF(val):
	return 0

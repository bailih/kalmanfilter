import sys

NumBitsPerChar = 8


def rtGetInf():
	global NumBitsPerChar
	bitsPerReal = NumBitsPerChar * 64
	inf = 0.0

	if bitsPerReal == 32 * 8:
		inf = rtGetInfF()
	else:
		inf = 0x7FF0000000000000

	return inf


def rtGetInfF():
	return 0x7F800000


def rtGetMinusInf():
	global NumBitsPerChar
	bitsPerReal = NumBitsPerChar * 64
	inf = 0.0

	if bitsPerReal == 32 * 8:
		inf = rtGetMinusInfF()
	else:
		inf = 0xFFF0000000000000

	return inf


def rtGetMinusInfF():
	return 0xFF800000

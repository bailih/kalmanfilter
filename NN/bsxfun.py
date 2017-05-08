'''imports here'''


def bsxfun(a, b):
	c = [None] * 4
	for k in range(4):
		c[k] = a[k] - b
	return c

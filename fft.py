from cmath import exp, pi
from math import cos

def fft(x):
	N = len(x)
	if N <= 1:	return x
	even = fft( x[ 0::2 ] )
	odd = fft( x[ 1::2 ] )
	return [ even[ k ] + exp( -2j * pi * k / N ) * odd[ k ] for k in xrange( N / 2 ) ] + \
		[ even[ k ] - exp( -2j * pi * k / N ) * odd[ k ] for k in xrange( N / 2 ) ]

#print fft( [ 1.0, 1.0, 1.0, 1.0, 0.0, 0.0, 0.0, 0.0 ] )

N = 8 
ww = 4000
rate = 10000
y = []
for i in range( N ):
	phase = ( i * 1.0 )/N * 2 * pi
	value = cos( phase ) + 1
	y.append( value )

print fft( y )

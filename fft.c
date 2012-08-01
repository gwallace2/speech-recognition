#include <assert.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/*computes the complex exponential of a complex number.*/
int cmplxExp( double* zr, double* zi, double xr, double xi ) {
	double factor = exp( xr );
	*zr = factor * cos( xi );
	*zi = factor * sin( xi );
	return 0;
}

/*
 * output:
 * 	matrix: the covariance matrix of the data set, which is always a Hermitian matrix
 *
 * input:
 *	vector: the set of vectors for which the covariance will be found
 *	m:	the number of vectors in the set
 *	n:	the number of elements in each vector
 */
int covariance( double** matrix, const double** vector, int m, int n ) {
	
	/*compute mean*/
	double* mean = ( double* )malloc( sizeof( double ) * m );
	int i;
	int j;
	int k;

	for( i = 0; i < n; i++ ) { /*iterating through each vector*/
		printf( "i = %d\n", i );
		for( j = 0; j < m; j++ ) {
			if( i == 0 ) { /*iterating through each element*/
				mean[ j ] = vector[ i ][ j ] / n;
			} else {
				mean[ j ] += vector[ i ][ j ] / n;
			}
		}
	}
	for( i = 0; i < m; i++ ) {
		for( j = 0; j < m; j++ ) {
			for( k = 0; k < n; k++ ) {
				if( k == 0 ) {
					matrix[ i ][ j ] = ( vector[ k ][ i ] - mean[ i ] ) * ( vector[ k ][ j ] - mean[ j ] ) / ( m - 1 );
				} else {
					matrix[ i ][ j ] += ( vector[ k ][ i ] - mean[ i ] ) * ( vector[ k ][ j ] - mean[ j ] ) / ( m - 1 );
				}
			}
		}
	}
	free( mean );
}

/* Pre-emphasis of a sound clip. Doing this compensates for the 6 decibel roll-off of the power spectrum, resulting in more clear
 * Convolves the sound clip with the the digital filter: [1,-0.97].
 *
 * inputs:
 *	input: the sound clip
 *	size:  the number of samples in the sound clip
 *
 * outputs:
 *	result: the resulting, pre-emphasized sound clip
 */
int preemp( double* result, const double* input, int size ) {
	int i;
	result = realloc( result, sizeof( double ) * size + 1 );
	for( i = 0; i < size + 1; i++ ) {
		if( i == 0 ) {
			result[ i ] = input[ i ];
		} else if( i == size ) {
			result[ i ] = -0.97 * input[ i - 1 ];
		} else {
			result[ i ] = input[ i ] - 0.97 * input[ i - 1 ];
		}
	}
}

/* Brackets a portion of a sound clip using a Hamming window.
 *
 * inputs:
 *	data: the sound clip
 *	m:    the starting point of the bracket
 *	n:    the length of the bracket
 * outputs:
 *	segment: the resulting sound segment
 */
int hmwindow( double* segment, const double* data, int m, int n, int size ) {
	int retval = 0;
	assert( m + n < size );
	segment = realloc( segment, n );
	if( segment == 0 ) {
		retval = 1;
	} else {
		int i;
		for( i = m; i < m + n; i++ ) {
			segment[ i ] = ( 0.54 - 0.46 * cos( 2 * M_PI * ( i - m ) / ( n - 1 ) ) ) * data[ i ];
		}
	}
	return retval;
}

/* Fast Fourier transform, using the Tukey-Cooley algorithm.
 *
 * inputs:
 *	xr:   real input vector
 *	xi:   imaginary input vector
 *	samples: number of samples in the clip
 *
 * outputs:
 *	real: real output vector
 *	imag: imaginary output vector
 *
 */
int fft( double* real, double* imag, const double* xr, const double* xi, int samples ) {
	int result;
	if( samples <= 1 ) {

		/*base case. Fourier transform of a single sample is itself*/
		memcpy( real, xr, samples );
		memcpy( imag, xi, samples );
	} else {

		/*declare variables*/
		double* evenSamplesr;
		double* evenSamplesi;
		double* oddSamplesr;
		double* oddSamplesi;
		double* evenResultsr;
		double* evenResultsi;
		double* oddResultsr;
		double* oddResultsi;
		int numEvenSamples;
		int i;

		/*split the samples into even and odd indexed elements*/
		if( samples % 2 == 0 ) {
			numEvenSamples = samples / 2;
		} else {
			numEvenSamples = samples / 2 + 1;
		}
		evenSamplesr = ( double* )malloc( sizeof( double ) * numEvenSamples );
		evenSamplesi = ( double* )malloc( sizeof( double ) * numEvenSamples );
		oddSamplesr = ( double* )malloc( sizeof( double ) * samples / 2 );
		oddSamplesi = ( double* )malloc( sizeof( double ) * samples / 2 );
		for( i = 0; i < samples; i++ ) {
			if( i % 2 == 0 ) {
				evenSamplesr[ i / 2 ] = xr[ i ];
				evenSamplesi[ i / 2 ] = xi[ i ];
			} else {
				oddSamplesr[ i / 2 ] = xr[ i ];
				oddSamplesi[ i / 2 ] = xi[ i ];
			}
		}

		/*take the Fourier transform of each set*/
		evenResultsr = ( double* )malloc( sizeof( double ) * numEvenSamples );
		evenResultsi = ( double* )malloc( sizeof( double ) * numEvenSamples );
		oddResultsr = ( double* )malloc( sizeof( double ) * samples / 2 );
		oddResultsi = ( double* )malloc( sizeof( double ) * samples / 2 );
		fft( evenResultsr, evenResultsi, evenSamplesr, evenSamplesi, samples / 2 );
		fft( oddResultsr, oddResultsi, oddSamplesr, oddSamplesi, samples / 2 );

		/*compute the real and imaginary results*/
		real = realloc( real, sizeof( double ) * samples );
		imag = realloc( imag, sizeof( double ) * samples );
		for( i = 0; i < samples; i++ ) {
			double zr;
			double zi;
			int k;
			cmplxExp( &zr, &zi, 0, -2 * M_PI * k / samples );
			if( i < samples / 2 ) {
				k = i;
				real[ i ] = evenSamplesr[ k ] + zr * oddSamplesr[ k ] - zi * oddSamplesi[ k ];
				imag[ i ] = evenSamplesi[ k ] + zr * oddSamplesi[ k ] + zi * oddSamplesr[ k ];
			} else {
				k = i % ( samples / 2 );
				real[ i ] = evenSamplesr[ k ] - zr * oddSamplesr[ k ] + zi * oddSamplesr[ k ];
				imag[ i ] = evenSamplesr[ k ] - zr * oddSamplesr[ k ] - zi * oddSamplesr[ k ];
			}
		}
		free( evenSamplesr );
		free( evenSamplesi );
		free( oddSamplesr );
		free( oddSamplesi );
		free( evenResultsr );
		free( evenResultsi );
		free( oddResultsr );
		free( oddResultsi );
	}
	return 0;
}

/*
 * Get the mel frequency from a given linear frequency.
 *
 * input:
 *	freq: linear frequency
 *
 * returns:
 *	the mel frequency
 */
double getMelFreq( double freq ) {
	return 2595 * log10( 1 + freq / 700 );
}

/*
 * Compute the total energy of the signal. This is a simple task, as it just
 * involves taking the sum of the square of all samples.
 *
 * inputs:
 *	signal: the signal in question
 *	numSamples: the number of samples in the signal
 *
 * returns:
 *	the energy
 */
double energy( const double* signal, int numSamples ) {
	int i;
	double result = 0.0;
	for( i = 0; i < numSamples; i++ ) {
		result = signal[ i ] * signal[ i ];
	}
	return result;
}

int main( int argc, char* argv[] ) {

	/*unit test for complex exponential*/
	double zr;
	double zi;
	clock_t start = clock() / ( CLOCKS_PER_SEC / 1000 );
	cmplxExp( &zr, &zi, 0, M_PI );
	printf( "exp( 1i*pi ) = %lf + i%lf\n", zr, zi );
	cmplxExp( &zr, &zi, 1, M_PI );
	printf( "exp( 1 + 1i*pi ) = %lf + i%lf\n", zr, zi );
	cmplxExp( &zr, &zi, 0, 0 );
	printf( "exp( 0 ) = %lf + i%lf\n", zr, zi );
	clock_t end = clock(); 
	printf( "took %d to compute complex exponential\n", ( start ) );

	/*unit test for covariance*/
	double* v1 = ( double* )malloc( sizeof( double ) * 3 );
	double* v2 = ( double* )malloc( sizeof( double ) * 3 );
	double* v3 = ( double* )malloc( sizeof( double ) * 3 );
	double* v4 = ( double* )malloc( sizeof( double ) * 3 );
	double* v5 = ( double* )malloc( sizeof( double ) * 3 );
	double** C = ( double** )malloc( sizeof( double* ) * 3 );
	int i;
	int j;
	for( i = 0; i < 3; i++ ) {
		C[ i ] = ( double* )malloc( sizeof( double ) * 3 );
	}
	v1[ 0 ] = 4.0;
	v1[ 1 ] = 3.5;
	v1[ 2 ] = 3.2;

	v2[ 0 ] = 3.0;
	v2[ 1 ] = 7.8;
	v2[ 2 ] = 10.4;

	v3[ 0 ] = -8.0;
	v3[ 1 ] = 4.0;
	v3[ 2 ] = 5.0;

	v4[ 0 ] = 1.6;
	v4[ 1 ] = 2.7;
	v4[ 2 ] = -1.5;

	v5[ 0 ] = 4.6;
	v5[ 1 ] = 3.4;
	v5[ 2 ] = 3.1;
	double** vv = ( double** )malloc( sizeof( double* ) * 5 );
	vv[ 0 ] = v1;
	vv[ 1 ] = v2;
	vv[ 2 ] = v3;
	vv[ 3 ] = v4;
	vv[ 4 ] = v5;
	covariance( C, vv, 3, 5 );
	for( i = 0; i < 3; i++ ) {
		printf( "[ " );
		for( j = 0; j < 3; j++ ) {
			printf( "%lf ", C[ i ][ j ] );
		}
		printf( "]\n" );
	}
	free( vv );
	for( i = 0; i < 3; i++ ) {
		free( C[ i ] );
	}
	free( C );
	free( v5 );
	free( v4 );
	free( v3 );
	free( v2 );
	free( v1 );

	/*Unit test for pre-emphasis*/
	double* clip = ( double* )malloc( sizeof( double ) * 1103 );
	double* result = ( double* )malloc( sizeof( double ) * 1104 );
	srand( time( 0 ) );
	for( i = 0; i < 1103; i++ ) {
		clip[ i ] = 2 * ( ( double )rand() / ( double )RAND_MAX ) - 1;
		printf( "clip[ %d ] = %lf\n", i, clip[ i ] );
	}
	preemp( result, clip, 1103 );
	printf( "results:\n" );
	for( i = 0; i < 1104; i++ ) {
		printf( "result[ %d ] = %lf\n", i, result[ i ] );
	}
	return 0;
}

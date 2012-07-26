#include <math.h>
#include <stdio.h>
#include <stdlib.h>

/*computes the complex exponential of a complex number.*/
int cmplxExp( double* zr, double* zi, double xr, double xi ) {
	double factor = exp( xr );
	*zr = factor * cos( xi );
	*zi = factor * sin( xi );
	return 0;
}

/*
 * output:
 * 	matrix: the covariance matrix, which is always a Hermitian matrix
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
	
}

int preemp( double* result, const double* input, int size ) {
	int i;
	result = ( double* )malloc( sizeof( double ) * size + 1 );
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

int fft( double* real, double* imag, const double* x, int samples ) {
	int result;
	if( samples <= 1 ) {
		real = 
	} else {
		double* evenSamples;
		double* oddSamples;
		int numEvenSamples;
		int i;
		if( samples % 2 == 0 ) {
			numEvenSamples = samples / 2;
		} else {
			numEvenSamples = samples / 2 + 1;
		}
		evenSamples = ( double* )malloc( sizeof( double ) * numEvenSamples );
		oddSamples = ( double* )malloc( sizeof( double ) * samples / 2 );
		for( i = 0; i < samples; i++ ) {
			
		}
	}
	return 0;
}

int main( int argc, char* argv[] ) {
	double zr;
	double zi;
	cmplxExp( &zr, &zi, 1, M_PI );
	printf( "exp( 1 + 1*pi ) = %lf + i%lf\n", zr, zi );
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
	return 0;
}

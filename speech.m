function [ ] = speech( speechfile, frameWidth, frameStep, numTriangles )
    [ clip fs ] = ( wavread( speechfile ) );
    N = ceil( fs * frameWidth / 1000 );
    M = ceil( fs * frameStep / 1000 );
    clip = clip';
    
    %pre-emphasis
    filt = [ 1 0.95 ];
    s = conv( clip, filt );
    
    %get the first frame using a Hamming window
    nrange = [0:1:N-1];
    w = 0.54 - 0.46 * cos( ( 2 * pi * nrange ) / ( N - 1 ) );
    N
    length( s )
    v0 = s(1:1:N) .* w;
    
    %get frequency spectrum of frame
    Vhat0 = fftshift( fft( v0 ) );
    
    %get triangle center frequencies, equally spaced on the mel scale
    %(should it include the first triangle?
    K = numTriangles;
    fmax = fs/2;
    fmaxmel = 2595 * log10( 1 + fmax/700 );
    delmel = fmaxmel/(K+1);
    krange = [0:1:K-1];
    fcmel = ( krange + 1 ) * delmel;
    fc = 700 * ( 10 .^ ( fcmel / 2595 ) - 1 );
    
    %
    W = linspace( 
end
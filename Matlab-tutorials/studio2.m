%% Solution to 10.19
%
% Tomas McKelvey, November 2007
%
% Spec. FS = 100
Fs = 100;
Fp = 15;
Fst = 20;
% a) Firpm design (equiripple design)

%%
h = firpm(37,[0 Fp/Fs Fst/Fs 1/2]*2,[1 1 0 0],[1 2]);
NN = 2^10;
F = (0:NN-1)/NN*Fs;
H = abs([fft(h.',NN)]);
plot(F,20*log10(H)); % Zoom in plot to check specifications

% Stopband ripple in decibel
20*log10(0.01)
% Passband ripple in decibel
20*log10(1.02)

%% There is an automatic function which provides an order estimate
% Note though that the order had to be increased a factor 2 to meet the
% specs.
[n,f,a,w] = firpmord( [15 20], [1 0], [0.02 0.01], 100);
h2 = firpm(n+2,f,a,w);
H = abs([fft(h.',NN),fft(h2.',NN)]);
plot(F,20*log10(H)); % Zoom in plot to check specifications

% b) eq 10.2.94 gives
Mhat1 = (-20*log10(sqrt(0.01*0.02))-13)/14.6/((20-15)/100) + 1
% while 10.2.95 is the equation used in firpmord


%% Now design using window technique
h3 = fir1(37,15/Fs*2,hamming(37+1));
H = abs([fft(h.',NN),fft(h3.',NN)]);
plot(F,20*log10(H)); % Zoom in plot to check specifications

% The high stop-band attenuation makes the transition region large and
% filter does not meet specs.

%% Try a chebyshev window
h4 = fir1(37,17.5/Fs*2,chebwin(37+1,33));
H = abs([fft(h.',NN),fft(h4.',NN)]);
plot(F,20*log10(H)); % Zoom in plot to check specifications
% Almost meets the specs.

%% d) IIR elliptic filter 

[b,a] = ellip(5,0.17,40,Fp/Fs*2);
H5 = freqz(b,a,F*2*pi/Fs);
plot(F,20*log10([H(:,1) H5.']))
%%
% check the poles of the filter
p = roots(a);
[p abs(p) angle(p)*180/pi]
plot(p,'*'); axis square;
axis([-1 1 -1 1]);
cc = exp(i*2*pi*(0:100)/100);
plot(cc); hold on ; %plot a circle 
plot(p,'*'); axis square;
axis([-1 1 -1 1]);
hold off

%% e) 
% The fir filter (order 37) has 38 coefficients
% To produce one sample output 38 multiplications and 37 additions are
% needed and 37 memory elements

% The elliptical IIR filter implemented in second order sections requires 
% two second order sections and one first order section. 
% eight multiplications 10 additions and 5 memory elements
% This particular filter has severeal coefficents which are one and thus
% does not need any multiplyer for those coefficients.

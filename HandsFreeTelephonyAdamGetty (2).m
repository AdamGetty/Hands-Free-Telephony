%% Stage 1 -----------------------------------------------------------------------------------------------
%Read in and process audio file to create signal x[n]
[x,Fs] = audioread('FarEndSignal.wav');			%Signal 'x[n]' is sampled at Fs = 8000Hz
Fs= 8000;
xf = x;
[y,N] = room(xf,201747002)	;					%Obtain y[n] signal. N parameter based on student registration number
y(numel(xf)) = 0;								%Extend y vector with zeros to match lengths.

%% Stage 2 -----------------------------------------------------------------------------------------------
%Create filter g[n], which is 28 fold expansion of g3[n] = (-1)^n*h3[n]
 
%%h3[n] Filter
h3 = fir1(23,0.4);            					%Created a window or order 23 with a cuttoff frequency of Fsample/2*0.4;

%%g3[n] Filter
nrange = (0:1:length(h3)-1);       				%Define a sample range for (-1)^n = cos(pi*n)
h5 = cos(pi*nrange);	
g3 = h5.*h3;          							%Multiply the signals together to satisfy g3 = (-1)^n*h3[n] = cos(pi*n)*fir1(23,0.4)

%N-form expansion of g3 to produce g[n]
g = zeros(1,N*length(g3));
g(1:N:end) = g3;                                %Create the N-fold expansion of h3 by 28
figure;
subplot(2,1,1);
freqz(g);
ylim([-100 10]);
title('Magnitude Response(dB) of g[n]');
subplot(2,1,2);
impz(g);
title('Impulse Response(dB) of g[n]'); 

%% 
%%Stage 3 ----------------------------------------------------------------------------------------------
%Pass signal xf[n] through filter g[n]

xf = filter(g,1,x);								%xf[n] %%This could also possibly be done with manual filtering from labs??
[y,N] = room(xf,201747002);


%% 
%%Stage 4-------------------------------------------------------------------------------------------------
%Create Filter h[n] where, h[n] = h3[28*n]

%Comb Filter Design for h[n]
h = zeros(1,N*length(h3));
h(1:N:end) = h3;			 					%Create the N-fold expansion of h3 by 28
figure;
subplot(2,1,1);
freqz(h);
ylim([-100 10]);
title('Magnitude Response(dB) of h[n]');
subplot(2,1,2);
impz(h);
title('Impulse Response(dB) of h[n]'); 

%%
%%Stage 5--------------------------------------------------------------------------------------------------
%Send Near end, noisy signal through h[n] filter.

yf = filter(h,1,y);

%%Uncomment to plot resultant signal after filtering
% %FIR Matlab Filter Figure
% t = 1:1:length(yf);					%Range for plot
% figure; plot(t,yf); 
% title('FIR Matlab Filter');
% xlabel('Time (s)'); ylabel('Amplitude'); %Create axis labels % Change from samples to seconds

%%
%%Stage 6-------------------------------------------------------------------------------------2------------
%Implementation of function without the use of Matlab functions
tic;							%Start Timing 

%%Manual Hamming Window Implementation                      
FL=24;  						%Length of Filter
n=0:1:FL-1; 					%Defines Sample Range
p = n-(FL-1)/2;                 %Range of sinc Function
fc=0.4;  						%Will have 6dB cuoff at 0.4 normalized frequency.
Z=sin(2*pi*(fc/2)*p)./(pi*p); 	%Define truncated Sinc function		

s = 2*pi*(n/(FL-1));      		%Frequency For Hamming Window
w = 0.54-0.46*cos(s);  			%Define Hamming window function 

h3manualfilt = Z.*w;   			%Multiplication of Hamming Window and sin function. (Giving filter coefficients)

%%Perform N-fold
hmanual = zeros(1,N*length(h3manualfilt));
hmanual(1:N:end) = h3manualfilt;
figure;
subplot(2,1,1);
freqz(hmanual);
ylim([-100 10]);
title('Magnitude Response(dB) of h[n](manual)');
subplot(2,1,2);
impz(hmanual);
title('Impulse Response(dB) of h[n](manual)'); 

%%Create iteration to perform filtering on signal 

%FIR Filter Implementation
L = length(y); 					% length of simulation
% data = y

%FIR Parameter Initializations
N = length(hmanual); 				% length of TDL
x_tdl = zeros(N,1); 				% TDL initialisation
y_manualfilt =	zeros(L,1);			% TDL Output initialisation
t = 1:1:length(y);					%Range for plot

for n = 1:L 						% iteration
x_tdl = [y(n); x_tdl(1:N-1)]; 		% update
y_manualfilt(n) = hmanual*x_tdl;	%Filter Output Calaculation /  Can calcuate as hman[1:N:length(hamn)]*x_tdl[1:N:length(x_tdl)] To reduce computation!!

end 							% end iteration
toc; 							%End and read timing Window
looptime = toc;

%%Uncomment to plot resultant signal after filtering
%%FIR Manual Filter Figure
%figure; plot(t,y_manualfilt); 
%title('FIR Manual Filter');
%xlabel('Time (s)'); ylabel('Amplitude'); %Create axis labels % Change from samples to seconds

%%
%%Stage 7---------------------------------------------------------------------------------------
%%Inspect y[n] for noise floor frequency
%Write y to file for information
filename = 'NearEndSignal.wav';
audiowrite('NearEndSignal.wav',y,Fs);
[ns,Fs] = audioread('NearEndSignal.wav');			%Signal 'y[n]' is sampled at Fs = 8000Hz
infons = audioinfo('NearEndSignal.wav');			%Extract information from audio file

%Create relevant time vector for plotting y
ty = 0:seconds(1/Fs):seconds(infons.Duration);		%Call to the infons data to extract signal duration
ty = ty(1:end-1);

%y[n] for time domain calculation of noise frequency (Used zoom function)
figure;
subplot(2,1,1);
plot(ty,y);
xlim([0 seconds(0.1)]); ylim([-0.25 0.25]);
title('Near End Signal for Noise floor frequency');
xlabel('Time(s)'); ylabel('y[n]');

%Run the filtered signal through a Fast Fourier Transform to determine the frequency of the noise floor

Y = fft(y);				

%Plot the fft of the signal
  
a1 = length(y);
freqd = Fs*(0:(a1/2))/a1;    %Define frequency domain

%P = two-sided spectrum. Compute single sided spectrum using P3.

P = abs(Y/a1);
P2 = P(1:a1/2+1);
P2(2:end-1) = 2*P2(2:end-1);

%Plot the Fourier Transform

%figure;
subplot(2,1,2);
plot(freqd,P2);
title('Spectrum of yf[n]'); 
xlabel('Frequency(Hz)'); 
ylabel('Magnitude');

%%
%%Stage 8-----------------------------------------------------------------------------------------------
%%Plot magnitude response of notch filter design where (Transfer function from logbook)
%%Design notch filter to eliminate noise floor using Fourier approximation (Use transfer function from logbook)
f0 = 857.1429;           	%Define Notch Cutoff  (Using Fast Fourier Trnsfomr)

%Co-efficients
r = 0.9;
b0 = 1; b1 = -2*cos(2*pi*f0/Fs); b2 = 1;
a0 = 1; a1 = -2*r*cos(2*pi*f0/Fs); a2 = r^2;

%Co-efficient Vector for Matlab processing
b = [b0;b1;b2];
a = [a0;a1;a2];

figure;
subplot(2,1,1);
freqz(b,a);
title('Magnitude Response(dB) of q[n]');
subplot(2,1,2);
impz(b,a);
title('Impulse Response(dB) of q[n]'); 

%%
%%Stage 9--------------------------------------------------------------------------------------

%Apply notch filter to y signal
yff = filter(b,a,yf);
sound(yff,8000)

% %%Uncomment for signal using Stage 6 manual filter
% %Applying Notch filter to manually filterd y signal.
% yffmanual = filter(b,a,y_manualfilt);
% sound(yffmanual, 8000);




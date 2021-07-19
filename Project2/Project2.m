% AUTHOR: KWM
% Reference: http://aaronscher.com/wireless_com_SDR/RTL_SDR_AM_spectrum_demod.html
% DATE: 10/3/18





%%%%%%%%%%%%%%%%%%%%%% BRINGING IT ALL TOGETHER %%%%%%%%%%%%%%%%%%%%
fc = 108; % MHz, carrier of KissFM Worcester
fs = 2.5; % MHz, sample rate of RTL_SDR used to collect 4 sec of data
plot_time = 0.002*2.5E6; % amount of time to plot of the 4 sec of data

fid = fopen('KissFM_4sec.dat','rb');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
y = fread(fid,'uint8=>double');
y = y-127.5;
y = y(1:2:end) + 1i*y(2:2:end); % convert to complex values
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

figure(1)
plot_FFT_IQ(y,1,plot_time,fs,fc,'Spectrum of observed signal');

% plot one shows a 1.78 MHz frequency shift from our 108 MHz expected
y_shifted=y.*transpose(exp(-1i*2*pi*(1.78E6)*[1:1:length(y)]/2.5E6)); 
figure(2)
plot_FFT_IQ(y_shifted,1,plot_time,fs,fc,'Spectrum of shifted signal'); 

% decimate
d_f = 8; % first decimation factor
d = decimate(y_shifted,d_f,'fir'); 
figure(3)
plot_FFT_IQ(d,1,plot_time/d_f,fs/d_f,fc,'Spectrum of decimated signal'); 

% Demodulate
[y_FM_demodulated] = FM_IQ_Demod(d); %d is the decimated signal
figure(4)
plot_FFT_IQ(y_FM_demodulated,1,.05*2.5E6/d_f,fs/d_f,0,'Spectrum of demodulated signal');

% Time to listen!
d_f2 = 10; % final decimation factor
df = decimate(y_FM_demodulated,d_f2,'fir'); % decimate to bring into audio range
figure(5)
plot_FFT_IQ(df,1,.05*2.5E6/d_f/d_f2,fs/d_f/d_f2,0,'Spectrum of audio signal');
sound(df,2.5E6/d_f/d_f2);






% AUTHOR: KWM

% Specify carrier frequency
w_c = 2.*pi.*5;
% Generate a discrete version of a random continuous analog
% waveform using a Uniform Random Number Generator and
% an interpolation function to smooth out the result
L = 100;  % Length of the overall transmission
M = 10;   % Upsampling factor for generating analog waveform
analog_wavefm = interp(rand(1,(L/M)),M);






plot(analog_wavefm)
% figure(1)
% subplot(2,1,1)
% plot()
% subplot(2,1,2)
% plot()








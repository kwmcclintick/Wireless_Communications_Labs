% AUTHOR: Alex Wyglinski, Kyle McClintick and Gigi Dong
% DATE: 10/3/18
close all;
clear;
clc;

% Define constants
L = 100; % Length of the overall transmission
M = 10; % Pulse duration for rectangular pulse train
N = 10; % Upsampling factor for generating analog waveform

% Generate a discrete version of a random continuous analog
% waveform using a Uniform Random Number Generator and
% an interpolation function to smooth out the result
analog_wavefm = interp((2*rand(1,(L/M))-1),M);
% Generate a rectangular pulse train of samples
impulsetrain_wavefm = reshape(ones(N,1)*rem(1:1:(L/N),2),[1,L]);

% Plotting Base signals
% figure(1)
% subplot(2,1,1)
% plot(analog_wavefm); title('Random Analog Waveform')
% ylim([-2 2]); ylabel('Signal Amplitude'); xlabel('Discrete Time (n)');
% subplot(2,1,2)
% stem(impulsetrain_wavefm); ylabel('Signal Amplitude');
% xlabel('Discrete Time (n)'); title('Impulse Train');

% PAM
natural_PAM = analog_wavefm .* impulsetrain_wavefm;
flattop_PAM = repelem(natural_PAM(1:N:end), N);

% plotting PAM
% figure(2)
% subplot(2,1,1)
% stem(natural_PAM); title('Naturally Sampled PAM');
% ylim([-2 2]); ylabel('Signal Amplitude'); xlabel('Discrete Time (n)');
% subplot(2,1,2)
% stem(flattop_PAM); ylabel('Signal Amplitude');
% xlabel('Discrete Time (n)'); title('Flat-Top PAM');

% Quantizing
q = zeros(1,L);
[ind,quantv] = quantiz(downsample(analog_wavefm,N), ...
    [-0.8 -0.6 -0.4 -0.2 0 0.2 0.4 0.6 0.8], ...
    [-0.9 -0.7 -0.5 -0.3 -0.1 0.1 0.3 0.5 0.7 0.9]);
for i=1:N-1
    [d, ix] = min(abs(quantv -   flattop_PAM(1+((i-1)*N))  ));
    q(1+((i-1)*N):i*N) = quantv(ix);
end
q = q .* impulsetrain_wavefm;
residual = flattop_PAM - q;

% 
figure(3)
subplot(5,1,1)
plot(analog_wavefm); title('Random Analog Waveform'); xlabel('Time');
ylabel('Signal Amplitude');
subplot(5,1,2)
stem(natural_PAM); title('Natural PAM'); xlabel('Discrete Time (n)');
ylabel('Signal Amplitude');
subplot(5,1,3)
stem(flattop_PAM); title('Flat-Top PAM'); xlabel('Discrete Time (n)');
ylabel('Signal Amplitude');
subplot(5,1,4)
stem(q); title('PCM'); xlabel('Discrete Time (n)');
ylabel('Signal Amplitude');
subplot(5,1,5)
stem(residual); ylabel('Residual Amplitdue'); xlabel('Discrete Time (n)');
title('Quantization Error');


% Generate your own line codes for the binary string `11001100'
L_lc = 20; % Line coding pulse duration
bin_str = round(rand(1,L));
upnrz1 = ones(1,L_lc);
upnrz0 = zeros(1,L_lc);
upnrz_wavefm = [];
for ind = 1:1:length(bin_str)
    if (bin_str(ind) == 1)
        upnrz_wavefm = [upnrz_wavefm upnrz1];
    else
        upnrz_wavefm = [upnrz_wavefm upnrz0];
    end
end


% Visualize 2 periods of the eye diagram
% eyediagram(upnrz_wavefm,2*L_lc,L_lc,L_lc/2);
% ylim([-0.5 1.5]);


% Bringing it all together
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Fs = 44100; % 44.1kHz sample rate, as most sounds is <20 kHz (Nyquist)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[notdone, Fs_1] = audioread('chickens.wav');
notdone = notdone(:,1);
%[notdone, Fs_3] = audioread('elbow.wav');

% Plot raw data
% figure(1)
% stem(notdone); title('Raw Mystery Soundfile'); xlabel('Discrete Time (n)');
% ylabel('Signal Amplitude');

% Define constants
l = length(notdone); % Length of the overall transmission
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n = 2; % downsample/upsample rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% encode the waveforms
encode_impulsetrain = reshape(ones(n,1)*rem(1:1:(l/n),2),[l,1]);
Vmax = max(notdone);
Vmin = min(notdone);

natural = notdone .* encode_impulsetrain;
% figure(2)
% stem(natural); title('Natural PAM'); xlabel('Discrete Time (n)');
% ylabel('Signal Amplitude');
flattop = repelem(natural(1:n:end), n);
% figure(3)
% stem(flattop); title('Flat-Top PAM'); xlabel('Discrete Time (n)');
% ylabel('Signal Amplitude');

% Quantization codebook
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
codes = 2^10; % number of codes in the codebook
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

codebook = linspace(Vmin, Vmax, codes);
partition = linspace(Vmin+abs(codebook(2)-codebook(1)), Vmax, codes-1);
[index,quant_file] = quantiz(flattop,partition,codebook); % Quantize.

% Convert to binary where Vmax is all zero, -Vmax is all 1, and each
% descending value increments the least sig bit
binary_stream = zeros(1, l/n*2);
nz_encoded = quant_file(quant_file ~= 0);
binary_assignments = de2bi(linspace(0,codes-1,codes));

for i=1:length(nz_encoded)
    curr_assign = (nz_encoded(i) == codebook);
    binary_stream(1+(i-1)*log2(codes):i*log2(codes)) = ...
        binary_assignments(curr_assign,1:end);
end
% figure(4)
% stem(binary_stream); title('Mystery Binary Stream'); xlabel('Discrete Time (n)');
% ylabel('Signal Amplitude');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         STUDENTS TASK BEGINS HERE GIVEN ASSUMPTIONS SUPPLIED
%                      binary_stream
%                      Vmax, Vmin to make codebook
%                      rules for making binary assignments
%         THEY MUST ADJUST
%                      codes, n, Fs
bin = load('binary_stream_1.mat');

% Convert back to amplitude levels from binary
decoded_binary = zeros(1, length(bin.binary_stream)/log2(codes));
for i=1:length(decoded_binary)
    curr_code = (bin.binary_stream(1+(i-1)*log2(codes):i*log2(codes)) == ...
        binary_assignments(:,1:end));
    decoded_binary(i) = codebook(sum(curr_code,2) == log2(codes));
end
Vmax = max(decoded_binary);
Vmin = min(decoded_binary);

% play files according to sampling rate
s = load('sb_jr_p1_file1.mat');
sound(s, Fs);

% audiowrite('mystery_3.wav',decoded_binary, Fs);










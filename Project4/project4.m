% AUTHOR: KWM
close all;
clear;

% Define simulation parameters
N_bits = 10000; % Total number of bits for entire transmission
f_c = 2; % Carrier frequency
f_s = 10; % Sampling frequency
f_d = 1; % Digital system sampling frequency
beta = 0.5; % Roll-off factor for raised cosine filter
sigma_n = 0.25; % Noise standard deviation

% STEP 1: Modulate BASK and 4-ASK waveforms
% Generate random binary data stream
bin_data = round(rand(1,N_bits));

% Assignment of binary sets to amplitude values
%ask
M = 4;
bin = reshape(bin_data,2,N_bits/2);
ampl_ask = zeros(1,N_bits/2);
ampl_ask(find( (bin(1,:) == 0)&(bin(2,:) == 0) )) = -3;
ampl_ask(find( (bin(1,:) == 1)&(bin(2,:) == 0) )) = -1;
ampl_ask(find( (bin(1,:) == 0)&(bin(2,:) == 1) )) = 1;
ampl_ask(find( (bin(1,:) == 1)&(bin(2,:) == 1) )) = 3;
%psk
ampl_psk = zeros(1,N_bits/2);
ampl_psk(find( (bin(1,:) == 0)&(bin(2,:) == 0) )) = complex(1,0);
ampl_psk(find( (bin(1,:) == 1)&(bin(2,:) == 0) )) = complex(0,1);
ampl_psk(find( (bin(1,:) == 0)&(bin(2,:) == 1) )) = complex(-1,0);
ampl_psk(find( (bin(1,:) == 1)&(bin(2,:) == 1) )) = complex(0,-1);
%qam
ampl_qam = zeros(1,N_bits/2);
ampl_qam(find( (bin(1,:) == 0)&(bin(2,:) == 0) )) = complex(-1,1);
ampl_qam(find( (bin(1,:) == 1)&(bin(2,:) == 0) )) = complex(-1,-1);
ampl_qam(find( (bin(1,:) == 0)&(bin(2,:) == 1) )) = complex(1,1);
ampl_qam(find( (bin(1,:) == 1)&(bin(2,:) == 1) )) = complex(1,-1);
%bpsk


% Apply raised cosine pulse shaping filter

rolloff = 0.35;     % Rolloff factor
span = 6;           % Filter span in symbols
sps = 8;            % Samples per symbol
rcos = rcosdesign(rolloff, span, sps);
ampl_ask_rcos = upfirdn(ampl_ask, rcos, sps);
ampl_psk_rcos = upfirdn(ampl_psk, rcos, sps);
ampl_qam_rcos = upfirdn(ampl_qam, rcos, sps);

% modulate to carrier frequency f_c
tx_ask_wavefm = ampl_ask_rcos.*...
        cos(2.*pi.*f_c.*(0:(1/f_s):((length(ampl_ask_rcos)-1)*(1/f_s))));
tx_psk_wavefm = ampl_psk_rcos.*...
        cos(2.*pi.*f_c.*(0:(1/f_s):((length(ampl_psk_rcos)-1)*(1/f_s))));
tx_qam_wavefm = ampl_qam_rcos.*...
        cos(2.*pi.*f_c.*(0:(1/f_s):((length(ampl_qam_rcos)-1)*(1/f_s))));
    
% STEP 2: Introduce Passband AWGN to transmission (both I and Q components)
noise_ask = (sigma_n/sqrt(2)).*randn(1,length(tx_ask_wavefm)).*...
    cos(2.*pi.*f_c.*(0:(1/f_s):((length(ampl_ask_rcos)-1)*(1/f_s))))+...
    (sigma_n/sqrt(2)).*randn(1,length(tx_ask_wavefm)).*...
    sin(2.*pi.*f_c.*(0:(1/f_s):((length(ampl_ask_rcos)-1)*(1/f_s))));
noise_psk = (sigma_n/sqrt(2)).*randn(1,length(tx_psk_wavefm)).*...
    cos(2.*pi.*f_c.*(0:(1/f_s):((length(ampl_psk_rcos)-1)*(1/f_s))))+...
    (sigma_n/sqrt(2)).*randn(1,length(tx_psk_wavefm)).*...
    sin(2.*pi.*f_c.*(0:(1/f_s):((length(ampl_psk_rcos)-1)*(1/f_s))));
noise_qam = (sigma_n/sqrt(2)).*randn(1,length(tx_qam_wavefm)).*...
    cos(2.*pi.*f_c.*(0:(1/f_s):((length(ampl_qam_rcos)-1)*(1/f_s))))+...
    (sigma_n/sqrt(2)).*randn(1,length(tx_qam_wavefm)).*...
    sin(2.*pi.*f_c.*(0:(1/f_s):((length(ampl_qam_rcos)-1)*(1/f_s))));

% introduce frequency offset to a signal to export for the project
% fo = 0.2;
% freq_offset_qam = tx_qam_wavefm.* ...
%     exp(-1i*2*pi*fo/f_s*[1:1:length(tx_qam_wavefm)]);
% rx_mystery_waveform = 2*(freq_offset_qam + noise_qam);
% rx_mystery_waveform = 2*(tx_qam_wavefm + noise_qam).* ...
%     exp(-1i*2*pi*fo/f_s*[1:1:length(tx_qam_wavefm)]);

% Waveforms to detect:
rx_ask_wavefm = 2*(tx_ask_wavefm + noise_ask);
rx_psk_wavefm = 2*(tx_psk_wavefm + noise_psk);
rx_qam_wavefm = 2*(tx_qam_wavefm + noise_qam);
scatter(real(rx_qam_wavefm), imag(rx_qam_wavefm))

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% MODULATION DETECTION %%%
% remove carrier ask

no_carrier_ask = rx_ask_wavefm.*cos(2.*pi.*f_c.*(0:(1/f_s):...
            ((length(ampl_ask_rcos)-1)*(1/f_s)))) + ...
            1i*rx_ask_wavefm.*sin(2.*pi.*f_c.*(0:(1/f_s):...
            ((length(ampl_ask_rcos)-1)*(1/f_s))));
% psk
no_carrier_psk = rx_psk_wavefm.*cos(2.*pi.*f_c.*(0:(1/f_s):...
            ((length(ampl_psk_rcos)-1)*(1/f_s)))) + ...
            1i*rx_psk_wavefm.*sin(2.*pi.*f_c.*(0:(1/f_s):...
            ((length(ampl_psk_rcos)-1)*(1/f_s))));
% qam
no_carrier_qam = rx_qam_wavefm.*cos(2.*pi.*f_c.*(0:(1/f_s):...
            ((length(ampl_qam_rcos)-1)*(1/f_s)))) + ...
            1i*rx_qam_wavefm.*sin(2.*pi.*f_c.*(0:(1/f_s):...
            ((length(ampl_qam_rcos)-1)*(1/f_s))));

% Remove double frequency terms of I and Q components via LPF
filt_coeffs = fir1(61,0.25);
lpfed_ask = conv(filt_coeffs,no_carrier_ask);
lpfed_psk = conv(filt_coeffs,no_carrier_psk);
lpfed_qam = conv(filt_coeffs,no_carrier_qam);
% truncate
lpfed_ask = lpfed_ask(floor(length(filt_coeffs)/2):1:end-...
    floor(length(filt_coeffs)/2));
lpfed_psk = lpfed_psk(floor(length(filt_coeffs)/2):1:end-...
    floor(length(filt_coeffs)/2));
lpfed_qam = lpfed_qam(floor(length(filt_coeffs)/2):1:end-...
    floor(length(filt_coeffs)/2));

% matched filter
downsampled_ask = upfirdn(lpfed_ask, rcos, 1, sps);
downsampled_psk = upfirdn(lpfed_psk, rcos, 1, sps);
downsampled_qam = upfirdn(lpfed_qam, rcos, 1, sps);

% truncate tails from convolution in upfirdn
downsampled_ask = downsampled_ask(7:1:end-6);
downsampled_psk = downsampled_psk(7:1:end-6);
downsampled_qam = downsampled_qam(7:1:end-6);

% create error vector magnitude object
evm = comm.EVM();

% decide which hypothesis to test here
ref = ampl_qam';

% lowest rmsEVM is the decided modulation scheme
rmsEVM_ask = evm(ref,downsampled_ask')
rmsEVM_psk = evm(ref,downsampled_psk')
rmsEVM_qam = evm(ref,downsampled_qam')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% DECODING TO OBTAIN ORIGINAL BITS %%%
% demodulate
demod_ask = pamdemod(downsampled_ask, M);
demod_psk = pskdemod(downsampled_psk, M);
demod_qam = qamdemod(downsampled_qam, M);

% back to binary from decimal
% 4-ask
ask_results = zeros(2,length(demod_ask));
ask_results(1, find( demod_ask == 0 )) = 0;
ask_results(2, find( demod_ask == 0 )) = 0;
ask_results(1, find( demod_ask == 1 )) = 1;
ask_results(2, find( demod_ask == 1 )) = 0;
ask_results(1, find( demod_ask == 2 )) = 0;
ask_results(2, find( demod_ask == 2 )) = 1;
ask_results(1, find( demod_ask == 3 )) = 1;
ask_results(2, find( demod_ask == 3 )) = 1;
correct_bits_ask = sum(sum(ask_results == bin));
BER_ask = (N_bits - correct_bits_ask)/N_bits

% 4-opsk
psk_results = zeros(2,length(demod_psk));
psk_results(1, find( demod_psk == 0 )) = 0;
psk_results(2, find( demod_psk == 0 )) = 0;
psk_results(1, find( demod_psk == 1 )) = 1;
psk_results(2, find( demod_psk == 1 )) = 0;
psk_results(1, find( demod_psk == 2 )) = 0;
psk_results(2, find( demod_psk == 2 )) = 1;
psk_results(1, find( demod_psk == 3 )) = 1;
psk_results(2, find( demod_psk == 3 )) = 1;
correct_bits_psk = sum(sum(psk_results == bin));
BER_psk = (N_bits - correct_bits_psk)/N_bits

% 4-qam
qam_results = zeros(2,length(demod_qam));
qam_results(1, find( demod_qam == 0 )) = 0;
qam_results(2, find( demod_qam == 0 )) = 0;
qam_results(1, find( demod_qam == 1 )) = 1;
qam_results(2, find( demod_qam == 1 )) = 0;
qam_results(1, find( demod_qam == 2 )) = 0;
qam_results(2, find( demod_qam == 2 )) = 1;
qam_results(1, find( demod_qam == 3 )) = 1;
qam_results(2, find( demod_qam == 3 )) = 1;
correct_bits_qam = sum(sum(qam_results == bin));
BER_qam = (N_bits - correct_bits_qam)/N_bits
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% received plots vs constellation
figure(1)
subplot(3,1,1)
hold on; title('4-ASK rx'); xlabel('I'); ylabel('Q');
scatter(real(rx_ask_wavefm), imag(rx_ask_wavefm),'b');
scatter(real(ampl_ask), imag(ampl_ask),'r', 'LineWidth', 5);
legend('Received','Constellation Map');
hold off;
subplot(3,1,2)
hold on; title('4-OPSK rx'); xlabel('I'); ylabel('Q');
scatter(real(rx_psk_wavefm), imag(rx_psk_wavefm),'b');
scatter(real(ampl_psk), imag(ampl_psk),'r', 'LineWidth', 5);
hold off;
subplot(3,1,3)
hold on; title('4-QAM rx'); xlabel('I'); ylabel('Q');
scatter(real(rx_qam_wavefm), imag(rx_qam_wavefm),'b');
scatter(real(ampl_qam), imag(ampl_qam),'r', 'LineWidth', 5);
hold off;

% downsampled vs constellation
figure(2)
subplot(3,1,1)
hold on; title('4-ASK downsampled'); xlabel('I'); ylabel('Q');
scatter(real(downsampled_ask), imag(downsampled_ask),'b');
scatter(real(ampl_ask), imag(ampl_ask),'r', 'LineWidth', 5);
legend('pre-demodulation','Constellation Map');
hold off;
subplot(3,1,2)
hold on; title('4-OPSK downsampled'); xlabel('I'); ylabel('Q');
scatter(real(downsampled_psk), imag(downsampled_psk),'b');
scatter(real(ampl_psk), imag(ampl_psk),'r', 'LineWidth', 5);
hold off;
subplot(3,1,3)
hold on; title('4-QAM downsampled'); xlabel('I'); ylabel('Q');
scatter(real(downsampled_qam), imag(downsampled_qam),'b');
scatter(real(ampl_qam), imag(ampl_qam),'r', 'LineWidth', 5);
hold off;





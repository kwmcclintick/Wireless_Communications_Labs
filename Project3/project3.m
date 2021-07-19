% AUTHOR: Kyle McClintick
clear;
clc;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Design the Phase Locked Loop here
Downsampling = 32;%downsampling factor...
Upsampling = 64;
K = 1;
A = 1/sqrt(Upsampling/Downsampling);
PhaseRecoveryLoopBandwidth = 0.01; % Normalized loop bandwidth for fine frequency compensation
PhaseRecoveryDampingFactor = 1; % Damping Factor for fine frequency compensation
PostFilterOversampling = Upsampling/Downsampling;
PhaseErrorDetectorGain = 2*K*A^2+2*K*A^2; % K_p for Fine Frequency Compensation PLL, determined by 2KA^2 (for binary PAM), QPSK could be treated as two individual binary PAM
PhaseRecoveryGain = 1; % K_0 for Fine Frequency Compensation PLL

theta = PhaseRecoveryLoopBandwidth/(PhaseRecoveryDampingFactor + ...
                0.25/PhaseRecoveryDampingFactor)/PostFilterOversampling;
d = 1 + 2*PhaseRecoveryDampingFactor*theta+theta*theta;
K1 = (4*PhaseRecoveryDampingFactor*theta/d)/...
                (PhaseErrorDetectorGain*PhaseRecoveryGain);
K2 = (4*theta*theta/d)/...
                (PhaseErrorDetectorGain*PhaseRecoveryGain);

FineFreqCompensator = QPSKFineFrequencyCompensator(...
                'ProportionalGain', K1, ...
                'IntegratorGain', K2, ...
                'DigitalSynthesizerGain', -1*PhaseRecoveryGain);% Fine Frequency Correction..
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call Phy Layer (DO NOT EDIT THIS!!!)
rx(Upsampling, Downsampling, A, K, FineFreqCompensator)










% QPSK receiver for over the air transmission..
%% Parameter initialization..
warning('off');
clear;
clc;
agc = comm.AGC;%automatic gain control
Downsampling = 32;%downsampling factor...
Upsampling = 64;
M = 4;%modulation order (QPSK)
Fs = 100e3;% 
Ts = 1/Fs;
FrameSize = 100;% Number of Modulated Symbols
MessageLength = 105;
BarkerLength = 13;
RCFiltSpan = 10;
Rolloff = 0.5;
RxBufferedFrames = 10;
% Square root raised cosine receive filter
rctFiltRx = comm.RaisedCosineReceiveFilter(...
  'Shape',                  'Square root', ...
  'RolloffFactor',          Rolloff, ...
  'FilterSpanInSymbols',    RCFiltSpan, ...
  'InputSamplesPerSymbol', Upsampling,...
  'DecimationFactor',Downsampling);
            % Visualize the impulse response
            %fvtool(rctFiltRx, 'Analysis', 'impulse');
K = 1;
A = 1/sqrt(Upsampling/Downsampling);
% Look into model for details for details of PLL parameter choice. Refer equation 7.30 of "Digital Communications - A Discrete-Time Approach" by Michael Rice. 
PhaseErrorDetectorGain = 2*K*A^2+2*K*A^2; % K_p for Fine Frequency Compensation PLL, determined by 2KA^2 (for binary PAM), QPSK could be treated as two individual binary PAM
PhaseRecoveryGain = 1; % K_0 for Fine Frequency Compensation PLL
TimingErrorDetectorGain = 2*K*A^2*2.7; % K_p for Timing Recovery PLL, determined by 2KA^2*2.7 (for binary PAM), QPSK could be treated as two individual binary PAM, 2.7 is for raised cosine filter with roll-off factor 0.5
TimingRecoveryGain = -1; % K_0 for Timing Recovery PLL, fixed due to modulo-1 counter structure
CoarseCompFrequencyResolution = 50; % Frequency resolution for coarse frequency compensation
PhaseRecoveryLoopBandwidth = 0.01; % Normalized loop bandwidth for fine frequency compensation
PhaseRecoveryDampingFactor = 1; % Damping Factor for fine frequency compensation
TimingRecoveryLoopBandwidth = 0.01; % Normalized loop bandwidth for timing recovery
TimingRecoveryDampingFactor = 1; % Damping Factor for timing recovery
PostFilterOversampling = Upsampling/Downsampling;
fc = 450e6;
Gain = 59;
Decimationfactor = 100e6/Fs;
USRPFs = 1/Fs;
USRPFrameLength = Upsampling*FrameSize*RxBufferedFrames;
FrameTime = USRPFrameLength/Fs;
StopTime = 30;
%% Setting up the system objects...
Buffer = dsp.Buffer(FrameSize*2,FrameSize);
bbc = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];%bipolar barker code
ubc = ((bbc+1)/2)';%unipolar barker code..
temp = (repmat(ubc,1,2))';%repeating barker sequence for I-Q components...
Header = temp(:);
bpsk = comm.QPSKModulator('PhaseOffset',pi/4,'BitInput',true);
ModulatedHeader = step(bpsk,Header);




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qpskdemod = comm.QPSKDemodulator('PhaseOffset',pi/4,'BitOutput',true);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




biterr = comm.ErrorRate;
Correlator = dsp.Crosscorrelator;
cdiagcorruptsignal = comm.ConstellationDiagram('SamplesPerSymbol',Upsampling,...
    'XLimits',[-1 1],'YLimits',[-1 1],'SymbolsToDisplaySource',...
    'Property','SymbolsToDisplay',60000,'Title','Corrupted Signal');
cdiagfreqcorr = comm.ConstellationDiagram('SamplesPerSymbol',Upsampling/Downsampling,...
    'XLimits',[-1 1],'YLimits',[-1 1],'SymbolsToDisplaySource',...
    'Property','SymbolsToDisplay',6000,'Title','Frequency Correction');
cdiagfreqcorr.Position(1) = cdiagcorruptsignal.Position(1)+100;
cdiagtimecorr =  comm.ConstellationDiagram('SamplesPerSymbol',1,...
    'XLimits',[-0.06 0.06],'YLimits',[-0.06 0.06],'SymbolsToDisplaySource',...
    'Property');
cdiagtimecorr.Position(1) = cdiagfreqcorr.Position(1)+500;
coarseFreqcomp = QPSKCoarseFrequencyCompensator('ModulationOrder',M,...
	'CoarseCompFrequencyResolution', CoarseCompFrequencyResolution,...
	'SampleRate', Fs, 'DownsamplingFactor', Downsampling);%coarse Frequency Correction...
%% Starting the receiver chain...
theta = PhaseRecoveryLoopBandwidth/(PhaseRecoveryDampingFactor + ...
                0.25/PhaseRecoveryDampingFactor)/PostFilterOversampling;
d = 1 + 2*PhaseRecoveryDampingFactor*theta+theta*theta;
K1 = (4*PhaseRecoveryDampingFactor*theta/d)/...
                (PhaseErrorDetectorGain*PhaseRecoveryGain);
K2 = (4*theta*theta/d)/...
                (PhaseErrorDetectorGain*PhaseRecoveryGain);
            OldOutput = complex(0); % used to store past value
FineFreqCompensator = QPSKFineFrequencyCompensator(...
                'ProportionalGain', K1, ...
                'IntegratorGain', K2, ...
                'DigitalSynthesizerGain', -1*PhaseRecoveryGain);% Fine Frequency Correction..
% Refer C.57 to C.61 in Michael Rice's "Digital Communications 
% - A Discrete-Time Approach" for K1 and K2
theta = TimingRecoveryLoopBandwidth/...
                (TimingRecoveryDampingFactor + ...
                0.25/TimingRecoveryDampingFactor)/PostFilterOversampling;
d = 1 + 2*TimingRecoveryDampingFactor*theta + theta*theta;
K1 = (4*TimingRecoveryDampingFactor*theta/d)/...
                (TimingErrorDetectorGain*TimingRecoveryGain);
K2 = (4*theta*theta/d)/...
                (TimingErrorDetectorGain*TimingRecoveryGain);
TimingRec = QPSKTimingRecovery('ProportionalGain', K1,...
                'IntegratorGain', K2, ...
                'PostFilterOversampling', PostFilterOversampling, ...
                'BufferSize', FrameSize);% Timing Recovery...
      
%% Setting up the radio now..
% radio = comm.SDRuReceiver('N210rx','Gain',Gain,...
%     'CenterFrequency',fc,...
%     'BasebandSampleRate',Fs,...
%     'IPAddress', '192.168.10.2',...
%     'OutputDataType','double','SamplesPerFrame',USRPFrameLength);

currentTime = 0;
len = uint32(0);
temp = [];
[Count,Delay,Phase] = deal(0);
FrameIndex = 0;
SyncIndex = 0;
SyncFlag = true;
testing = [];
TotalData = [];
FrameCount = 0;

saved_signal = load('savedata.mat');
for j=1:length(saved_signal)
    corruptSignal = saved_signal(j,:);
	%Keep accesssing the SDRu system object output until it is valid...
	while len <= 0
		[corruptSignal, len] = step(radio);
    end
	if len > 0
        % Apply automatic gain control to the signal...
        AGCSignal = (1/sqrt(Upsampling))*step(agc,corruptSignal);
        % Pass the signal through square root raised cosine received
        % filter
        FiltSignal = step(rctFiltRx,AGCSignal);
        % Coarsely compensate for the frequency offset..
        coarseCompSignal = step(coarseFreqcomp,FiltSignal);
        coarseCompBuffer = complex(zeros(size(coarseCompSignal))); 
        timingRecBuffer = complex(zeros(size(coarseCompSignal)));
%         cdiagcorruptsignal(FiltSignal);
        for i=1:length(coarseCompSignal)
            % Scalar processing for fine frequency compensation and timing
            % recovery
            fineCompSignal = step(FineFreqCompensator,[OldOutput coarseCompSignal(i)]);
            coarseCompBuffer(i) = fineCompSignal;
            OldOutput = fineCompSignal;
%             cdiagfreqcorr(fineCompSignal);
            % Timing recovery of the received signal...
            [dataOut, isDataValid, timingRecBuffer(i)] = step(TimingRec, fineCompSignal);
%             cdiagtimecorr(dataOut);
            if isDataValid
                %Decoding the received data...
                rxData = step(Buffer,dataOut);
                % Get a frame of data aligned on the frame boundary
                Data = rxData(Delay+1:Delay+length(rxData)/2);
                % Phase Estimation..
                y = mean(conj(ModulatedHeader).*Data(1:BarkerLength));

                %compensating for the phase offset...
                if Data(1)~=0
                    ShiftedData = Data.*exp(-1j*Phase);
                else
                    ShiftedData = complex(zeros(size(Data)));
                end
                
                
                
                
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                % Demodulate the phase recovered data
                demodOut = step(qpskdemod,ShiftedData);
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                
                
                %receivedBits = vitdec(demodOut,trellis,tbl,'cont','hard');
                % Recover the message from the data..
                Received = demodOut(BarkerLength*log2(4)+1:FrameSize*log2(4));
                
                % Turbo decode the demodulated signal. Because the bit mapping from the
                % demodulator is opposite that expected by the turbo decoder, the
                % decoder input must use the inverse of demodulated signal.
               
                receivedmsg = Received(1:MessageLength);
                %Finding delay to achieve frame synchronization...
                z = abs(step(Correlator, ModulatedHeader, dataOut));
                [~,ind] = max(z);
                Delay = mod(length(dataOut)-ind,(length(dataOut)-1));
                %phase ambiguity correction...
                Phase = round(angle(y)*2/pi)/2*pi;
                % Print received frame and estimate the received frame index..
                msgascii = bin2text(receivedmsg');
                disp(msgascii);
%                 FrameCount = FrameCount+1;
%                 TotalData = [TotalData;receivedmsg];
            end
        end
    end
%     currentTime = currentTime+FrameTime;
    len = uint32(0);
end










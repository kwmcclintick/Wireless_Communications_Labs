% AUTHOR: Travis Collins and Kuldeep Gill with edits by Kyle McClintick
% for the purpose of being a ECE3311 project
% DATE: 11/7/18
% Description: Applies AGC to the incoming signal (collected by a PLUTO
% SDR at a 450 MHz carrier, 70 kHz sample rate), then a matched filter
% (Raised Cosine), 

function rx(Upsampling, Downsampling, A, K, FineFreqCompensator)
    %% Parameter initialization..
    warning('off');

    agc = comm.AGC;%automatic gain control
    M = 4;%modulation order (QPSK)
    Fs = 70e3;% 
    FrameSize = 100;% Number of Modulated Symbols
    MessageLength = 105;
    BarkerLength = 13;
    RCFiltSpan = 10;
    Rolloff = 0.5;

    % Square root raised cosine receive filter
    rctFiltRx = comm.RaisedCosineReceiveFilter(...
      'Shape',                  'Square root', ...
      'RolloffFactor',          Rolloff, ...
      'FilterSpanInSymbols',    RCFiltSpan, ...
      'InputSamplesPerSymbol', Upsampling,...
      'DecimationFactor',Downsampling);
                % Visualize the impulse response
                %fvtool(rctFiltRx, 'Analysis', 'impulse');

    % Look into model for details for details of PLL parameter choice. Refer equation 7.30 of "Digital Communications - A Discrete-Time Approach" by Michael Rice. 
    TimingErrorDetectorGain = 2*K*A^2*2.7; % K_p for Timing Recovery PLL, determined by 2KA^2*2.7 (for binary PAM), QPSK could be treated as two individual binary PAM, 2.7 is for raised cosine filter with roll-off factor 0.5
    TimingRecoveryGain = -1; % K_0 for Timing Recovery PLL, fixed due to modulo-1 counter structure
    CoarseCompFrequencyResolution = 50; % Frequency resolution for coarse frequency compensation
    TimingRecoveryLoopBandwidth = 0.01; % Normalized loop bandwidth for timing recovery
    TimingRecoveryDampingFactor = 1; % Damping Factor for timing recovery

    %% Setting up the system objects...
    Buffer = dsp.Buffer(FrameSize*2,FrameSize);
    bbc = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];%bipolar barker code
    ubc = ((bbc+1)/2)';%unipolar barker code..
    temp = (repmat(ubc,1,2))';%repeating barker sequence for I-Q components...
    Header = temp(:);
    bpsk = comm.QPSKModulator('PhaseOffset',pi/4,'BitInput',true);
    ModulatedHeader = step(bpsk,Header);


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Instead of demodulating from the carrier frequency to the baseband,
    % demodulate and filter in three steps
    qpskdemod = comm.QPSKDemodulator('PhaseOffset',pi/4,'BitOutput',true);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    Correlator = dsp.Crosscorrelator;
    coarseFreqcomp = QPSKCoarseFrequencyCompensator('ModulationOrder',M,...
        'CoarseCompFrequencyResolution', CoarseCompFrequencyResolution,...
        'SampleRate', Fs, 'DownsamplingFactor', Downsampling);%coarse Frequency Correction...

    %% Starting the receiver chain...
    % Refer C.57 to C.61 in Michael Rice's "Digital Communications 
    % - A Discrete-Time Approach" for K1 and K2
    PostFilterOversampling = Upsampling/Downsampling;
    theta = TimingRecoveryLoopBandwidth/...
                    (TimingRecoveryDampingFactor + ...
                    0.25/TimingRecoveryDampingFactor)/PostFilterOversampling;
    d = 1 + 2*TimingRecoveryDampingFactor*theta + theta*theta;
    ProportionalGain = (4*TimingRecoveryDampingFactor*theta/d)/...
                    (TimingErrorDetectorGain*TimingRecoveryGain);
    IntegratorGain = (4*theta*theta/d)/...
                    (TimingErrorDetectorGain*TimingRecoveryGain);
    PostFilterOversampling = Upsampling/Downsampling;
    TimingRec = QPSKTimingRecovery('ProportionalGain', ProportionalGain,...
                    'IntegratorGain', IntegratorGain, ...
                    'PostFilterOversampling', PostFilterOversampling, ...
                    'BufferSize', FrameSize);% Timing Recovery...

    len = uint32(0);
    temp = [];
    [Count,Delay,Phase] = deal(0);
    OldOutput = complex(0); % used to store past value

    PLUTO_signal = load('savedata'); % import PLUTO transmission
    for z=1:length(PLUTO_signal)
        corruptSignal = PLUTO_signal.savedata(z,:)';

        % Apply automatic gain control to the signal...
        AGCSignal = (1/sqrt(Upsampling))*step(agc,corruptSignal);
        % Pass the signal through square root raised cosine received
        % filter
        FiltSignal = step(rctFiltRx,AGCSignal);
        % Coarsely compensate for the frequency offset..
        coarseCompSignal = step(coarseFreqcomp,FiltSignal);

        % initialize
        coarseCompBuffer = complex(zeros(size(coarseCompSignal))); 
        timingRecBuffer = complex(zeros(size(coarseCompSignal)));
        for i=1:length(coarseCompSignal)
            % Scalar processing for fine frequency compensation and timing
            % recovery

            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % Fine frequency correction occurs here (PLL)
            fineCompSignal = step(FineFreqCompensator,[OldOutput coarseCompSignal(i)]);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            



            coarseCompBuffer(i) = fineCompSignal;
            OldOutput = fineCompSignal;
            % Timing recovery of the received signal...
            [dataOut, isDataValid, timingRecBuffer(i)] = step(TimingRec, fineCompSignal);

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

                          
                % Demodulate
                demodOut = step(qpskdemod,ShiftedData);



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
            end
        end
        len = uint32(0);
    end
    
    
    
     figure(1)
    scatter(real(coarseCompSignal), imag(coarseCompSignal));
    title('Coarse Frequency Corrected IQ Plot');
    xlim([-2 2]); ylim([-2 2]);
    xlabel('In-Phase'); ylabel('Quadriture');



    figure(2)
    scatter(real(dataOut), imag(dataOut));
    title('PLL Corrected IQ Plot');
    lim = max(abs(fineCompSignal))
    xlim([-2 2]); ylim([-2 2]);
    xlabel('In-Phase'); ylabel('Quadriture');

    figure(3)
    scatter(real(ShiftedData), imag(ShiftedData));
    title('Phase Corrected QPSK IQ Plot');
    xlim([-2 2]); ylim([-2 2]);
    xlabel('In-Phase'); ylabel('Quadriture');
    
    
end




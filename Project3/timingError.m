
%% General system details
sampleRateHz = 1e6; % Sample rate
samplesPerSymbol = 2;
frameSize = 2^10;
numFrames = 100;
numSamples = numFrames*frameSize; % Samples to simulate
modulationOrder = 2;
filterUpsample = 2;
filterSymbolSpan = 4;

%% Visuals
cdPre = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
    'Name','Baseband');
cdPost = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
    'SymbolsToDisplaySource','Property',...
    'SymbolsToDisplay',frameSize/2,...
    'Name','Baseband with Timing Offset');
cdPre.Position(1) = 50;
cdPost.Position(1) = cdPre.Position(1)+cdPre.Position(3)+10;% Place side by side

%% Impairments
snr = 15;
timingOffset = filterUpsample*0.1; % Samples

%% Generate symbols
data = randi([0 modulationOrder-1], numSamples*2, 1);
%mod = comm.QPSKModulator();
mod = comm.BPSKModulator();
modulatedData = mod.step(data);

%% Add TX/RX Filters
TxFlt = comm.RaisedCosineTransmitFilter(...
    'OutputSamplesPerSymbol', filterUpsample,...
    'FilterSpanInSymbols', filterSymbolSpan);

RxFlt = comm.RaisedCosineReceiveFilter(...
    'InputSamplesPerSymbol', filterUpsample,...
    'FilterSpanInSymbols', filterSymbolSpan,...
    'DecimationFactor', filterUpsample/2);

%% Add noise source
chan = comm.AWGNChannel( ...
    'NoiseMethod',  'Signal to noise ratio (SNR)', ...
    'SNR',          snr, ...
    'SignalPower',  1, ...
    'RandomStream', 'mt19937ar with seed');

%% Add delay
varDelay = dsp.VariableFractionalDelay;

%% Setup visualization object(s)
sa = dsp.SpectrumAnalyzer('SampleRate',sampleRateHz,'ShowLegend',true);

%% Model of error
% Add timing offset to baseband signal

filteredData = [];%zeros(length(modulatedData)*2,1);

for k=1:frameSize:(numSamples - frameSize)
    
    timeIndex = (k:k+frameSize-1).';
    
    % Filter signal
    filteredTXData = step(TxFlt, modulatedData(timeIndex));
    
    % Pass through channel
    noisyData = step(chan, filteredTXData);
    
    % Time delay signal
    offsetData = step(varDelay, noisyData, timingOffset); % Variable delay
    
    % Filter signal
    %filteredData(timeIndex) = step(RxFlt, offsetData);
    nFD = step(RxFlt, offsetData);
    filteredData = [filteredData; nFD]; %#ok<AGROW>
    
    % Visualize Error
    step(cdPre,noisyData);step(cdPost,nFD);pause(0.01); %#ok<*UNRCH>
end

%% Relate configuration from TX side
clear mod
inputDataFull = filteredData;
InputFrameLen = frameSize;

%% Setup
SamplesPerSymbol = 2;
SPS = SamplesPerSymbol;
MaxOutputFrameLen = ceil(InputFrameLen*11/SPS/10);
alpha = 0.5;

b = @(Mu) [alpha*Mu^2 - alpha*Mu;...
    -alpha*Mu^2 - (1-alpha)*Mu + 1;...
    -alpha*Mu^2 + (1+alpha)*Mu;...
    alpha*Mu^2 - alpha*Mu];

LoopFilterState = 0;
LoopPreviousInput = 0;
Trigger = false;
TriggerHistory = false(1, SPS);
Mu = 0;
M1Counter = 0;
InterpFilterState = complex(zeros(3, 1),zeros(3, 1));
TEDBuffer = complex(zeros(1, SPS), zeros(1, SPS));
maxOutputSize = ceil(InputFrameLen*11/double(SamplesPerSymbol)/10);
SymbolHolder = complex(zeros(maxOutputSize, 1), zeros(maxOutputSize, 1));

%% Visuals
cdPre = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
    'Name','Baseband');
cdPost = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
    'SymbolsToDisplaySource','Property',...
    'SymbolsToDisplay',frameSize/2,...
    'Name','Baseband with Timing Offset');
cdBuiltin = comm.ConstellationDiagram('ReferenceConstellation', [-1 1],...
    'SymbolsToDisplaySource','Property',...
    'SymbolsToDisplay',frameSize/2,...
    'Name','Baseband with Timing Offset');

%% Define loop gains
DetectorGain = 2.7;
zeta = 1/sqrt(2);
NormalizedLoopBandwidth = 0.1;

% Calculate
Kp = DetectorGain;
K0 = -1;
theta = NormalizedLoopBandwidth/SPS/(zeta + 0.25/zeta);
d  = (1 + 2*zeta*theta + theta^2) * K0 * Kp;
ProportionalGain = (4*zeta*theta) /d;
IntegratorGain   = (4*theta*theta)/d;

%% Setup interpolator
filt = dsp.FIRFilter;
filt.NumeratorSource='Input port';
timeCorrect = zeros(size(inputDataFull));
index = 1;
buffer = [0 0 0];

LoopFilter = dsp.IIRFilter( ...
    'Structure', 'Direct form II transposed', ...
    'Numerator', [1 0], 'Denominator', [1 -1]);

%% Run Zero Crossing Implementation
NumTrigger = 0;
for k=1:frameSize:(numSamples - frameSize)

    timeIndex = (k:k+frameSize-1).';
    inputData = inputDataFull(timeIndex); % Grab frame

    %% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Call Synchronizer (sample by sample on frame)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    for i = 1 : InputFrameLen % Process input frame sample-by-sample

        %% Using Trigger count as output input count
        NumTrigger = NumTrigger + Trigger;

        %% Interpolator
        % Piecewise parabolic interpolator

        % Update coefficents based on new Mu
        intOut = step(filt,inputData(i),fliplr(b(Mu).'));

        if Trigger % Interpolation output as symbols
            SymbolHolder(NumTrigger) = intOut;
        end

        %% Timing error detector (TED)
        % For zero-crossing, TED wants to make every other sample time
        % aligned with a zero crossing in the eye diagram.  Therefore, for
        % zero-cross to operate we must have 2 SPS.

        % TED calculation occurs on a Trigger
        if Trigger && all(~TriggerHistory(2:end))

            % The above condition allows TED update after a skip. If we want
            % TED update to happen only at regular strobings, need to check
            % TriggerHistory(1) == true in addition to the condition above.

            % Calculate the midsample from interpolator output
            t1 = TEDBuffer(end/2 + 1 - rem(SPS,2));
            t2 = TEDBuffer(end/2 + 1);
            midSample = (t1+t2)/2;

            % Rice Notation
            % x -> real(in)
            % y -> imag(in)
            % a0 -> sign(x)
            % a1 -> sign(y)

            % Data aided method
            e = real(midSample) * (sign(real(TEDBuffer(1))) - sign(real(intOut))) + ...
                imag(midSample) * (sign(imag(TEDBuffer(1))) - sign(imag(intOut)));
        else
            e = 0;
        end

        % Handle bit stuffing/skipping
        switch sum([TriggerHistory(2:end), Trigger])
            case 0
                % No Trigger in awhile
            case 1
                % Update buffer (Normal)
                TEDBuffer = [TEDBuffer(2:end), intOut];
            otherwise % > 1
                % Stuff condition
                TEDBuffer = [TEDBuffer(3:end), 0, intOut];
        end

        %% Loop filter
        % The output v essentially adjust the period of the Trigger

        %loopFiltOut = LoopPreviousInput + LoopFilterState;
        %v = e*ProportionalGain + loopFiltOut;
        %LoopFilterState = loopFiltOut;
        %LoopPreviousInput = e*IntegratorGain;

        v = LoopFilter(e*IntegratorGain);
        v = e*ProportionalGain + v;

        %% Interpolation controller
        % We want to get a Trigger signal every SPS samples, which should
        % happen when we acquire lock
        % Essentially a Trigger should occur at the start of a symbol,
        % therefore based on the start of the symbol we will then utilize
        % the next sample (which should also be within that symbol),
        % interpolate across both with the appropriate delay and produce an
        % output
        %
        % Modulo-1 counter interpolation controller which generates/updates
        % Trigger signal and fractional interpolation interval.
        %
        % Trigger - Trigger edge found and to output data
        % Mu - fractional interpolation interval

        W = v + 1/SPS; % W should be small or == SPS when locked

        TriggerHistory = [TriggerHistory(2:end), Trigger];

        % Determine if we have an underflow, aka we are at the start edge
        % of a sample
        Trigger = (M1Counter < W);

        % With a new underflow we must update interpolator coefficients
        if Trigger
            Mu = M1Counter / W;
        end

        % Update counter
        M1Counter = mod(M1Counter - W, 1);

    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Output 1 SPS
    y = SymbolHolder(1:NumTrigger, 1);
    timeCorrect(index:index+length(y)-1) = y;
    index = index + length(y);

end

%% Remove zeros
timeCorrect = timeCorrect(1:index-1);

%% Visualize
for k=1:frameSize:(index - frameSize)
    
    timeIndex = (k:k+frameSize-1).';
    step(cdPre,timeCorrect(timeIndex));pause(0.01);
end




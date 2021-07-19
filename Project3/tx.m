% QPSK Transmitter with barker sequence..
%% Parameter Initialization...
Fs = 200e3;%sampling frequency
Upsampling = 64;%Upsampling factor..
RCFiltSpan = 10;% Filter span for 10 symbols
Rolloff = 0.5;
fc = 450e6;
Gain = 0;%Maximum Gain allowed...
%% Bits Generation with Barker Sequence..
message = 'Hello World 000';
bbc = [+1 +1 +1 +1 +1 -1 -1 +1 +1 -1 +1 -1 +1];%bipolar barker code
ubc = ((bbc+1)/2)';%unipolar barker code..
temp = (repmat(ubc,1,2))';%repeating barker sequence for I-Q components...
Header = temp(:);
data = [Header;text2bin(message)';zeros(1,69)'];

%% Modulating the data
bpsk = comm.QPSKModulator('BitInput',true,'PhaseOffset',pi/4);
rctFiltTx = comm.RaisedCosineTransmitFilter(...
  'Shape',                  'Square root', ...
  'RolloffFactor',          Rolloff, ...
  'FilterSpanInSymbols',    RCFiltSpan, ...
  'OutputSamplesPerSymbol', Upsampling);
% Visualize the impulse response

% fvtool(rctFiltTx, 'Analysis', 'impulse')
modulateddata = bpsk(data);
transmittedsignal = rctFiltTx(modulateddata);
% Visualization..
cdiag = comm.ConstellationDiagram('Name','Transmitted Signal','SamplesPerSymbol',1,...
    'ShowGrid',true,'ShowReferenceConstellation',false,...
    'XLimits',[-2 2],'YLimits',[-2 2]);
cdiag(transmittedsignal);%BPSK constellation..
%% Transmitting packets over the air..
radio  = comm.SDRuTransmitter('N210tx','Gain',Gain,...
    'CenterFrequency',fc,...
    'BasebandSampleRate',Fs,...
    'IPAddress', '192.168.10.2');
currenttime = 0;
% Starting the transmission...
while true
   radio(transmittedsignal);
   disp('Frame Transmitted');
%    currenttime = currenttime+Frametime;
end




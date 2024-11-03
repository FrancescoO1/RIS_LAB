function [risElementCoeff,w,Htot] = calculateRISCoeff(enableRIS,risCh,carrier)
    
    % Initialize RIS element coefficients (gain = 1 and random phase)
    numRISElements = prod(risCh.RISSize);
    theta = 2*pi*rand(1,numRISElements); % uniformly distributed phases in [0, 2*pi]
    risElementCoeff = exp(1i*theta);
    % If the RIS algorithm is disabled, this is the value used in the
    % simulation

    [G,h] = getRISChannelFreqResponse(risCh,carrier);

    H = h*diag(risElementCoeff)*G;

    % Calculate precoding weights using MRT (maximum ratio transmission)
    w = H'/norm(H);

    if enableRIS
        numIter = 10;  % Number of iterations
        Htot = zeros(numIter+1,1); % Total channel H*wf. Used to check convergence
        Htot(1) = H*w;
        for n = 1:numIter
            % Need to calculate B = h*diag(G*w), where w is the transmitter
            % precoding vector, G is the transmitter to RIS channel matrix and
            % h is the RIS to receiver channel matrix. B is used to calculate
            % the new RIS phase values theta as -angle(B)
            B = h*diag(G*w);

            % Calculate the new phase vector phi that compensates for the phase changes in hr, G and w
            theta = -angle(B);

            % New RIS coefficients
            risElementCoeff = exp(1i*theta); 

            % New combined channel matrix H
            H = h*diag(risElementCoeff)*G;            

            % Get new weights w based on new channel matrix H
            w = H'/norm(H);

            % Overall channel response for this iteration. Used to check
            % convergence
            Htot(n+1) = H*w;
        end
    else
        Htot = H*w;
    end
end

function [G,h] = getRISChannelFreqResponse(risCh,carrier)
    % Calculate the overall RIS channel response averaging the 
    % channel response, that is, one channel matrix for the whole
    % bandwidth, and not per resource element. 
    [TxRISGrid,RISRxGrid] = channelResponse(risCh,carrier);
    
    h = zeros(size(RISRxGrid,3),size(RISRxGrid,4));
    RISRxGrid = mean(RISRxGrid,[1 2]); % assume flat in time and freq
    h(:,:) = RISRxGrid(1,1,:,:);

    G = zeros(size(TxRISGrid,3),size(TxRISGrid,4));
    TxRISGrid = mean(TxRISGrid,[1 2]); % assume flat in time and freq
    G(:,:) = TxRISGrid(1,1,:,:);
end

function [txWaveform,pdschSymbols,pdschIndices,dmrsSymbols,dmrsIndices,waveformInfo,pdschGrid] ...
    = generateTxWaveform(carrier,pdsch,wtx,txPower)

    nTxAnts = length(wtx);

    % PDSCH and PDSCH DM-RS
    [pdschIndices,pdschInfo] = hpre6GPDSCHIndices(carrier,pdsch);
    pdschBits = randi([0 1],pdschInfo.G,1);
    pdschSymbols = hpre6GPDSCH(carrier,pdsch,pdschBits);
    
    dmrsSymbols = hpre6GPDSCHDMRS(carrier,pdsch);
    dmrsIndices = hpre6GPDSCHDMRSIndices(carrier,pdsch);
    
    % PDSCH precoding
    pdschSymbolsPrecoded = pdschSymbols*wtx;
    
    % Grid
    pdschGrid = hpre6GResourceGrid(carrier,nTxAnts);
    
    [~,pdschAntIndices] = nrExtractResources(pdschIndices,pdschGrid);
    pdschGrid(pdschAntIndices) = pdschSymbolsPrecoded;
    
    % PDSCH DM-RS precoding and mapping
    for p = 1:size(dmrsSymbols,2)
        [~,dmrsAntIndices] = nrExtractResources(dmrsIndices(:,p),pdschGrid);
        pdschGrid(dmrsAntIndices) = pdschGrid(dmrsAntIndices) + dmrsSymbols(:,p)*wtx(p,:);
    end
    
    % OFDM Modulation
    [txWaveform,waveformInfo] = hpre6GOFDMModulate(carrier,pdschGrid);

    % Scale power of transmitted signal to desired value.
    txWaveform = scalePower(txWaveform,txPower);
end

function  waveform = scalePower(waveform,desiredPower)
    % Scale input WAVEFORM to achieve the desiredPower in dBm
    K = sqrt((1/mean(rms(waveform).^2))*10^((desiredPower-30)/10));
    waveform = K*waveform;
end

function info = hpre6GOFDMInfo(carrier,varargin)

    narginchk(1,7);

    % Validate carrier input
    mustBeA(carrier,'pre6GCarrierConfig');

    % Parse options
    fcnName = 'hpre6GOFDMInfo';
    optNames = {'Nfft','Windowing','CarrierFrequency'};
    opts = nr5g.internal.parseOptions(fcnName,optNames,varargin{:});

    % Get OFDM information
    info = nrOFDMInfo(carrier,opts);

end

rng('default')

fc = 28e9; % Carrier frequency (Hz)
carrier = pre6GCarrierConfig("NSizeGrid",50,"SubcarrierSpacing",120);
pdsch = pre6GPDSCHConfig("Modulation","16QAM","NumLayers",1,"PRBSet",0:carrier.NSizeGrid-1);

ofdmInfo = hpre6GOFDMInfo(carrier);
fs = ofdmInfo.SampleRate;

% Transmit array
txArraySize = [4 2 2 1 1]; % [M N P Mg Ng]
Ntx =  prod(txArraySize(1:3));

% RIS parameters 
ris.Enable = true;
ris.Size = [8 4 2];
numRISElements = prod(ris.Size);

% Wavelength
c = physconst("lightspeed");
lambda = c/fc; 

% RIS element size
ris.dx = lambda/5;
ris.dy = lambda/5;

ris.A = 0.8; % Amplitude of reflection coefficients of each RIS element

risCh = hpre6GRISChannel("SampleRate",ofdmInfo.SampleRate,"RISSize",ris.Size,"CarrierFrequency",fc);
risCh.TransmitAntennaArray.Size = txArraySize;
risCh.ReceiveAntennaArray.Size = [1 1 1 1 1]; % Receive array must be a single antenna in this example
risCh.Seed = randi([0,(2^32)-1]);             % Random seed

% Calculate the overall channel delay
chInfo = risCh.info;
channelDelay = chInfo.TxRISChInfo.MaximumChannelDelay + chInfo.RISRXChInfo.MaximumChannelDelay;

txPower = 20;                                % dBm. Average transmit power per antenna
% Noise and interference parameters
noiseFigure = 7;                             % dB
thermalNoiseDensity = -174;                  % dBm/Hz
rxInterfDensity = -165.7;                    % dBm/Hz

% Calculate the corresponding noise power
totalNoiseDensity = 10*log10(10^((noiseFigure+thermalNoiseDensity)/10)+10^(rxInterfDensity/10));
BW = 12*carrier.NSizeGrid*carrier.SubcarrierSpacing*1e3;
noisePower = totalNoiseDensity+10*log10(BW); % dBm
N = 10^((noisePower-30)/10);                 % Linear

d_max = 50;                          % Distanza massima dal RIS
grid_size = 20;                        % Numero di punti della griglia (grid_size x grid_size)
x = linspace(-d_max, d_max, grid_size);
y = linspace(-d_max, d_max, grid_size);
[X, Y] = meshgrid(x, y);               % Creazione della griglia
Z = zeros(grid_size);                  % Matrice per memorizzare la capacità di Shannon

% Distanza tra RIS e Tx
dt = 105; % distanza fissa tra RIS e Tx

for i = 1:grid_size
    for j = 1:grid_size
        % Distanza tra RIS e Rx
        dr = sqrt(X(i,j)^2 + Y(i,j)^2);

        % Path loss (power gain)
        PL = ((4*pi*dt*dr)/(prod(ris.Size(1:3))*ris.dx*ris.dy*ris.A)).^2;

        [risElementCoeff,w,overallChannel] = calculateRISCoeff(ris.Enable,risCh,carrier);

        w = w.'; % transmitter requires w to be of the form number of layers by number of transmit antennas

        [txWaveform,pdschSym,pdschInd,dmrsSym,dmrsInd,~,pdschGrid] = generateTxWaveform(carrier,pdsch,w,txPower);

        % Pad with zeros to flush full slot from the channel
        txWaveform = [txWaveform; zeros(channelDelay,size(txWaveform,2))]; 

        % Send signal through RIS assisted channel
        rxWaveformNoisless = risCh(txWaveform,risElementCoeff);

        % Apply path loss
        rxWaveformNoisless = (1/sqrt(PL))* rxWaveformNoisless;

        % Add noise
        noise = sqrt(N)*randn(size(rxWaveformNoisless),"like",rxWaveformNoisless);

        % SNR
        SNR = 10*log10(mean(rms(rxWaveformNoisless).^2)/mean(rms(noise).^2));

        % Capacità di Shannon (Hz*bit/s)
        Z(i,j) = BW * log2(1 + 10^(SNR / 10));  % Conversione da dB in scala lineare
    end
end

y = flip(y);
figure;
h = heatmap(x, y, Z);
h.Colormap = jet;
h.ColorbarVisible = 'on';
xlabel('Distanza X (m)');
ylabel('Distanza Y (m)');
title('Capacità di Shannon nel territorio');
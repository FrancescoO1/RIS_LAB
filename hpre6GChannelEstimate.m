%hpre6GChannelEstimate Practical channel estimation
%   [H,NVAR,INFO] = hpre6GChannelEstimate(...) performs channel estimation,
%   returning channel estimate H, noise variance estimate NVAR and
%   information structure INFO. H is a K-by-N-by-R-by-P array where K is
%   the number of subcarriers, N is the number of OFDM symbols, R is the
%   number of receive antennas and P is the number of reference signal
%   ports. NVAR is scalar indicating the measured variance of additive
%   white Gaussian noise on the received reference symbols. INFO is a
%   structure containing the field:
%   AveragingWindow - a 2-element row vector [F T] indicating the number 
%                     of adjacent reference symbols in the frequency
%                     direction F and time direction T over which averaging
%                     was performed prior to interpolation
%
%   [H,NVAR,INFO] = hpre6GChannelEstimate(CARRIER,RXGRID,REFIND,REFSYM)
%   performs channel estimation on the received resource grid RXGRID using
%   reference symbols REFSYM whose locations are given by REFIND.
%
%   CARRIER is an extended carrier configuration object, <a 
%   href="matlab:help('pre6GCarrierConfig')"
%   >pre6GCarrierConfig</a>.
%   Only this object property is relevant for this function:
%
%   CyclicPrefix      - Cyclic prefix ('normal', 'extended')
%
%   RXGRID is an array of size K-by-L-by-R. K is the number of subcarriers,
%   given by CARRIER.NSizeGrid * 12. L is the number of OFDM symbols in one
%   slot, given by CARRIER.SymbolsPerSlot.
%
%   REFIND and REFSYM are the reference signal indices and symbols,
%   respectively. REFIND is an array of 1-based linear indices addressing a
%   K-by-L-by-P resource array. P is the number of reference signal ports
%   and is inferred from the range of values in REFIND. Only nonzero
%   elements in REFSYM are considered. Any zero-valued elements in REFSYM
%   and their associated indices in REFIND are ignored.
%
%   [H,NVAR,INFO] = hpre6GChannelEstimate(...,NAME,VALUE,...) specifies
%   additional options as NAME,VALUE pairs:
%
%   'CDMLengths'      - A 2-element row vector [FD TD] specifying the 
%                       length of FD-CDM and TD-CDM despreading to perform.
%                       A value of 1 for an element indicates no CDM and a
%                       value greater than 1 indicates the length of the
%                       CDM. For example, [2 1] indicates FD-CDM2 and no
%                       TD-CDM. The default is [1 1] (no orthogonal
%                       despreading)
%
%   'AveragingWindow' - A 2-element row vector [F T] specifying the number
%                       of adjacent reference symbols in the frequency
%                       domain F and time domain T over which to average
%                       prior to interpolation. F and T must be odd or
%                       zero. If F or T is zero, the averaging value is
%                       determined automatically from the estimated SNR
%                       (calculated using NVAR). The default is [0 0]
%
%   Example:
%   % Create a resource grid containing the PDSCH DM-RS and pass it through
%   % a TDL-C channel. Estimate the channel response and compare it with
%   % the perfect channel estimator.
%
%   carrier = pre6GCarrierConfig;
%   carrier.NSizeGrid = 330;
%   carrier.SubcarrierSpacing = 120;
%   pdsch = pre6GPDSCHConfig;
%   pdsch.PRBSet = 0:carrier.NSizeGrid-1;
%   dmrsInd = hpre6GPDSCHDMRSIndices(carrier,pdsch);
%   dmrsSym = hpre6GPDSCHDMRS(carrier,pdsch);
%   nTxAnts = pdsch.NumLayers;
%   txGrid = hpre6GResourceGrid(carrier,nTxAnts);
%   txGrid(dmrsInd) = dmrsSym;
%
%   [txWaveform,ofdmInfo] = hpre6GOFDMModulate(carrier,txGrid);
%
%   channel = nrTDLChannel;
%   channel.NumTransmitAntennas = nTxAnts;
%   channel.NumReceiveAntennas = 1;
%   channel.SampleRate = ofdmInfo.SampleRate;
%   channel.DelayProfile = 'TDL-C';
%   channel.DelaySpread = 100e-9;
%   channel.MaximumDopplerShift = 20;
%   chInfo = info(channel);
%   maxChDelay = chInfo.MaximumChannelDelay;
%   [rxWaveform,pathGains,sampleTimes] = channel([txWaveform; zeros(maxChDelay,nTxAnts)]);
%   
%   refGrid = hpre6GResourceGrid(carrier,nTxAnts);
%   refGrid(dmrsInd) = dmrsSym;
%   offset = hpre6GTimingEstimate(carrier,rxWaveform,refGrid);
%   rxWaveform = rxWaveform(1+offset:end,:);
%
%   rxGrid = hpre6GOFDMDemodulate(carrier,rxWaveform);
%
%   [H,nVar,estInfo] = hpre6GChannelEstimate(carrier,rxGrid,dmrsInd,dmrsSym);
%
%   pathFilters = getPathFilters(channel);
%   H_ideal = hpre6GPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);
%
%   figure;
%   subplot(1,2,1);
%   imagesc(abs(H));
%   xlabel('OFDM symbol');
%   ylabel('Subcarrier');
%   title('Practical estimate magnitude');
%   subplot(1,2,2);
%   imagesc(abs(H_ideal));
%   xlabel('OFDM symbol');
%   ylabel('Subcarrier');
%   title('Perfect estimate magnitude');
%
%   See also hpre6GTimingEstimate, hpre6GPerfectChannelEstimate, 
%   pre6GCarrierConfig.

%   Copyright 2023 The MathWorks, Inc.

function [H,nVar,info] = hpre6GChannelEstimate(carrier,rxGrid,refInd,refSym,varargin)

    narginchk(4,8);

    % Validate carrier input
    mustBeA(carrier,'pre6GCarrierConfig');

    % Parse options
    fcnName = 'hpre6GChannelEstimate';
    optNames = {'CDMLengths','AveragingWindow'};
    opts = nr5g.internal.parseOptions(fcnName,optNames,varargin{:});

    % Get channel estimate, noise estimate, and information structure
    [H,nVar,info] = nrChannelEstimate(carrier,rxGrid,refInd,refSym,opts);

end

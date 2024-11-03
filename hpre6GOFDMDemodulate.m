%hpre6GOFDMDemodulate OFDM demodulation
%   GRID = hpre6GOFDMDemodulate(CARRIER,WAVEFORM) performs OFDM
%   demodulation of the time domain waveform, WAVEFORM, given extended
%   carrier configuration object CARRIER.
%
%   CARRIER is an extended carrier configuration object, <a 
%   href="matlab:help('pre6GCarrierConfig')"
%   >pre6GCarrierConfig</a>.
%   Only these object properties are relevant for this function:
%
%   SubcarrierSpacing - Subcarrier spacing in kHz (default 15)
%   CyclicPrefix      - Cyclic prefix ('normal', 'extended')
%   NSizeGrid         - Number of resource blocks in carrier resource
%                       grid (default 52)
%   NSlot             - Slot number
%
%   WAVEFORM is a T-by-R matrix where T is the number of time-domain
%   samples in the waveform and R is the number of receive antennas. Only
%   complete OFDM symbols in WAVEFORM are demodulated. Any additional
%   samples corresponding to part of an OFDM symbol are discarded.
%
%   GRID is an array of size K-by-L-by-R where K is the number of
%   subcarriers and L is the number of OFDM symbols.
%
%   GRID = hpre6GOFDMDemodulate(...,NAME,VALUE) specifies additional
%   options as NAME,VALUE pairs to allow control over the OFDM
%   demodulation:
%
%   Nfft                 - Desired number of FFT points to use in the OFDM
%                          demodulator. If absent or set to [], a default 
%                          value is selected based on other parameters, see
%                          <a href="matlab: doc('nrOFDMDemodulate')"
%                          >nrOFDMDemodulate</a> for details
%   CarrierFrequency     - Carrier frequency (in Hz) to calculate the phase
%                          decompensation applied for each OFDM symbol 
%                          (denoted f_0 in TS 38.211 Section 5.4). Default 
%                          is 0
%   CyclicPrefixFraction - Starting position of OFDM symbol demodulation 
%                          (FFT window position) within the cyclic prefix.
%                          Specified as a fraction of the cyclic prefix, in
%                          the range [0,1], with 0 representing the start
%                          of the cyclic prefix and 1 representing the end
%                          of the cyclic prefix. Default is 0.5
%
%   Note that for the numerologies specified in TS 38.211 Section 4.2, 
%   extended cyclic prefix length is only applicable for 60 kHz subcarrier
%   spacing.
%
%   Example:
%   % Configure carrier for 330 resource blocks and 120 kHz subcarrier 
%   % spacing, consistent with a 500 MHz bandwidth
%   carrier = pre6GCarrierConfig;
%   carrier.NSizeGrid = 330;
%   carrier.SubcarrierSpacing = 120;
%
%   % Configure PDSCH and create PDSCH DM-RS symbols and indices
%   pdsch = pre6GPDSCHConfig;
%   pdsch.NumLayers = 2;
%   sym = hpre6GPDSCHDMRS(carrier,pdsch);
%   ind = hpre6GPDSCHDMRSIndices(carrier,pdsch);
%
%   % Create a carrier resource grid and map PDSCH DM-RS symbols
%   txGrid = hpre6GResourceGrid(carrier,pdsch.NumLayers);
%   txGrid(ind) = sym;
%
%   % Perform OFDM modulation
%   [txWaveform,info] = hpre6GOFDMModulate(carrier,txGrid);
%
%   % Apply a simple 2-by-1 channel
%   H = [0.6; 0.4];
%   rxWaveform = txWaveform * H;
%
%   % Perform OFDM demodulation
%   rxGrid = hpre6GOFDMDemodulate(carrier,rxWaveform);
%
%   See also pre6GCarrierConfig, hpre6GOFDMInfo, hpre6GOFDMModulate,
%   hpre6GResourceGrid.

%   Copyright 2023 The MathWorks, Inc.

function grid = hpre6GOFDMDemodulate(carrier,waveform,varargin)

    narginchk(2,8);

    % Validate carrier input
    mustBeA(carrier,'pre6GCarrierConfig');

    % Parse options
    fcnName = 'hpre6GOFDMDemodulate';
    optNames = {'CarrierFrequency','Nfft','CyclicPrefixFraction'};
    opts = nr5g.internal.parseOptions(fcnName,optNames,varargin{:});
    opts.SampleRate = [];

    % Perform OFDM demodulation
    grid = nrOFDMDemodulate(carrier,waveform,opts);

end

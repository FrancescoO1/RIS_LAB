%hpre6GPerfectChannelEstimate perfect channel estimation
%   H = hpre6GPerfectChannelEstimate(...) performs perfect channel
%   estimation, producing a perfect channel estimate H, by reconstructing
%   the channel impulse response from information about the propagation
%   channel and then performing OFDM demodulation. H is a
%   K-by-N-by-Nr-by-Nt array where K is the number of subcarriers, N is the
%   number of OFDM symbols, Nr is the number of receive antennas and Nt is
%   the number of transmit antennas.
%
%   H = hpre6GPerfectChannelEstimate(CARRIER,PATHGAINS,PATHFILTERS,OFFSET,SAMPLETIMES)
%   performs perfect channel estimation by reconstructing the channel
%   impulse response from the channel path gains array PATHGAINS and path
%   filter impulse response matrix PATHFILTERS, and then performing OFDM
%   demodulation according to the extended carrier configuration given by
%   CARRIER. OFFSET specifies the timing offset and SAMPLETIMES specifies
%   the sample times of the channel snapshots.
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
%   NStartGrid        - Start of carrier resource grid relative to CRB 0
%                       (default 0)
%   NSlot             - Slot number
%
%   PATHGAINS must be an array of size Ncs-by-Np-by-Nt-by-Nr, where Ncs is
%   the number of channel snapshots and Np is the number of paths. The
%   times of the channel snapshots are given by the SAMPLETIMES input (see
%   below).
%
%   PATHFILTERS must be a matrix of size Nh-by-Np where Nh is the number of
%   impulse response samples.
%
%   OFFSET specifies the timing offset, an integer number of samples
%   indicating where the OFDM demodulation will start on the reconstructed
%   waveform.
%
%   SAMPLETIMES specifies the sample times of the channel snapshots.
%   SAMPLETIMES must be of size Ncs-by-1 and specifies the time of
%   occurrence of each channel snapshot (the 1st dimension of PATHGAINS).
%   Ensure that the channel snapshots span at least one slot. The function
%   performs channel estimation for each complete slot.
%   H = hpre6GPerfectChannelEstimate(...,NAME,VALUE) specifies additional
%   options as NAME,VALUE pairs to allow control over the OFDM demodulation
%   of the channel impulse responses:
%
%   Nfft                 - Desired number of FFT points to use in the OFDM
%                          demodulator. If absent or set to [], a default 
%                          value is selected based on other parameters, see
%                          <a href="matlab: doc('nrOFDMDemodulate')"
%                          >nrOFDMDemodulate</a> for details
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
%   % Configure a TDL-C channel with 100 ns delay spread and plot the 
%   % estimated channel magnitude response for the first receive antenna.
%
%   carrier = pre6GCarrierConfig;
%   carrier.NSizeGrid = 330;
%   carrier.SubcarrierSpacing = 120;
%   ofdmInfo = hpre6GOFDMInfo(carrier);
%   
%   tdl = nrTDLChannel;
%   tdl.DelayProfile = 'TDL-C';
%   tdl.DelaySpread = 100e-9;
%   tdl.MaximumDopplerShift = 300;
%   tdl.SampleRate = ofdmInfo.SampleRate;
%   
%   T = tdl.SampleRate * 1e-3;
%   tdlInfo = info(tdl);
%   Nt = tdlInfo.NumTransmitAntennas;
%   in = complex(randn(T,Nt),randn(T,Nt));
%
%   [~,pathGains,sampleTimes] = tdl(in);
%   pathFilters = getPathFilters(tdl);
%
%   offset = nrPerfectTimingEstimate(pathGains,pathFilters)
%   hest = hpre6GPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes);
%   size(hest)
%
%   figure;
%   surf(abs(hest(:,:,1)));
%   shading('flat');
%   xlabel('OFDM symbols');
%   ylabel('Subcarriers');
%   zlabel('|H|');
%   title('Channel magnitude response');
%
%   See also hpre6GChannelEstimate, hpre6GTimingEstimate,
%   pre6GCarrierConfig.

%   Copyright 2023 The MathWorks, Inc.

function H = hpre6GPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes,varargin)

    narginchk(5,9);

    % Validate carrier input
    mustBeA(carrier,'pre6GCarrierConfig');

    % Parse options
    fcnName = 'hpre6GPerfectChannelEstimate';
    optNames = {'Nfft','CyclicPrefixFraction'};
    opts = nr5g.internal.parseOptions(fcnName,optNames,varargin{:});
    opts.SampleRate = [];

    % Perform perfect channel estimation
    H = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes,opts);

    % For large subcarrier spacings, the longer cyclic prefix which occurs 
    % every 0.5 ms can be longer than the FFT used in OFDM demodulation. 
    % This requires special handling because of how the perfect channel
    % estimator utilizes OFDM demodulation
    ofdminfo = nr5g.internal.OFDMInfo(carrier,opts);
    N = size(H,2);
    cpLengths = nr5g.internal.OFDMInfoRelativeNSlot(ofdminfo,carrier.NSlot,N);
    cpLengths = cpLengths(1:N);
    cpFFTs = floor(cpLengths/ofdminfo.Nfft);
    idx = find(cpFFTs > 1,1,'first');
    if (~isempty(idx))
        delta_t = (cpLengths(idx) + ofdminfo.Nfft) / ofdminfo.SampleRate;
        H2 = nrPerfectChannelEstimate(carrier,pathGains,pathFilters,offset,sampleTimes + delta_t,opts);
        H(:,idx,:,:) = H2(:,idx+1,:,:);
    end

end

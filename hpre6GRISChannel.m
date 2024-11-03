classdef hpre6GRISChannel < matlab.System
%hpre6GRISChannel RIS (Reflective Intelligent Surface) channel with CDL
%   CHAN = hpre6GRISChannel creates a RIS CDL channel. This channel
%   includes:
%   * a CDL channel to model the propagation between the transmitter
%   and the RIS
%   * the RIS, modeled by an amplitude and phase changed applied by each
%   RIS element
%   * a second CDL channel to model the propagation between the RIS and the
%   receiver.
%
%   CHAN = hpre6GRISChannel(Name,Value) creates a RIS CDL channel object,
%   CHAN, with the specified property Name set to the specified Value. You
%   can specify additional name-value pair arguments in any order as
%   (Name1,Value1,...,NameN,ValueN).
%   
%   Step method syntax for ChannelFiltering set to true:
%
%   Y = step(CHAN,X,COEFFS) filters the input signal X through a RIS CDL
%   MIMO fading channel, applies the RIS element coefficients COEFFS, and
%   returns the result in Y. The input X can be a double or single
%   precision data type scalar, vector, or 2-D matrix. X is of size
%   Ns-by-Nt, where Ns is the number of samples and Nt is the number of
%   transmit antennas. COEFFS contains the RIS element coefficients and is
%   a row vector with Nris values, where Nris is the number of RIS
%   elements. Y is the output signal of size Ns-by-Nr, where Nr is the
%   number of receive antennas. Y contains values of the same type as the
%   input signal X.
% 
%   System objects may be called directly like a function instead of using
%   the step method. For example, y = step(obj,x) and y = obj(x) are
%   equivalent.
%
%   hpre6GRISChannel methods:
%
%   step                   - Filter input signal through a RIS CDL MIMO fading channel
%                            (see above)
%   release                - Allow property value and input characteristics changes
%   clone                  - Create a RIS CDL channel object with same property values
%   isLocked               - Locked status (logical)
%   <a href="matlab:help hpre6GRISChannel/reset">reset</a>                  - Reset states of filters, and random stream if the
%                            RandomStream property is set to 'mt19937ar with seed'
%   <a href="matlab:help hpre6GRISChannel/infoImpl">info</a>                   - Return characteristic information about the RIS 
%                            channel
%   channelResponse        - Return the channel response in the frequency
%                            domain for both CDL channels
%
%   hpre6GRISChannel properties:
%
%   DelayProfileTxRIS         - Delay profile of CDL channel between transmitter and RIS
%   DelayProfileRISRx         - Delay profile of CDL channel between RIS and receiver
%   DelaySpreadTxRIS          - Desired delay spread (s) of transmitter to RIS channel
%   DelaySpreadRISRx          - Desired delay spread (s) of RIS to receiver channel
%   SampleRate                - Input signal sample rate (Hz)
%   TransmitAntennaArray      - Transmit antenna array characteristics
%   TransmitArrayOrientation  - Orientation of the transmit antenna array
%   RISSize                   - RIS size 
%   RISOrientation            - Orientation of the RIS
%   ReceiveAntennaArray       - Receive antenna array characteristics
%   ReceiveArrayOrientation   - Orientation of the receive antenna array
%   MaximumDopplerShiftTxRIS  - Maximum Doppler shift (Hz) of transmitter to RIS channel
%   MaximumDopplerShiftRISRx  - Maximum Doppler shift (Hz) of RIS to receiver channel
%   CarrierFrequency          - Carrier frequency (Hz)
%   Seed                      - Initial seed of mt19937ar random number stream

%   Copyright 2023 The MathWorks, Inc.

    properties (Dependent, Nontunable)

        %DelayProfileTxRIS Delay profile of CDL channel between transmitter and RIS
        %   Specify the CDL delay profile as one of 'CDL-A', 'CDL-B',
        %   'CDL-C', 'CDL-D' or 'CDL-E'.
        %
        %   The default value of this property is 'CDL-A'.
        DelayProfileTxRIS

        %DelayProfileRISRx Delay profile of CDL channel between RIS and receiver
        %   Specify the CDL delay profile as one of 'CDL-A', 'CDL-B',
        %   'CDL-C', 'CDL-D' or 'CDL-E'.
        %
        %   The default value of this property is 'CDL-A'.
        DelayProfileRISRx


        %SampleRate Input signal sample rate (Hz)
        %   Specify the sample rate of the input signal in Hz as a double
        %   precision, real, positive scalar.
        %
        %   The default value of this property is 15.36e6 Hz.
        SampleRate

        %TransmitAntennaArray Transmit antenna array characteristics
        %   Structure or Phased Array System Toolbox antenna array object
        %   specifying the transmit antenna array. The default is a
        %   structure containing the following fields:
        %   Size                - Size of antenna array [M,N,P,Mg,Ng]. M
        %                         and N are the number of rows and columns
        %                         in the antenna array. P is the number of
        %                         polarizations (1 or 2). Mg and Ng are the
        %                         number of row and column array panels
        %                         respectively. The defaults are 
        %                         [2,2,2,1,1].
        %   ElementSpacing      - Element spacing in wavelengths expressed
        %                         as [lambda_v lambda_h dg_v dg_h]
        %                         representing the vertical and horizontal
        %                         element spacing and the vertical and
        %                         horizontal panel spacing respectively.
        %                         The panel spacing is measured from the
        %                         center of the panels. The defaults are
        %                         [0.5 0.5 1.0 1.0].
        %   PolarizationAngles  - Polarization angles [theta rho] in
        %                         degrees applicable when P = 2. The
        %                         defaults are [45 -45] degrees.
        %   Element             - Antenna element radiation pattern. One of
        %                         'isotropic' or '38.901' (see TR 38.901 
        %                         Section 7.3). The default value is 
        %                         '38.901'.
        %   PolarizationModel   - Model describing how to determine the
        %                         radiation field patterns based on a 
        %                         defined radiation power pattern (see
        %                         TR 38.901 Section 7.3.2). One of
        %                         'Model-1' or 'Model-2'. The default value
        %                         is 'Model-2'.
        %
        % The antenna array elements are mapped to the input waveform
        % channels (columns) in the order that a 5-D array of size
        % M-by-N-by-P-by-Mg-by-Ng is linearly indexed (across the
        % dimensions first to last). The size of the array is given by
        % TransmitAntennaArray.Size = [M,N,P,Mg,Ng]. For example, an
        % antenna array of size [4,8,2,2,2] has the first M = 4 channels
        % mapped to the first column of the first polarization angle of the
        % first panel. The next M (equals 4) antennas are mapped to the
        % next column and so on, such that the first M*N (equals 32)
        % channels are mapped to the first polarization angle of the
        % complete first panel. Then the next 32 channels are mapped in the
        % same fashion to the second polarization angle for the first
        % panel. Subsequent sets of M*N*P (equals 64) channels are then
        % mapped to the remaining panels, panel rows first then panel
        % columns.
        TransmitAntennaArray

        %TransmitArrayOrientation Orientation of the transmit antenna array
        % Mechanical orientation of the transmit antenna array [alpha; beta;
        % gamma] in degrees (bearing, downtilt, slant). The default values
        % [0; 0; 0] indicate that the broadside direction of the array
        % points to the positive x-axis.
        TransmitArrayOrientation

        %ReceiveAntennaArray Receive antenna array characteristics
        %   Structure or Phased Array System Toolbox antenna array object
        %   specifying the receive antenna array. The default is a
        %   structure containing the following fields:
        %   Size                - Size of antenna array [M,N,P,Mg,Ng]. M
        %                         and N are the number of rows and columns
        %                         in the antenna array. P is the number of
        %                         polarizations (1 or 2). Mg and Ng are the
        %                         number of row and column array panels
        %                         respectively. The defaults are 
        %                         [1,1,2,1,1].
        %   ElementSpacing      - Element spacing in wavelengths expressed
        %                         as [lambda_v lambda_h dg_v dg_h]
        %                         representing the vertical and horizontal
        %                         element spacing and the vertical and
        %                         horizontal panel spacing respectively.
        %                         The panel spacing is measured from the
        %                         center of the panels. The defaults are
        %                         [0.5 0.5 0.5 0.5].
        %   PolarizationAngles  - Polarization angles [theta rho] in
        %                         degrees applicable when P is set to 2.
        %                         The defaults are [0 90] degrees.
        %   Element             - Antenna element radiation pattern. One of
        %                         'isotropic' or '38.901' (see TR 38.901
        %                         Section 7.3). The default value is 
        %                         'isotropic'.
        %   PolarizationModel   - Model describing how to determine the
        %                         radiation field patterns based on a 
        %                         defined radiation power pattern (see
        %                         TR 38.901 Section 7.3.2). One of
        %                         'Model-1' or 'Model-2'. The default value
        %                         is 'Model-2'.
        %
        % The antenna array elements are mapped to the output waveform
        % channels (columns) in the order that a 5-D array of size
        % M-by-N-by-P-by-Mg-by-Ng is linearly indexed (across the
        % dimensions first to last). The size of the array is given by
        % ReceiveAntennaArray.Size = [M,N,P,Mg,Ng]. For example, an antenna
        % array of size [4,8,2,2,2] has the first M (equals 4) channels
        % mapped to the first column of the first polarization angle of the
        % first panel. The next M (equals 4) antennas are mapped to the
        % next column and so on, such that the first M*N (equals 32)
        % channels are mapped to the first polarization angle of the
        % complete first panel. Then the next 32 channels are mapped in the
        % same fashion to the second polarization angle for the first
        % panel. Subsequent sets of M*N*P (equals 64) channels are then
        % mapped to the remaining panels, panel rows first then panel
        % columns.
        ReceiveAntennaArray

        %ReceiveArrayOrientation Orientation of the receive antenna array
        % Mechanical orientation of the receive antenna array [alpha; beta;
        % gamma] in degrees (bearing, downtilt, slant). The default values
        % [0; 0; 0] indicate that the broadside direction of the array
        % points to the positive x-axis.
        ReceiveArrayOrientation

        %RISSize RIS size
        %   Size of RIS array [M,N,P]. M and N are the number of rows
        %   and columns in the RIS array. P is the number of polarizations
        %   (1 or 2). The defaults are [8 4 2]. Elements in the RIS are
        %   spaced by half lambda. Each element in the RIS has an isotropic
        %   radiation characteristic.
        RISSize

        %RISOrientation Orientation of the RIS
        % Mechanical orientation of the RIS [alpha; beta; gamma] in degrees
        % (bearing, downtilt, slant). The default values [0; 0; 0] indicate
        % that the broadside direction of the array points to the positive
        % x-axis.
        RISOrientation

        %Seed Initial seed of mt19937ar random number stream
        %   Specify the initial seed of a mt19937ar random number generator
        %   algorithm as a double precision, real, nonnegative integer
        %   scalar. The Seed reinitializes the mt19937ar random number
        %   stream in the reset method. Each of the two CDL channels in the
        %   RIS channel has an mt19937ar random number stream. The stream
        %   corresponding to the channel between the transmitter and the
        %   RIS is initialized with the provided seed value. The random
        %   number stream corresponding to the channel between the RIS and
        %   the receiver is initialized with the value Seed+1.
        %
        %   The default value of this property is 73. 
        Seed

        %CarrierFrequency Carrier frequency (Hz)
        %   Specify the carrier frequency in Hertz as a scalar.
        %   
        %   The default value of this property is 4 GHz.
        CarrierFrequency

        %DelaySpreadTxRIS Desired delay spread (s) of transmitter to RIS channel
        %   Specify the desired RMS delay spread in seconds as a scalar.
        %   This applies to the transmitter to RIS channel.
        %
        %   The default value of this property is 1e-12 s.
        DelaySpreadTxRIS

        %DelaySpreadRISRx Desired delay spread (s) of RIS to receiver 
        %   Specify the desired RMS delay spread in seconds as a scalar.
        %   This applies to the RIS to receiver channel.
        %
        %   The default value of this property is 1e-12 s.
        DelaySpreadRISRx

        %MaximumDopplerShiftTxRIS Maximum Doppler shift (Hz) between transmitter and RIS
        %   Specify the maximum Doppler shift between the transmitter and
        %   the RIS. This applies to all channel paths in Hertz as a double
        %   precision, real, nonnegative scalar.
        % 
        %   When you set the MaximumDopplerShiftTxRIS to 0, the CDL channel
        %   between transmitter and RIS remains static for the entire
        %   input. To generate a new channel realization, call the reset
        %   function on the hpre6GRISChannel object.
        %
        %   The default value of this property is 0 Hz. 
        MaximumDopplerShiftTxRIS

        %MaximumDopplerShiftRISRx Maximum Doppler shift (Hz) between RIS and receiver
        %   Specify the maximum Doppler shift between the RIS and the
        %   receiver. This applies to all channel paths in Hertz as a
        %   double precision, real, nonnegative scalar.
        % 
        %   When you set the MaximumDopplerShiftRISRx to 0, the CDL channel
        %   between RIS and receiver remains static for the entire input.
        %   To generate a new channel realization, call the reset function
        %   on the hpre6GRISChannel object.
        %
        %   The default value of this property is 5 Hz. 
        MaximumDopplerShiftRISRx
    end

    properties (Access = private)
        Ch1 % CDL channel between transmitter and RIS
        Ch2 % CDL channel between RIS and receiver
    end

    properties (Access = private)
        PathGains1
        PathGains2
        SampleTimes1
        SampleTimes2
        PathFilters1
        PathFilters2
    end

    methods
        % hpre6GRISChannel constructor
        function obj = hpre6GRISChannel(varargin)
            obj.Ch1 = nrCDLChannel;
            obj.Ch1.SampleRate = 15360000;

            % Tx array
            txArrayCh1 = struct(Size=[2 2 2 1 1],ElementSpacing=[0.5 0.5 1 1],PolarizationAngles=[45 -45], ...
                Element='38.901',PolarizationModel='Model-2');
            obj.Ch1.TransmitAntennaArray = txArrayCh1;

            % ch1 rx antenna array is the RIS
            rxArrayCh1 = struct(Size=[8 4 2 1 1],ElementSpacing=[0.5 0.5 1 1],PolarizationAngles=[45 -45], ...
                Element='isotropic',PolarizationModel='Model-2');
            obj.Ch1.ReceiveAntennaArray = rxArrayCh1;

            obj.Ch1.NormalizeChannelOutputs = false;
            obj.Ch1.DelayProfile = 'CDL-A';
            obj.Ch1.DelaySpread = 1e-12; % flat channel
            obj.Ch1.CarrierFrequency = 4e9;
            obj.Ch1.MaximumDopplerShift = 0;
            obj.Ch1.Seed = 73;

            obj.Ch2 = nrCDLChannel;
            obj.Ch2.SampleRate = obj.Ch1.SampleRate;

            % Rx array
            rxArrayCh2 = struct(Size=[1,1,1,1,1],ElementSpacing=[0.5 0.5 0.5 0.5],PolarizationAngles=[0 90], ...
                Element='isotropic',PolarizationModel='Model-2');
            obj.Ch2.ReceiveAntennaArray = rxArrayCh2;

            obj.Ch2.NormalizeChannelOutputs = false;
            % ch2 tx antenna array is the RIS
            obj.Ch2.TransmitAntennaArray = obj.Ch1.ReceiveAntennaArray;
            obj.Ch2.TransmitArrayOrientation = obj.Ch1.ReceiveArrayOrientation;

            obj.Ch2.DelayProfile = 'CDL-A';
            obj.Ch2.DelaySpread = obj.Ch1.DelaySpread;
            obj.Ch2.CarrierFrequency = obj.Ch1.CarrierFrequency;
            obj.Ch2.MaximumDopplerShift = 5;
            obj.Ch2.Seed = obj.Ch1.Seed+1;

            setProperties(obj,nargin,varargin{:});
        end

        function [h1,h2,offset1,offset2] = channelResponse(obj,carrier)
        %channelResponse Get channel frequency response 
        %   [HTXRIS,HRISRX,OFFSETTXRIS,OFFSETRISRX] =
        %   channelResponse(obj,CARRIER) returns the frequency domain
        %   channel responses HTxRIS and HRISRx for both CDL channels for
        %   the parameters in the |pre6GCarrierConfig| object CARRIER.
        %   HTXRIS is a K-by-N-by-Nris-by-Nt array where K is the number of
        %   subcarriers, N is the number of OFDM symbols, Nris is the
        %   number of RIS elements, and Nt is the number of transmit
        %   antennas. HTXRIS is the frequency domain response of the CDL
        %   channel between the transmitter and the RIS.
        %
        %   HRISRX is a K-by-N-by-Nr-by-Nris array where Nr is the number
        %   of receive antennas. HRISRX is the frequency domain response of
        %   the CDL channel between the RIS and the receiver.
        %
        %   OFFSETTXRIS and OFFSETRISRX are the timing offsets, in samples,
        %   of the transmitter to RIS, and RIS to receiver channels,
        %   respectively. These offsets indicate the number of samples by
        %   which the peak of the channel impulse response is delayed.

            % In case data has not been sent through the channel
            if isempty(obj.PathGains1) || isempty(obj.PathGains2) || isempty(obj.SampleTimes1) || isempty(obj.SampleTimes2)
                Ntx = prod(obj.TransmitAntennaArray.Size(1:3));
                numRISElements = prod(obj.RISSize(1:3));
                x = zeros(1,Ntx,'like',single(1i));
                risElementCoeff = ones(1,numRISElements)/numRISElements; % normalize diagonal matrix
                % If data has not been sent through the channel, the path
                % gains do not exist. Pass one sample into the channel to
                % generate the path gains.
                [~] = obj.stepImpl(x,risElementCoeff);
                obj.releaseImpl(); % here we assume single input data, step may be called with double input data
            end

            if isempty(obj.PathFilters1) || isempty(obj.PathFilters2)
                [obj.PathFilters1,obj.PathFilters2] = getPathFilters(obj);
            end
            
            offset1 = nrPerfectTimingEstimate(obj.PathGains1,obj.PathFilters1);
            offset2 = nrPerfectTimingEstimate(obj.PathGains2,obj.PathFilters2);
            TxRISGrid = hpre6GPerfectChannelEstimate(carrier,obj.PathGains1,obj.PathFilters1,offset1,obj.SampleTimes1);
            RISRxGrid = hpre6GPerfectChannelEstimate(carrier,obj.PathGains2,obj.PathFilters2,offset2,obj.SampleTimes2);

            h1 = TxRISGrid;
            h2 = RISRxGrid;
        end
        
        function set.DelayProfileTxRIS(obj,val)
            obj.Ch1.DelayProfile = val;
        end

        function delayProfileTxRIS= get.DelayProfileTxRIS(obj)
            delayProfileTxRIS = obj.Ch1.DelayProfile;
        end

        function set.DelayProfileRISRx(obj,val)
            obj.Ch2.DelayProfile = val;
        end

        function delayProfileRISRx= get.DelayProfileRISRx(obj)
            delayProfileRISRx = obj.Ch2.DelayProfile;
        end

        function set.CarrierFrequency(obj,val)
            obj.Ch1.CarrierFrequency = val;
            obj.Ch2.CarrierFrequency = val;
        end

        function carrierFreq = get.CarrierFrequency(obj)
            carrierFreq = obj.Ch1.CarrierFrequency;
        end

        function set.MaximumDopplerShiftTxRIS(obj,val)
            obj.Ch1.MaximumDopplerShift = val;
        end

        function maximumDopplerShiftTxRIS = get.MaximumDopplerShiftTxRIS(obj)
            maximumDopplerShiftTxRIS = obj.Ch1.MaximumDopplerShift;
        end

        function set.MaximumDopplerShiftRISRx(obj,val)
            obj.Ch2.MaximumDopplerShift = val;
        end

        function maximumDopplerShiftRISRx = get.MaximumDopplerShiftRISRx(obj)
            maximumDopplerShiftRISRx = obj.Ch2.MaximumDopplerShift;
        end

        function set.DelaySpreadTxRIS(obj,val)
            obj.Ch1.DelaySpread = val;
        end

        function delaySpread = get.DelaySpreadTxRIS(obj)
            delaySpread = obj.Ch1.DelaySpread;
        end

        function set.DelaySpreadRISRx(obj,val)
            obj.Ch2.DelaySpread = val;
        end

        function delaySpread = get.DelaySpreadRISRx(obj)
            delaySpread = obj.Ch2.DelaySpread;
        end

        function set.Seed(obj,val)
            obj.Ch1.Seed = val;
            obj.Ch2.Seed = val+1;
        end

        function s1 = get.Seed(obj)
            s1 = obj.Ch1.Seed;
        end

        function set.SampleRate(obj,val)
            obj.Ch1.SampleRate = val;
            obj.Ch2.SampleRate = val;
        end

        function sampleRate= get.SampleRate(obj)
            sampleRate = obj.Ch1.SampleRate;
        end

        function set.TransmitAntennaArray(obj,val)
            obj.Ch1.TransmitAntennaArray = val;
        end

        function txAntennaArray = get.TransmitAntennaArray(obj)
            txAntennaArray = obj.Ch1.TransmitAntennaArray;
        end

        function set.TransmitArrayOrientation(obj,val)
            obj.Ch1.TransmitArrayOrientation = val;
        end

        function val = get.TransmitArrayOrientation(obj)
            val = obj.Ch1.TransmitArrayOrientation;
        end  

        function set.ReceiveAntennaArray(obj,val)
            obj.Ch2.ReceiveAntennaArray = val;
        end

        function rxAntennaArray = get.ReceiveAntennaArray(obj)
            rxAntennaArray = obj.Ch2.ReceiveAntennaArray;
        end

        function set.ReceiveArrayOrientation(obj,val)
            obj.Ch2.ReceiveArrayOrientation = val;
        end

        function val = get.ReceiveArrayOrientation(obj)
            val = obj.Ch2.ReceiveArrayOrientation;
        end

        function set.RISSize(obj,val)
            arraySize = [val 1 1]; % only one panel needed in the RIS
            obj.Ch1.ReceiveAntennaArray.Size = arraySize;
            obj.Ch2.TransmitAntennaArray.Size = arraySize;
        end

        function risSize = get.RISSize(obj)
            risSize = obj.Ch1.ReceiveAntennaArray.Size(1:3);
        end

        function set.RISOrientation(obj,val)
            obj.Ch1.ReceiveArrayOrientation = val;
            obj.Ch2.TransmitArrayOrientation = val;
        end

        function val = get.RISOrientation(obj)
            val = obj.Ch1.ReceiveArrayOrientation;
        end

    end

    methods (Access = protected)
        function setupImpl(obj)
            % Perform one-time calculations, such as computing constants
            [obj.PathFilters1,obj.PathFilters2] = getPathFilters(obj);
        end

        function y = stepImpl(obj,in,elementCoeffs)
            [risRx,obj.PathGains1,obj.SampleTimes1] = obj.Ch1(in);
            risTx = risRx.*elementCoeffs;
            [y,obj.PathGains2,obj.SampleTimes2] = obj.Ch2(risTx);
        end

        function resetImpl(obj)
            % Initialize / reset discrete-state properties
            obj.Ch1.reset();
            obj.Ch2.reset();
        end

        function releaseImpl(obj)
            release(obj.Ch1);
            release(obj.Ch2);
        end

        % hpre6GRISCHannel infoImpl method
        function s = infoImpl(obj)
        %info Returns characteristic information about the RIS channel
        %   S = info(CHAN) returns a structure containing characteristic
        %   information, S, about the both CDL fading channels within the
        %   RIS channel:
        %   TxRISChInfo         - Information about the CDL channel between
        %                         the transmitter and the RIS.
        %   RISRxChInfo         - Information about the CDL channel between
        %                         the RIS and the receiver.
        %
        %   Each of these structures contains the following fields:
        % 
        %   ClusterTypes        - A row cell array of character vectors,
        %                         indicating the type of each cluster in
        %                         the delay profile ('LOS',
        %                         'SubclusteredNLOS', 'NLOS')
        %   PathDelays          - A row vector providing the delays of the
        %                         discrete channel paths, in seconds. These
        %                         values include the effect of the desired
        %                         delay spread scaling, and desired
        %                         K-factor scaling if enabled. 
        %   AveragePathGains    - A row vector of the average gains of the
        %                         discrete path or cluster, in dB. These
        %                         values include the effect of K-factor
        %                         scaling if enabled.
        %   AnglesAoD           - A row vector of the Azimuth of Departure
        %                         angles of the clusters in degrees.
        %   AnglesAoA           - A row vector of the Azimuth of Arrival
        %                         angles of the clusters in degrees.
        %   AnglesZoD           - A row vector of the Zenith of Departure
        %                         angles of the clusters in degrees.
        %   AnglesZoA           - A row vector of the Zenith of Arrival
        %                         angles of the clusters in degrees.
        %   KFactorFirstCluster - K-factor of first cluster of delay
        %                         profile, in dB. If the first cluster of
        %                         the delay profile follows a Laplacian
        %                         rather than Rician distribution,
        %                         KFactorFirstCluster is -Inf.
        %   ClusterAngleSpreads - A row vector of cluster-wise RMS angle
        %                         spreads [C_ASD C_ASA C_ZSD C_ZSA] (deg)
        %   XPR                 - A scalar of cross polarization power
        %                         ratio (dB) or NaN when DelayProfile =
        %                         'Custom' and XPR is specified as a
        %                         matrix.
        %   NumTransmitAntennas - Number of transmit antennas.
        %   NumInputSignals     - Number of input signals to CDL channel.
        %   NumReceiveAntennas  - Number of receive antennas.
        %   NumOutputSignals    - Number of output signals of CDL channel.
        %   ChannelFilterDelay  - Channel filter delay in samples.
        %   MaximumChannelDelay - Maximum channel delay in samples. This
        %                         delay consists of the maximum propagation
        %                         delay and the ChannelFilterDelay.

            s.TxRISChInfo = obj.Ch1.info();
            s.RISRXChInfo = obj.Ch2.info();
            
        end

        % nrCDLChannel saveObjectImpl method
        function s = saveObjectImpl(obj)

            % Save public properties
            s = saveObjectImpl@matlab.System(obj);

            % Save private properties
            s.Ch1 = obj.Ch1;
            s.Ch2 = obj.Ch2;
            s.PathGains1 = obj.PathGains1;
            s.PathGains2 = obj.PathGains2;
            s.SampleTimes1 = obj.SampleTimes1;
            s.SampleTimes2 = obj.SampleTimes2;
            s.PathFilters1 = obj.PathFilters1;
            s.PathFilters2 = obj.PathFilters2;
            
        end

        % nrCDLChannel loadObjectImpl method
        function loadObjectImpl(obj,s,wasLocked)

            % Load private properties
            obj.Ch1 = s.Ch1;
            obj.Ch2 = s.Ch2;
            obj.PathGains1 = s.PathGains1;
            obj.PathGains2 = s.PathGains2;
            obj.SampleTimes1 = s.SampleTimes1;
            obj.SampleTimes2 = s.SampleTimes2;
            obj.PathFilters1 = s.PathFilters1;
            obj.PathFilters2 = s.PathFilters2;

            loadObjectImpl@matlab.System(obj,s,wasLocked);                        
            
        end

    end


    methods (Access = private)
        function [pathFiltersTxRIS,pathFiltersRISRx] = getPathFilters(obj)
        %getPathFilters Get path filter impulse responses
        %   [FTxRIS,FRISRX] = getPathFilters(obj) returns two double
        %   precision real matrices of size Nh-by-Np where Nh is the number
        %   of impulse response samples and Np is the number of paths.
        %   Matrices FTxRIS and FRISRX correspond to each of the two CDL
        %   channels of the RIS CDL channel (Tx to RIS, and RIS to RX).
        %   Each column of FTxRIS and FRISRX contains the filter impulse
        %   response for each path of the delay profile. This information
        %   facilitates reconstruction of a perfect channel estimate when
        %   used in conjunction with the PATHGAINS output of the step
        %   method. These filters do not change once the object is created,
        %   therefore it only needs to be called once.
            pathFiltersTxRIS = obj.Ch1.getPathFilters();
            pathFiltersRISRx = obj.Ch2.getPathFilters();
        end
    end
end

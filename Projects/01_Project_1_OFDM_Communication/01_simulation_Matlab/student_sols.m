%NO_PFILE
function [funs, student_id] = student_sols()
%STUDENT_SOLS Contains all student solutions to problems.

% ----------------------------------------
%               STEP 1
% ----------------------------------------
% Set to your birthdate / the birthdate of one member in the group.
% Should a numeric value of format YYYYMMDD, e.g.
% student_id = 19900101;
% This value must be correct in order to generate a valid secret key.
student_id = 19930425;


% ----------------------------------------
%               STEP 2
% ----------------------------------------
% Your task is to implement the following skeleton functions.
% You are free to use any of the utility functions located in the same
% directory as this file.


	function z = add_cyclic_prefix(x,Ncp)  %#ok<*INUSD>
		% Adds (prepends) a Ncp long cyclic prefix to the ofdm block x.
		x = x(:);   %#ok<*NASGU> % Ensure x is a column vector        
        z = [x(end-(Ncp-1):end); x(:)];
    end

    function x = remove_cyclic_prefix(z,Ncp)
        % Removes a Ncp long cyclic prefix from the ofdm package z
        z = z(:);   % Ensure z is a column vector
        x = z(Ncp+1:end);
    end

    function symb = bits2qpsk(bits)
        % Encode bits as qpsk symbols 
        % ARGUMENTS:
        % bits = array of bits. Numerical values converted as:
        %   zero -> zero
        %   nonzero -> one
        % Must be of even length!
        % OUTPUT:
        % symb = complex array of qpsk symbols encoding the bits. Will contain
        % length(bits)/2 elements. Valid output symbols are
        % 1/sqrt(2)*(+/-1 +/- i). Symbols grouped by pairs of bits, where
        % the first corresponds the real part of the symbol while the
        % second corresponds to the imaginary part of the symbol.
        
        % Convert bits vector of +/- 1
        bits = double(bits);
        bits = bits(:);
        bits(bits ~= 0) = 1;
        bits(bits == 0) = -1;

        if rem(length(bits),2) == 1
            error('bits must be of even length');
        end
    
        % s(1) = ... ( b(1)    + i*b(2) )
        % s(2) = ... ( b(3)    + i*b(4) )
        % s(3) = ... ( b(5)    + i*b(6) )
        % s(k) = ... ( b(2k-1) + i*b(2k) )
        symb = zeros(numel(bits)/2,1);
        for k=1: numel(symb)
            symb(k) = 1/sqrt(2)* ( bits(2*k - 1) + 1i * bits(2*k) );
        end
    end

    function bits  = qpsk2bits(x)
        % Convert qpsk symbols to bits.
        % Output will be a vector twice as long as the input x, with values
        % 0 or 1.
        x = x(:);
        bits = false(2*length(x),1);
        % TODO: finish this function
        % Note: you only need to check which quadrant of the complex plane
        % the symbol lies in in order to map it to a pair of bits. The
        % first bit corresponds to the real part of the symbol while the
        % second bit corresponds to the imaginary part of the symbol.
        
        for k=1:length(x)
            bits(2*k-1) = real(x(k)) > 0;
            bits(2*k)   = imag(x(k)) > 0;
        end
        
        % Ensure output is of correct type
        % zero value -> logical zero
        % nonzero value -> logical one
        bits = logical(bits);
    end

    function [rx, evm, ber, symbs] = sim_ofdm_known_channel(tx, h, N_cp, snr, sync_err)
        % Simulate OFDM signal transmission/reception over a known channel.
        %
        % -----------------------------------------------------------------
        % NOTE: THIS FUNCTION WILL NOT BE SELF-TESTED!
        % It will be up to you to study the output from this function and
        % determine if the results are correct or not.
        % -----------------------------------------------------------------
        %
        % Arguments:
        %   tx          Bits to transmit [-]
        %   h           Channel impulse response [-]
        %   N_cp        Cyclic prefix length [samples]
        %   snr         Channel signal/noise ration to apply [dB]
        %   sync_err    Reciever synchronization error [samples]
        % Outputs:
        %   rx          Recieved bits [-]
        %   evm         Error vector magnitude (see below) [-]
        %   ber         Bit error rate (see below) [-]
        %   symbs       Structure containing fields:
        %       .tx         Transmitted symbols
        %       .rx_pe      Recieved symbols, pre-equalization
        %       .rx_e       Recieved symbols, post-equalization
        %
        % In this function, you will now fully implement a simulated
        % base-band OFDM communication scheme. The relevant steps in this
        % are:
        %   - Get a sequence of bits to transmit
        %   - Convert the bits to OFDM symbols
        %   - Create an OFDM block from the OFDM symbols
        %   - Add a cyclic prefix
        %   - Simulate the transmission and reception over the channel using the
        %   simulate_baseband_channe function.
        %   - Remove the cyclic prefix from the recieved message.
        %   - Equalize the recieved symbols by the channel gain
        %   - Convert the equalized symbols back to bits
        %   - Compare the recieved bits/symbols to the transmitted bits/symbols.
        %
        % If you have implemented the skeleton functions earlier in this
        % file then this function will be very simple as you can call your
        % functions to perform the needed tasks.
		
		warning('Note that this function is _not_ self-tested. It is up to you to study the output any verify that it is correct! You can remove this warning if you wish.');
        
        % Ensure inputs are column vectors
        tx = tx(:);
        h = h(:);
        
        % Convert bits to QPSK symbols
        % x = 0; %TODO: This line is missing some code!
        x = bits2qpsk(tx);
        
        symbs.tx = x;   % Store transmitted symbols for later
        
        % Number of symbols in message
        N = length(x);

        plot_FT(x,'s',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);
        
        % Create OFDM time-domain block using IDFT 
        z = ifft(x);
        
        plot_data_time(z,'z',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);

        % Add cyclic prefix to create OFDM package
        zcp = add_cyclic_prefix(z,N_cp); %TODO: This line is missing some code!
        
        plot_data_time(zcp,'zcp',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);
        
        % Send package over channel
        ycp = simulate_baseband_channel(zcp, h, snr, sync_err);
        
        plot_data_time(ycp,'ycp',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);
        
        % Only keep the first N+Ncp received samples. Consider why ycp is longer
        % than zcp, and why we only need to save the first N+Ncp samples. This is
        % important to understand.
        % ANSWER: Since we are treating the signal as periodic, we can
        %         remove the initial and final Ncp values.
        ycp = ycp(1:N+N_cp);            

        % Remove cyclic prefix
        y = remove_cyclic_prefix(ycp,N_cp);

        plot_data_time(y,'y',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);

        % Convert to frequency domain using DFT
        r = fft(y);

        plot_FT(r,'r',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);
        
        symbs.rx_pe = r; % Store symbols for later

        % Remove effect of channel by equalization. Here, we can do this by
        % dividing r (which is in the frequency domain) by the channel gain (also
        % in the frequency domain).
        H = fft(h, N);
        r_eq = conj(H).*r./abs(H).^2; %TODO: This line is missing some code!
        
        plot_FT(r_eq,'hat s',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);
        
        symbs.rx_e = r_eq; %Store symbols for later

        % Calculate the quality of the received symbols.
        % The error vector magnitude (EVM) is one useful metric.
        evm = norm(x - r_eq)/sqrt(N);

        % Convert the recieved symsbols to bits
        rx = qpsk2bits(r_eq); %TODO: This line is missing some code!

        % Calculate the bit error rate (BER).
        % This indicates the relative number of bit errors.
        % Typically this will vary from 0 (no bit errors) to 0.5 (half of all
        % receieved bits are different, which is the number we'd expect if we
        % compare two random bit sequences).
        ber = 1-sum(rx == tx)/length(rx);
              
        plot_FT(r_eq./x,'H',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);
    end

    function txFrame = concat_packages(txPilot,txData)
        % Concatenate two ofdm blocks of equal size into a frame
        txPilot = txPilot(:);
        txData = txData(:);
        if(length(txData) ~= length(txPilot))
            error('Pilot and data are not of the same length!');
        end
        txFrame = [txPilot; txData]; %TODO: This line is missing some code!
    end

    function [rxPilot, rxData] = split_frame(rxFrame)
        % Split an ofdm frame into 2 equal ofdm packages
        rxFrame = rxFrame(:);
        if rem(length(rxFrame),2) > 0
            error('Vector z must have an even number of elements'); 
        end
        N = length(rxFrame);
        rxPilot = rxFrame(1:N/2); %TODO: This line is missing some code!
        rxData = rxFrame(N/2+1:end); %TODO: This line is missing some code!
    end

    function [rx, evm, ber, symbs] = sim_ofdm_unknown_channel(tx, h, N_cp, snr, sync_err)
        % Simulate OFDM signal transmission/reception over an unknown
        % channel.
        %
        % -----------------------------------------------------------------
        % NOTE: THIS FUNCTION WILL NOT BE SELF-TESTED!
        % It will be up to you to study the output from this function and
        % determine if the results are correct or not.
        % -----------------------------------------------------------------
        %
        % Arguments:
        %   tx          Structure with fields:
        %     .p        Pilot bits to transmit
        %     .d        Data bits to transmit
        %   h           Channel impulse response [-]
        %   N_cp        Cyclic prefix length [samples]
        %   snr         Channel signal/noise ration to apply [dB]
        %   sync_err    Reciever synchronization error [samples]
        % Outputs:
        %   rx          Recieved bits [-]
        %   evm         Error vector magnitude (see below) [-]
        %   ber         Bit error rate (see below) [-]
        %   symbs       Structure containing fields:
        %       .tx         Transmitted symbols
        %       .rx_pe      Recieved symbols, pre-equalization
        %       .rx_e       Recieved symbols, post-equalization
        %
        %
        % This function is similar to the known-channel problem, but with
        % the added complexity of requiring to estimate the channel
        % response. The relevant steps to perform here are:
        %   - Get a sequence of pilot and data bits to transmit
        %   - Convert the pilot and data bits to OFDM symbols
        %   - Create an OFDM block from the OFDM symbols for the pilot and
        %   data
        %   - Add a cyclic prefix to the pilot and data
        %   - Concatenate the pilot and data blocks to create an entire
        %   OFDM frame
        %   - Simulate the transmission and reception over the channel using the
        %   simulate_baseband_channe function.
        %   - Split the recieved message into a recieved pilot and data
        %   segment
        %   - Remove the cyclic prefixes from the recieved messages
        %   - Estimate the channel gain from the pilot block
        %   - Equalize the recieved data symbols by the channel gain
        %   - Convert the equalized symbols back to bits
        %   - Compare the recieved bits/symbols to the transmitted bits/symbols.

		warning('Note that this function is _not_ self-tested. It is up to you to study the output any verify that it is correct! You can remove this warning if you wish.');
		
        % Ensure inputs are column vectors
        tx.d = tx.d(:);
        tx.p = tx.p(:);
        h = h(:);
        
        % Convert bits to QPSK symbols
        x.p = bits2qpsk(tx.p); %TODO: This line is missing some code!
        x.d = bits2qpsk(tx.d); %TODO: This line is missing some code!

        symbs.tx = x.d;   % Store transmitted data symbols for later

        % Number of symbols in message
        N = length(x.d);
        if length(x.d) ~= length(x.p)
           error('Pilot and data messages must be of equal length'); 
        end

        % Create OFDM time-domain block using IDFT
        z.p = ifft(x.p); %TODO: This line is missing some code!
        z.d = ifft(x.d); %TODO: This line is missing some code!

        % Add cyclic prefix to create OFDM package
        zcp.p = add_cyclic_prefix(z.p,N_cp); %TODO: This line is missing some code!
        zcp.d = add_cyclic_prefix(z.d,N_cp); %TODO: This line is missing some code!
        
        % Concatenate the messages
        tx_frame = concat_packages(zcp.p,zcp.d); %TODO: This line is missing some code!
                
        plot_data_time(tx_frame,'z',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);

        % Send package over channel
        rx_frame = simulate_baseband_channel(tx_frame, h, snr, sync_err);
        
        plot_data_time(rx_frame,'y',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);
        
        % As before, only keep the first samples
        rx_frame = rx_frame(1:2*(N+N_cp));
        
        % Split frame into packages
        ycp = struct();
        [ycp.p, ycp.d] = split_frame(rx_frame); %TODO: This line is missing some code!
        
        % Remove cyclic prefix
        y.p = remove_cyclic_prefix(ycp.p,N_cp); %TODO: This line is missing some code!
        y.d = remove_cyclic_prefix(ycp.d,N_cp); %TODO: This line is missing some code!

        % Convert to frequency domain using DFT
        r.p = fft(y.p); %TODO: This line is missing some code!
        r.d = fft(y.d); %TODO: This line is missing some code!
        symbs.rx_pe = r.d; % Store symbols for later
        
        % Esimate channel
        H = r.p./x.p; %TODO: This line is missing some code!
        
        plot_FT(H,'H',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);
        
        plot_FT(r.d,'r_d',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);

        % Remove effect of channel on the data package by equalization.
        r_eq = conj(H).*r.d./abs(H).^2; %TODO: This line is missing some code!
        
        plot_FT(r_eq,'hat s',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);

        symbs.rx_e = r_eq; %Store symbols for later

        % Calculate the quality of the received symbols.
        % The error vector magnitude (EVM) is one useful metric.
        evm = norm(x.d - r_eq)/sqrt(N);

        % Convert the recieved symsbols to bits
        rx = qpsk2bits(r_eq); %TODO: This line is missing some code!

        % Calculate the bit error rate (BER).
        % This indicates the relative number of bit errors.
        % Typically this will vary from 0 (no bit errors) to 0.5 (half of all
        % receieved bits are different, which is the number we'd expect if we
        % compare two random bit sequences).
        ber = 1-sum(rx == tx.d)/length(rx);
    end

    function z = frame_interpolate(x,L,hlp)
        % Interpolate (upsample) a signal x by factor L, with an optionally
        % configurable lowpass filter.
        % Arguments:
        %   x   Signal to interpolate, length N
        %   L   Upsampling factor
        %   hlp FIR filter coefficents for lowpass filter, length Nh
        %       If not supplied, a default filter will be used with length
        %       62.
        % Returns:
        %   z   Interpolated signal of length N*L + Nh-1
        %
        
        if nargin < 3       % Default filter design
            SBscale = 1.7;  % Factor for stop band position
            Nfir = 61;      % The filter length if Nfir + 1
            hlp = firpm(Nfir, [0 1/L 1/L*SBscale 1], [1 1 0 0]);
        end
        
        % Make x, hlp column vectors
        x = x(:);
        hlp = hlp(:);
        
        figure('Color','white','Position',[376.2000e+000   408.2000e+000   644.8000e+000   218.4000e+000]);
        freqz(hlp);
        
        % Get the length of the input signal
        N = length(x);
        
        % Preallocate vector for upsampled, unfiltered, signal
        zup = zeros((N)*L,1);
        
        % Upsample by a factor L, i.e. insert L-1 zeros after each original
        % sample
        zup(1:L:end) = x; %TODO: This line is missing some code!
        
        % Apply the LP filter to the upsampled (unfiltered) signal.
        z = conv(zup,hlp); %TODO: This line is missing some code!
%         z = zup;
    end

    function z = frame_decimate(x,L,hlp)
        % Decimate (downsample) a signal x by factor L, with an optionally
        % configurable lowpass filter.
        % Arguments:
        %   x   Signal to decimate, length N
        %   L   Downsampling factor
        %   hlp FIR filter coefficents for lowpass filter, length Nh
        %       If not supplied, a default filter will be used with length
        %       61.
        % Returns:
        %   z   Interpolated signal of length N*L + Nh-1
        
        if nargin < 3       % Default filter design
            SBscale = 1.7;  % Factor for stop band position
            Nfir = 61;      % The filter length if Nfir + 1
            hlp = firpm(Nfir, [0 1/L 1/L*SBscale 1], [1 1 0 0]);
        end
        
        % Make x, hlp column vectors
        x = x(:);
        hlp = hlp(:);
        
        figure('Color','white','Position',[376.2000e+000   408.2000e+000   644.8000e+000   218.4000e+000]);
        freqz(hlp);
        
        % New Code
        N = length(x); % This N = Norignal-signal / M (here L is use)
        % Apply the lowpass filter to avoid aliasing when decimating
        xf = conv(x,hlp);   %TODO: This line is missing some code!
        
        % Downsample by keeping samples [1, 1+L, 1+2*L, ...]
        z = xf(1:L:end); %TODO: This line is missing some code!
    end

    function z = frame_modulate(x, theta)
       % Modulates a signal of length N with a modulation frequency theta.
       % Arguments:
       %    x       Signal to modulate of length N
       %    theta   Normalized modulation frequency
       % Outputs:
       %    z       Modulated signal
       
       % Make x a column vector
       x = x(:);
       
       N = length(x);
       
       % Generate vector of sample indices
       n = (0:N-1);
       n = n(:);
       
       % Modulate x by multiplying the samples with the complex exponential
       % exp(i * 2 * pi * theta * n)
       z = x.*exp(1i * 2 * pi * theta * n); %TODO: This line is missing some code!
    end

    function [rx, evm, ber, symbs] = sim_ofdm_audio_channel(tx, N_cp, snr, sync_err, f_s, f_c, L)
        % Simulate modulated OFDM signal transmission/reception over an
        % audio channel. This fairly accurately simulates the physical
        % channel of audio between a loudspeaker and a microphone.
        %
        % -----------------------------------------------------------------
        % NOTE: THIS FUNCTION WILL NOT BE SELF-TESTED!
        % It will be up to you to study the output from this function and
        % determine if the results are correct or not.
        % -----------------------------------------------------------------
        %
        % Arguments:
        %   tx          Structure with fields:
        %     .p        Pilot bits to transmit
        %     .d        Data bits to transmit
        %   N_cp        Cyclic prefix length [samples]
        %   snr         Channel signal/noise ration to apply [dB]
        %   f_s         The up-sampled sampling frequency [Hz]
        %   f_c         The modulation carrier frequency [Hz]
        %   L           The upsampling/downsampling factor [-]
        % Outputs:
        %   rx          Recieved bits [-]
        %   evm         Error vector magnitude (see below) [-]
        %   ber         Bit error rate (see below) [-]
        %   symbs       Structure containing fields:
        %       .tx         Transmitted symbols
        %       .rx_pe      Recieved symbols, pre-equalization
        %       .rx_e       Recieved symbols, post-equalization
        %
        %
        % This function is similar to the unknown-channel problem, but with
        % the added complexity of requiring to interpolate and modulate the
        % signal before transmission, followed by demodulation and
        % decimation on reception. The relevant steps to perform here are:
        %   - Get a sequence of pilot and data bits to transmit
        %   - Convert the pilot and data bits to OFDM symbols
        %   - Create an OFDM block from the OFDM symbols for the pilot and
        %   data
        %   - Add a cyclic prefix to the pilot and data
        %   - Concatenate the pilot and data blocks to create an entire
        %   OFDM frame
        %   - Interpolate the signal to a higher sample-rate
        %   - Modulate the signal, thereby moving it from the base-band to
        %   being centered about the modulation frequency.
        %   - Simulate the transmission and reception over the channel using the
        %   simulate_baseband_channe function.
        %   - Demodulate the signal, moving the recieved signal back to the
        %   base-band
        %   - Decimate the signal, reducing the sample-rate back to the
        %   original rate.
        %   - Split the recieved message into a recieved pilot and data
        %   segment
        %   - Remove the cyclic prefixes from the recieved messages
        %   - Estimate the channel gain from the pilot block
        %   - Equalize the recieved data symbols by the channel gain
        %   - Convert the equalized symbols back to bits
        %   - Compare the recieved bits/symbols to the transmitted bits/symbols.

		warning('Note that this function is _not_ self-tested. It is up to you to study the output any verify that it is correct! You can remove this warning if you wish.');
		
        % Ensure input is a column vector
        tx.d = tx.d(:);
        tx.p = tx.p(:);
        
        % Convert bits to QPSK symbols
        x.p = bits2qpsk(tx.p); %TODO: This line is missing some code!
        x.d = bits2qpsk(tx.d); %TODO: This line is missing some code!

        symbs.tx = x.d;   % Store transmitted data symbols for later

        % Number of symbols in message
        N = length(x.d);
        if length(x.d) ~= length(x.p)
           error('Pilot and data messages must be of equal length'); 
        end

        % Create OFDM time-domain block using IDFT
        z.p = ifft(x.p); %TODO: This line is missing some code!
        z.d = ifft(x.d); %TODO: This line is missing some code!

        % Add cyclic prefix to create OFDM package
        zcp.p = add_cyclic_prefix(z.p,N_cp); %TODO: This line is missing some code!
        zcp.d = add_cyclic_prefix(z.d,N_cp); %TODO: This line is missing some code!
        
        % Concatenate the messages
        tx_frame = concat_packages(zcp.p,zcp.d); %TODO: This line is missing some code!

        plot_FT_fs(fft(tx_frame),'z', f_s/L, 'z', [105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);
        
        % Increase the sample rate by interpolation
        tx_frame_us = frame_interpolate(tx_frame,L); %TODO: This line is missing some code!
        
        plot_FT_fs(fft(tx_frame_us),'z_{i}',f_s, 'z_i',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);
        
        % Modulate the upsampled signal
        tx_frame_mod = frame_modulate(tx_frame_us,f_c/f_s); %TODO: This line is missing some code!

        plot_FT_fs(fft(tx_frame_mod),'z_{m}',f_s, 'z_m',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);
        
        % Discard the imaginary part of the signal for transmission over a
        % scalar channel (simulation of audio over air)
        tx_frame_final = real(tx_frame_mod);
        
        plot_FT_fs(fft(tx_frame_final),'RE(z_{m})',f_s, 're_z_m',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);

        % Send package over channel
        [rx_frame_raw, rx_idx] = simulate_audio_channel(tx_frame_final, f_s, snr, sync_err);

        % Discard data before/after package
        rx_frame_raw = rx_frame_raw(rx_idx:rx_idx + length(tx_frame_final));
        
        plot_FT_fs(fft(rx_frame_raw),'y_{rec}',f_s, 'y_rec',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);
        
        % Demodulate to bring the signal back to the baseband
        rx_frame_us = frame_modulate(rx_frame_raw, -f_c/f_s); %TODO: This line is missing some code!
        
        plot_FT_fs(fft(rx_frame_us),'y_{ib}',f_s, 'y_ib',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);
        
        % Decimate the signal to bring the sample rate back to the original
        rx_frame = frame_decimate(rx_frame_us, L);
        
        plot_FT_fs(fft(rx_frame),'y',f_s/L, 'y',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);
        
        % Discard samples beyond OFDM frame
        rx_frame = rx_frame(1:2*(N+N_cp));

        % Split frame into packages
        ycp = struct();
        [ycp.p, ycp.d] = split_frame(rx_frame); %TODO: This line is missing some code!
        
        % Remove cyclic prefix
        y.p = remove_cyclic_prefix(ycp.p,N_cp); %TODO: This line is missing some code!
        y.d = remove_cyclic_prefix(ycp.d,N_cp); %TODO: This line is missing some code!

        % Convert to frequency domain using DFT
        r.p = fft(y.p); %TODO: This line is missing some code!
        r.d = fft(y.d); %TODO: This line is missing some code!
        symbs.rx_pe = r.d; % Store symbols for later
        
        % Esimate channel
        H = r.p./x.p; %TODO: This line is missing some code!
        
        plot_FT_fs(H,'H_{eff}',f_s/L, 'H_eff',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);
        
        plot_FT_fs(r.d,'r',f_s/L,'r',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);

        % Remove effect of channel on the data package by equalization.
        r_eq = conj(H).*r.d./abs(H).^2; %TODO: This line is missing some code!
        
        plot_FT_fs(r_eq,'hat s',f_s/L,'hat_s',[105.8000e+000   173.8000e+000     1.0888e+003   312.8000e+000]);

        symbs.rx_e = r_eq; %Store symbols for later

        % Calculate the quality of the received symbols.
        % The error vector magnitude (EVM) is one useful metric.
        evm = norm(x.d - r_eq)/sqrt(N);

        % Convert the recieved symsbols to bits
        rx = qpsk2bits(r_eq); %TODO: This line is missing some code!

        % Calculate the bit error rate (BER).
        % This indicates the relative number of bit errors.
        % Typically this will vary from 0 (no bit errors) to 0.5 (half of all
        % receieved bits are different, which is the number we'd expect if we
        % compare two random bit sequences).
        ber = 1-sum(rx == tx.d)/length(rx); 
    end



% Generate structure with handles to functions
funs.add_cyclic_prefix = @add_cyclic_prefix;
funs.remove_cyclic_prefix = @remove_cyclic_prefix;
funs.bits2qpsk = @bits2qpsk;
funs.qpsk2bits = @qpsk2bits;
funs.sim_ofdm_known_channel = @sim_ofdm_known_channel;
funs.concat_packages = @concat_packages;
funs.split_frame = @split_frame;
funs.sim_ofdm_unknown_channel = @sim_ofdm_unknown_channel;

funs.frame_interpolate = @frame_interpolate;
funs.frame_decimate = @frame_decimate;
funs.frame_modulate = @frame_modulate;
funs.sim_ofdm_audio_channel = @sim_ofdm_audio_channel;


% This file will return a structure with handles to the functions you have
% implemented. You can call them if you wish, for example:
% funs = student_sols();
% some_output = funs.some_function(some_input);

end
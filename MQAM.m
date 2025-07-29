function [EbN0_dB, SER, theo_SER] = simulate_MQAM(M, d, show_plot)
% Simulate M-QAM transmission over AWGN with RRC pulse shaping
% Inputs:
%   M         - Modulation order (e.g., 16, 64)
%   d         - Minimum distance between adjacent constellation points
%   show_plot - Boolean flag to plot SER results (true/false)
%
% Outputs:
%   EbN0_dB   - Eb/N0 values in dB
%   SER       - Simulated Symbol Error Rate
%   theo_SER  - Theoretical Symbol Error Rate

    L = 1e6;                 % Number of transmitted symbols
    f_ovsamp = 8;            % Oversampling factor
    span = 4;                % Span of RRC filter
    rolloff = 0.5;           % Rolloff factor
    EbN0_dB = 2:2:18;        % SNR values in dB
    EbN0_lin = 10.^(EbN0_dB / 10); % Linear SNR

    levels = (-sqrt(M)+1)*d/2 : d : (sqrt(M)-1)*d/2;
    Eavg = mean(levels.^2) * d^2 / 2; % Average symbol energy
    Eb = Eavg / log2(M);              % Energy per bit

    % Generate random symbols
    I = levels(randi(length(levels), L, 1));
    Q = levels(randi(length(levels), L, 1));
    sk = I + 1j * Q;

    % Upsample
    s_up = upsample(sk, f_ovsamp);

    % RRC pulse shaping
    pc = rcosdesign(rolloff, span, f_ovsamp, "sqrt");
    pc = pc / norm(pc);         % Normalize filter energy
    pcmatch = pc;               % Matched filter = same as shaping (time-reversal not needed with symmetric filters)
    xrcos = conv(s_up, pc);     % Transmit signal

    noise = randn(1, length(xrcos)) + 1j * randn(1, length(xrcos));
    SER = zeros(size(EbN0_dB));

    for i = 1:length(EbN0_dB)
        SNR = EbN0_lin(i);
        sigma = sqrt(Eb / (2 * SNR));    % Noise std dev
        awgn = sigma * noise;

        received = xrcos + awgn;

        delay = span * f_ovsamp;  % Delay = span * oversampling
        z = conv(received, pcmatch);
        z = z(delay+1 : f_ovsamp : delay + L*f_ovsamp);  % Sample matched output

        % Hard decision detector
        I_hat = zeros(1, length(z));
        Q_hat = zeros(1, length(z));
        r_hat = zeros(1, length(z));
        for k = 1:length(z)
            [~, idxI] = min(abs(real(z(k)) - levels));
            [~, idxQ] = min(abs(imag(z(k)) - levels));
            I_hat(k) = levels(idxI);
            Q_hat(k) = levels(idxQ);
            r_hat(k) = I_hat(k) + 1j * Q_hat(k);
        end
        SER(i) = mean(sk ~= r_hat);
    end

    % Theoretical SER
    arg_erfc = sqrt((3 * log2(M) / (2 * (M - 1))) * EbN0_lin);
    theo_SER = 1 - (1 - (1 - 1/sqrt(M)) * erfc(arg_erfc)).^2;

    % Plot results
    if show_plot
        figure;
        semilogy(EbN0_dB, SER, 'bo-'); grid on; hold on;
        semilogy(EbN0_dB, theo_SER, 'r--');
        xlabel('E_b/N_0 (dB)');
        ylabel('Symbol Error Rate');
        title(sprintf('%d-QAM SER over AWGN with RRC pulse shaping', M));
        legend('Simulated SER', 'Theoretical SER');
    end
end

[EbN0,SER1,SER2]=simulate_MQAM(64,2,true);
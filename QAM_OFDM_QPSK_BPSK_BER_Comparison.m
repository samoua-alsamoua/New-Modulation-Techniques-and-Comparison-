% Project: New Modulation Techniques and Comparison
% Objective: Implement and compare advanced modulation techniques (QAM, OFDM) with traditional ones (BPSK, QPSK)

clc;
clear;
close all;

%% Parameters
M_qam = [16, 64]; % Modulation orders for QAM (16-QAM, 64-QAM)
snr_range = 0:2:20; % SNR range in dB
num_bits = 1e5; % Number of bits to simulate
EbNo = snr_range + 10*log10(log2(M_qam(1))); % Convert SNR to Eb/No for QAM
ber_bpsk = zeros(size(snr_range)); % BER for BPSK
ber_qpsk = zeros(size(snr_range)); % BER for QPSK
ber_16qam = zeros(size(snr_range)); % BER for 16-QAM
ber_64qam = zeros(size(snr_range)); % BER for 64-QAM
ber_ofdm = zeros(size(snr_range)); % BER for OFDM

%% OFDM Parameters
num_subcarriers = 64; % Number of subcarriers
cyclic_prefix_length = 16; % Length of cyclic prefix
fft_size = num_subcarriers; % FFT size
mod_order_ofdm = 4; % Modulation order for OFDM (QPSK)

%% Simulation Loop
for i = 1:length(snr_range)
    snr = snr_range(i);
    
    %% BPSK Modulation
    data_bpsk = randi([0 1], num_bits, 1); % Generate random binary data
    mod_bpsk = pskmod(data_bpsk, 2); % BPSK modulation
    rx_bpsk = awgn(mod_bpsk, snr, 'measured'); % Add noise
    demod_bpsk = pskdemod(rx_bpsk, 2); % BPSK demodulation
    [~, ber_bpsk(i)] = biterr(data_bpsk, demod_bpsk); % Calculate BER
    
    %% QPSK Modulation
    data_qpsk = randi([0 1], num_bits, 1); % Generate random binary data
    mod_qpsk = pskmod(data_qpsk, 4, pi/4); % QPSK modulation
    rx_qpsk = awgn(mod_qpsk, snr, 'measured'); % Add noise
    demod_qpsk = pskdemod(rx_qpsk, 4, pi/4); % QPSK demodulation
    [~, ber_qpsk(i)] = biterr(data_qpsk, demod_qpsk); % Calculate BER
    
    %% 16-QAM Modulation
    data_16qam = randi([0 1], num_bits, 1); % Generate random binary data
    mod_16qam = qammod(data_16qam, M_qam(1)); % 16-QAM modulation
    rx_16qam = awgn(mod_16qam, snr, 'measured'); % Add noise
    demod_16qam = qamdemod(rx_16qam, M_qam(1)); % 16-QAM demodulation
    [~, ber_16qam(i)] = biterr(data_16qam, demod_16qam); % Calculate BER
    
    %% 64-QAM Modulation
    data_64qam = randi([0 1], num_bits, 1); % Generate random binary data
    mod_64qam = qammod(data_64qam, M_qam(2)); % 64-QAM modulation
    rx_64qam = awgn(mod_64qam, snr, 'measured'); % Add noise
    demod_64qam = qamdemod(rx_64qam, M_qam(2)); % 64-QAM demodulation
    [~, ber_64qam(i)] = biterr(data_64qam, demod_64qam); % Calculate BER
    
    %% OFDM Modulation
    data_ofdm = randi([0 1], num_bits, 1); % Generate random binary data

    % Ensure the total number of bits is divisible by the number of subcarriers
    num_symbols_per_carrier = floor(num_bits / (log2(mod_order_ofdm) * fft_size)); % Number of OFDM symbols
    total_bits_needed = num_symbols_per_carrier * log2(mod_order_ofdm) * fft_size; % Total bits for OFDM
    data_ofdm = data_ofdm(1:total_bits_needed); % Truncate data to fit OFDM structure

    % Modulate data using QPSK (or any other modulation scheme)
    data_symbols = qammod(data_ofdm, mod_order_ofdm); % QPSK modulation for OFDM

    % Reshape into OFDM symbols (each row corresponds to one OFDM symbol)
    data_symbols = reshape(data_symbols, fft_size, []).'; % Reshape into OFDM symbols
    
    % Perform IFFT to generate OFDM symbols
    ofdm_symbols = ifft(data_symbols, fft_size, 2);
    
    % Add cyclic prefix
    cp = ofdm_symbols(:, end-cyclic_prefix_length+1:end); % Extract cyclic prefix
    ofdm_symbols_with_cp = [cp, ofdm_symbols]; % Concatenate cyclic prefix
    
    % Flatten the OFDM symbols for transmission
    tx_ofdm = reshape(ofdm_symbols_with_cp.', [], 1);
    
    % Add noise
    rx_ofdm = awgn(tx_ofdm, snr, 'measured');
    
    % Remove cyclic prefix and perform FFT
    rx_ofdm = reshape(rx_ofdm, [], size(ofdm_symbols_with_cp, 2)); % Reshape back
    rx_ofdm_no_cp = rx_ofdm(:, cyclic_prefix_length+1:end); % Remove cyclic prefix
    rx_symbols = fft(rx_ofdm_no_cp, fft_size, 2); % Perform FFT
    
    % Demodulate and calculate BER
    rx_data = qamdemod(rx_symbols, mod_order_ofdm);
    rx_data = reshape(rx_data.', [], 1);
    [~, ber_ofdm(i)] = biterr(data_ofdm, rx_data); % Calculate BER
end

%% Plot BER vs SNR
figure;
semilogy(snr_range, ber_bpsk, 'b-o', 'LineWidth', 1.5); hold on;
semilogy(snr_range, ber_qpsk, 'r-s', 'LineWidth', 1.5);
semilogy(snr_range, ber_16qam, 'g-d', 'LineWidth', 1.5);
semilogy(snr_range, ber_64qam, 'm-^', 'LineWidth', 1.5);
semilogy(snr_range, ber_ofdm, 'k-p', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for Different Modulation Techniques');
legend('BPSK', 'QPSK', '16-QAM', '64-QAM', 'OFDM', 'Location', 'Best');
hold off;

%% Spectral Efficiency Comparison
spectral_efficiency = [1, 2, log2(M_qam(1)), log2(M_qam(2)), log2(mod_order_ofdm)];
modulation_types = {'BPSK', 'QPSK', '16-QAM', '64-QAM', 'OFDM'};

figure;
bar(spectral_efficiency, 'FaceColor', [0.2 0.4 0.6]);
set(gca, 'XTickLabel', modulation_types);
xlabel('Modulation Technique');
ylabel('Spectral Efficiency (bits/symbol)');
title('Spectral Efficiency Comparison');
grid on;

%% Computational Complexity Analysis
fprintf('Computational Complexity:\n');
fprintf('BPSK: Low\n');
fprintf('QPSK: Low\n');
fprintf('16-QAM: Moderate\n');
fprintf('64-QAM: High\n');
fprintf('OFDM: High (due to FFT/IFFT operations)\n');
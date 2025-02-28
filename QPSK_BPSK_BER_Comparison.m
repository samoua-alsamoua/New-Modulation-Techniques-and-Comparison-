% Project: Comparison Between BPSK and QPSK Modulation Techniques
% Objective: Compare BER performance of BPSK and QPSK under AWGN channel

clc;
clear;
close all;

%% Parameters
snr_range = 0:2:20; % SNR range in dB
num_bits = 1e5; % Number of bits to simulate
ber_bpsk = zeros(size(snr_range)); % BER for BPSK
ber_qpsk = zeros(size(snr_range)); % BER for QPSK

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
    mod_qpsk = pskmod(data_qpsk, 4, pi/4); % QPSK modulation with phase offset
    rx_qpsk = awgn(mod_qpsk, snr, 'measured'); % Add noise
    demod_qpsk = pskdemod(rx_qpsk, 4, pi/4); % QPSK demodulation
    [~, ber_qpsk(i)] = biterr(data_qpsk, demod_qpsk); % Calculate BER
end

%% Plot BER vs SNR
figure;
semilogy(snr_range, ber_bpsk, 'b-o', 'LineWidth', 1.5); hold on;
semilogy(snr_range, ber_qpsk, 'r-s', 'LineWidth', 1.5);
grid on;
xlabel('SNR (dB)');
ylabel('Bit Error Rate (BER)');
title('BER vs SNR for BPSK and QPSK');
legend('BPSK', 'QPSK', 'Location', 'Best');
hold off;

%% Spectral Efficiency Comparison
fprintf('Spectral Efficiency:\n');
fprintf('BPSK: %.2f bits/symbol\n', 1); % BPSK transmits 1 bit per symbol
fprintf('QPSK: %.2f bits/symbol\n', 2); % QPSK transmits 2 bits per symbol
%% Slippage Detection

clc;
clear;

tactile_raw_data_1_preprocess = dlmread('tactile_fine2_1.dat');

tactile_raw_data_1 = tactile_raw_data_1_preprocess(100:300,:); % choose plot data range

[M_row(1),N_col(1)] = size(tactile_raw_data_1);

Force_resultant_1 = zeros(1,M_row(1));
for I = 1:M_row(1)
    cor_1(I) = corr(tactile_raw_data_1(1,:)',tactile_raw_data_1(I,:)');
    if(isnan(cor_1(I)))
        cor_1(I) = 0;
    end
    for J = 1:N_col(1)
            if(tactile_raw_data_1(I,J) ~= 0)
                Force_resultant_1(I) = (tactile_raw_data_1(I,J)*1.6575 + 77.0588) * 61 * (0.34)^2 * 10^(-4) + Force_resultant_1(I);
            end
    end
end

period = 4; % Calculate FFT every period samples
NFFT = 2*nextpow2(period);
for flag_count = 1:floor(M_row(1)/period)
    cor_1_fft_process = abs(fft(cor_1(1+(flag_count-1)*period:period+(flag_count-1)*period),NFFT));
    cor_frequency_second_comp_1(flag_count) = cor_1_fft_process(2); % first is the DC component, second is the AC component
end

sampling_rate = 20;
font_size = 16;
size_fft = size(cor_frequency_second_comp_1);
time = 1/sampling_rate:1/sampling_rate:M_row(1)/sampling_rate;
time_fft = 1/5:1/5:size_fft(2)/5;
figure(1)
subplot(311)
for I = 1:size_fft(2)
    if(cor_frequency_second_comp_1(I) > 0.005)
        stem(time_fft(I),cor_frequency_second_comp_1(I),'Color','red');hold on;
    else
        stem(time_fft(I),cor_frequency_second_comp_1(I),'Color','blue');hold on;
    end
end
grid on;
% stem(time_fft,cor_frequency_second_comp_1);grid on;
set(gca,'LineWidth',2)
set(gca,'FontSize',font_size);
title('FFT-First Frequency Component');
xlabel('Time/s');
ylabel('Amplitude');


subplot(312)
plot(time,cor_1,'LineWidth',2); grid on;
set(gca,'LineWidth',2)
set(gca,'FontSize',font_size);
title('Correlation Coefficient');
xlabel('Time/s');
ylim([-0.5,1.2]);


subplot(313)
plot(time,Force_resultant_1,'LineWidth',2); grid on;
set(gca,'LineWidth',2)
set(gca,'FontSize',font_size);
title('Grasp Force');
ylabel('Force/N');
xlabel('Time/s');
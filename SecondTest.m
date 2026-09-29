clear; close all; clc;

DAS_points = 500;

distance = 1.65;
c = 300000000;

data = readmatrix('ranging_data.csv'); 
x = data(:,1); 
y = data(:,2:4);

time_diff = x(2) - x(1);
d = x.*c;

in_sig = data(:,2);
a1_sig = data(:,3);
a2_sig = data(:,4);


figure
plot(d,y, 'LineWidth', 1.2) 
xlabel('Distance')
ylabel('Signal')
legend('In signal','Ant 1','Ant 2','Location','southeast')
grid on


Fs = 1 / time_diff;   
N = length(in_sig);
nfft = 2^nextpow2(N);
f = (0:(nfft/2))' * (Fs / nfft);

% FFTs
IN = fft(in_sig, nfft);
A1 = fft(a1_sig, nfft);
A2 = fft(a2_sig, nfft);


mag_in = abs(IN(1:nfft/2+1)) / N;
mag_a1 = abs(A1(1:nfft/2+1)) / N;
mag_a2 = abs(A2(1:nfft/2+1)) / N;

phase_in = unwrap(angle(IN(1:nfft/2+1)));
phase_a1 = unwrap(angle(A1(1:nfft/2+1)));
phase_a2 = unwrap(angle(A2(1:nfft/2+1)));


% plot magnitude
figure
plot(f(1:200), mag_in(1:200), 'LineWidth',1.2); hold on
plot(f(1:200), mag_a1(1:200), 'LineWidth',1.2)
plot(f(1:200), mag_a2(1:200), 'LineWidth',1.2)
xlabel('Frequency (Hz)')
ylabel('Magnitude')
legend('in\_sig','a1\_sig','a2\_sig','Location','best')
grid on

% plot phase
figure
plot(f(1:500), phase_in(1:500), 'LineWidth',1.2); hold on
plot(f(1:500), phase_a1(1:500), 'LineWidth',1.2)
plot(f(1:500), phase_a2(1:500), 'LineWidth',1.2)
xlabel('Frequency (Hz)')
ylabel('Phase  (rads)')
legend('in\_sig','a1\_sig','a2\_sig','Location','best')
grid on

samples = length(a1_sig);


b = samples / DAS_points;
shift_time = time_diff * b;
summed_energy = zeros(DAS_points,1);


for i = 1:DAS_points
    
    k = round(i*b)/2;
    %shift_sig1 = zeros(size(a1_sig));
    %shift_sig2 = zeros(size(a2_sig));

    shift_sig1 = circshift(a1_sig, -k);
    shift_sig2 = circshift(a2_sig, -(DAS_points-k));

    

    summed1 = shift_sig1 + shift_sig2;
    squared1 = summed1 .^ 2;

    
    e1 = sum(squared1);
    
    
    
    summed_energy(i) = e1;
end

step = (shift_time * c) / 2;




sig_dis = (0:(DAS_points-1))' * (shift_time * c) / 2;



figure
plot(sig_dis,summed_energy, 'LineWidth', 1.2)  

title('Summed Energy')

xlabel('Distance') 
ylabel('Energy');
grid on

ave_energy = mean(summed_energy);
abs_energy = abs(summed_energy - ave_energy);

figure
plot(sig_dis,abs_energy, 'LineWidth', 1.2)  

title('Distance from Mean Energy')
xlabel('Distance') 
ylabel('Energy');
grid on


window_len = 20;
sm_abs_energy = movmean(abs_energy, window_len);


figure
plot(sig_dis, sm_abs_energy, 'LineWidth',1.2)
title('Smoothed Distance from Mean Energy')
xlabel('Distance')
ylabel('Energy')
grid on


heatmat = repmat(summed_energy(:),1,DAS_points);

figure
imagesc(sig_dis, 1:1, heatmat')
set(gca,'YDir','normal')
xlabel('Distance') 

title('Summed Energy Heatmap')
colorbar
colormap(jet)


heatmat = repmat(abs_energy(:),1,DAS_points);

figure
imagesc(sig_dis, 1:1, heatmat')
set(gca,'YDir','normal')
xlabel('Distance') 

title('Distance from Mean Energy Heatmap')
colorbar
colormap(jet)


heatmat = repmat(sm_abs_energy(:),1,DAS_points);

figure
imagesc(sig_dis, 1:1, heatmat')
set(gca,'YDir','normal')
xlabel('Distance') 

title('Distance from Mean Energy Heatmap')
colorbar
colormap(jet)
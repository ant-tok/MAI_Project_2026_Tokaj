clear; close all; clc;
% Add more layers, maybe model arm, look at TOF, make git

fs_permittivity = 8.854187817e-12; % permittivity of free space 
fs_permeability = pi*4e-7; % permeability of free space 


rel_permittivity_1 = 1; 
rel_permittivity_2 = 4; % relative permittivity of reflector

reflector_dist = 3.0;
zmax = 5.0;


frequency_spacing = 1e7;
max_freq = 3e9; 
no_tones = max_freq/frequency_spacing;
freq_list = (1:no_tones)'*frequency_spacing;

% Values for pulse
f0 = 1e9; %Carrier frequency
tp = 1e-9; % Envelope width



%c0 = 1/sqrt(fs_permittivity*fs_permeability); 
c0 = 3e8;


v1 = c0 * 1/sqrt(rel_permittivity_1); 
imped_1 = sqrt(fs_permeability/(fs_permittivity*rel_permittivity_1));

v2 = c0 * 1/sqrt(rel_permittivity_2); 
imped_2 = sqrt(fs_permeability/(fs_permittivity*rel_permittivity_2));

Reflection_coef = (imped_2 - imped_1)/(imped_2 + imped_1);
Transmission_coef = 2*imped_2/(imped_2 + imped_1);

fprintf('Gamma (Reflection Coefficient) = %+.4f,  Tau (Trasnmission Coefficient) = %.4f\n', Reflection_coef, Transmission_coef);


% Vectors of wave numbers
k1 = 2*pi*freq_list/v1;   
k2 = 2*pi*freq_list/v2;


wavelength = v1/f0;
dist_step = wavelength/40;
z = (0:dist_step:zmax); 
Nz = numel(z);

m1 = (z < reflector_dist);  
m2 = ~m1;

Ei = exp(-1j*k1*z).*m1;
Hi = -Ei/imped_1;

Er = Reflection_coef*exp(-1j*k1*(2*reflector_dist - z)).*m1;
Hr = Er/imped_1;

Et = Transmission_coef*exp(-1j*k1*reflector_dist).*exp(-1j*k2*(z - reflector_dist)).*m2;
E = Ei + Er + Et;

Ht = -Et/imped_2;
H  = Hi + Hr  + Ht;



Y = E(:, 1);


Frequency_Weights = (tp*sqrt(pi)/2) * exp(-1j*2*pi*freq_list*3*tp);

Frequency_Weights = Frequency_Weights .* exp(-(pi*(freq_list - f0)*tp).^2);   % Gaussian bell around f0

response = Frequency_Weights.*Y;

Ifft_Length = 2^12;

time_axis = (0:Ifft_Length-1)'/(Ifft_Length*frequency_spacing);  


synth = @(S) real(ifft([0; Ifft_Length*frequency_spacing*S; zeros(Ifft_Length-2*no_tones-1,1); conj(flipud(Ifft_Length*frequency_spacing*S))]));
full_synth = @(S) real(ifft([zeros(1,Nz); Ifft_Length*frequency_spacing*S; zeros(Ifft_Length-2*no_tones-1, Nz);conj(flipud(Ifft_Length*frequency_spacing*S))]));

recieved = synth(Y);
transmitted = synth(Frequency_Weights); 
response_signal = synth(response);

response_full = Frequency_Weights .* E; 
Etz = full_synth((response_full));


figure;

plot(freq_list/1e9, abs(Frequency_Weights), 'LineWidth', 1.2); grid on
xlabel('Frequency (GHz)'); ylabel('|weight|');
title('Pulse Spectrum');


figure;
subplot(3,1,1);

plot(time_axis*1e9, recieved); grid on
xlim([0 30]); xlabel('t (ns)'); ylabel('E (V/m)');
title('Recieved Signal in Time');

subplot(3,1,2);

plot(time_axis*1e9, transmitted); grid on
xlim([0 30]); xlabel('t (ns)'); ylabel('E (V/m)');
title('Pulse');

subplot(3,1,3);

plot(time_axis*1e9, response_signal); grid on
xlim([0 30]); xlabel('t (ns)'); ylabel('E at Rx (V/m)');
title('Recieved Conjugated with pulse');


figure;
imagesc(z, time_axis*1e9, Etz);
axis xy;
colormap(jet);
colorbar;
ylim([0 30]); xlabel('distance (m)'); ylabel('time (ns)');
title('E at distance/time');
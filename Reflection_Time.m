clear; close all; clc;

fs_permittivity = 8.854187817e-12; % permittivity of free space 
fs_permeability = pi*4e-7; % permeability of free space 


rel_permittivity_1 = 1; 
rel_permittivity_2 = 4; % relative permittivity of reflector

reflector_dist = 3.0;
reciever_dist  = 0;
zmax = 5.0; 


c0 = 1/sqrt(fs_permittivity*fs_permeability); 

v1 = c0 * 1/sqrt(rel_permittivity_1); 
imped_1 = sqrt(fs_permeability/(fs_permittivity*rel_permittivity_1));

v2 = c0 * 1/sqrt(rel_permittivity_2); 
imped_2 = sqrt(fs_permeability/(fs_permittivity*rel_permittivity_2));


Reflection_coef = (imped_2 - imped_1)/(imped_2 + imped_1);
Transmission_coef = 2*imped_2/(imped_2 + imped_1);

fprintf(['Reflection Coefficient = %+.4f,  Transmission Coefficient= %.4f  ' ...
    '    (1+Gamma = %.4f)\n'], Reflection_coef, Transmission_coef, 1+Reflection_coef);

f0 = 1e9; 
tp = 1e-9;     

%Pulse centred on zero
pulse = @(t) exp(-((t-3*tp)/tp).^2).*cos(2*pi*f0*(t-3*tp));


wavelength = v1/f0; 

dz = wavelength/50;     

z  = (0:dz:zmax); 

Sim_Time = 3e-8;  

time_step = Sim_Time/(length(z)-1);


t = (0:time_step:Sim_Time)'; 


[Z, T] = meshgrid(z, t);
m1 = (Z < reflector_dist); 
m2 = ~m1;  
Ei = pulse(T - Z/v1).*m1;
Hi = -Ei/imped_1;


Er = Reflection_coef*pulse(T -(2*reflector_dist-Z)/v1).*m1;
Hr = Er/imped_1;

E = Ei + Er;
H = Hi + Hr;


Et = Transmission_coef*pulse(T - reflector_dist/v1 - (Z - reflector_dist)/v2).*m2;     % Eq. (8.86)
E  = E + Et;

Ht = -Et/imped_2;
H  = H  + Ht;


figure;
imagesc(z, t*1e9, E); set(gca,'YDir','normal'); colormap(jet); colorbar

xlabel('distance (m)'); ylabel('time (ns)'); title('E at distance/time');


figure;
plot(t*1e9, E(:, 1), 'LineWidth', 1.2); grid on
xlabel('t (ns)'); ylabel('E_y at Rx (V/m)');
title('Signal at reciever, E');

figure;
plot(t*1e9, H(:, 1), 'LineWidth', 1.2); grid on
xlabel('t (ns)'); ylabel('H_y at Rx (A/m)');
title('Signal at reciever, H');
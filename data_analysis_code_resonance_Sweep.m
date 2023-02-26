% data analysis file for the resonance simulation data from waveguide
% solver 
clear; clc;
%% import and store the data from the txt file 
data = importdata("OptimizeListFull_2.txt");
Q_sim = data.data(:,1);
V_sim = data.data(:,2); 
resonance_f = data.data(:,3);
resonance_wavelength = (3e8)./resonance_f;
target_wavelength = 0.916*(10^(-6)); % target wavelength in m 
wavelength_detune = resonance_wavelength - target_wavelength;
%% plot the Q 
figure; hold off;
plot(Q_sim);
xlabel("Sim run"); ylabel("simulated Q");
%% plot the mode volume 
figure; hold off;
plot(V_sim);
xlabel("Sim run"); ylabel("simulated mode volume");
%% plot the wavelength detune 
figure; hold off;
plot(wavelength_detune);
xlabel("Sim run"); ylabel("simulated wavelength detune");
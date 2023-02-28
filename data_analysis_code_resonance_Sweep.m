% data analysis file for the resonance simulation data from waveguide
% solver 
clear; clc;
%% import and store the data from the txt file 
data = importdata("OptimizeListFull_waveguide_Sweep_char.csv");
cellNum = data(:,1);
Q_sim = data(:,2); 
resonance_f = data(:,4);
resonance_wavelength = (3e8)./resonance_f;
wavelength_detune = data(:,5);
calc_fitness = data(:,6)
%% plot the Q 
figure; hold off;
plot(cellNum,Q_sim);
xlabel("Weak mirror cell num"); ylabel("simulated Q");
%% plot the mode volume 
figure; hold off;
plot(cellNum,wavelength_detune);
xlabel("Weak mirror cell num"); ylabel("simulated detuning (m)");
%% plot the wavelength detune 
figure; hold off;
plot(cellNum,calc_fitness);
xlabel("Weak mirror cell num"); ylabel("simulated fitness");
% data analysis file for the resonance simulation data from waveguide
% solver 
clear; clc;
%% import and store the data from the txt file 
data = importdata("SiC_500nm_rib_unitcell_testSweep_t5.txt");
dielectric_band = data(:,6);
air_band = data(:,7);
mg = data(:,5);
bg = data(:,4);
bg_mg_rat = data(:,4);
%% plot the mg of all the trials
target_freq = 327.3e12;
wavelength_detuning = ((3e8)/target_freq-(3e8)./mg).*1e9;
%% plot the detuning of the mg from the target
figure; plot(wavelength_detuning);
%% plot the bandgap
figure; plot(bg);
%% plot the fitness 
wavelength_tolerance = 5;
wavelength_pen = exp(-((wavelength_detuning)./wavelength_tolerance).^2);
fitness = bg_mg_rat.*wavelength_pen;
figure; plot(fitness);



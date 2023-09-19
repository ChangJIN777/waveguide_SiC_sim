% data analysis file for the resonance simulation data from waveguide
% solver 
clear; clc;
%% import and store the data from the txt file 
data = importdata("SiC_500nm_rib_unitcell_testSweep_t3.txt");
mg = data(:,5);
%% plot the mg of all the trials
target_freq = 327.3e12;
wavelength_detuning = ((3e8)/target_freq-(3e8)./mg).*1e9;
plot(wavelength_detuning);
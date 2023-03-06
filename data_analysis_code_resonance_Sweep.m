% data analysis file for the resonance simulation data from waveguide
% solver 
clear; clc;
%% import and store the data from the txt file 
data = importdata("OptimizeListFull_resonance_sweep_v7.csv");
a = data(:,1);
sim_run = linspace(1,length(a),length(a));
d = data(:,2); 
w = data(:,3);
t = data(:,4);
Q_sim = data(:,5);
Q_sc = data(:,6);
Q_wvg = data(:,7);
V_mode = data(:,8);
wavelength_detune = data(:,9);
calc_fitness = data(:,10);
wavelength_detune = abs(wavelength_detune);
wavelength_detune = wavelength_detune * (10^9); % converting to nm
%% plot the calculated fitness
figure; hold off;
plot(sim_run,calc_fitness);
xlabel("simulation run"); ylabel("simulated fitness");
%% plot the wavelength detune
figure; hold off;
plot(sim_run,wavelength_detune);
xlabel("simulation run"); ylabel("simulated detuning (nm)");
%% plot the detuning vs the lattice constant
figure; hold off;
plot(a,wavelength_detune,'ro');
xlabel("lattice constant"); ylabel("wavelength detuning (nm)");
%% plot the Q_waveguide
figure; hold off;
plot(sim_run,Q_wvg);
xlabel("simulation run"); ylabel("simulated Q waveguide");
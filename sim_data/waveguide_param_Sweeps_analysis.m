% data analysis file for the resonance simulation data from waveguide
% solver 
clear; clc;
%% import and store the data from the txt file 
data = importdata("OptimizeListFull_with_waveguide_test_sweep_t3.csv");
weak_mirror_prefac = data(:,1);
sim_run = linspace(1,length(weak_mirror_prefac),length(weak_mirror_prefac));
waveguide_prefac = data(:,2); 
Q_sim = data(:,3);
Q_sc = data(:,4);
Q_wvg = data(:,5);
V_mode = data(:,6);
%% plot the Q 
figure; hold off;
plot(weak_mirror_prefac,Q_wvg,'ro');
xlabel("weak mirror prefactor"); ylabel("simulated Q waveguide");
%% plot the wavelength detune
figure; hold off;
plot(waveguide_prefac,Q_wvg,'ro');
xlabel("waveguide prefactor"); ylabel("simulated Q waveguide");
%% plot the Q_scatter
figure; hold off;
plot(sim_run,Q_wvg,'ro');
xlabel("simulation run"); ylabel("simulated Q waveguide");
%% plot the Q_waveguide
figure; hold off;
plot(sim_run,Q_wvg);
xlabel("simulation run"); ylabel("simulated Q waveguide");
%% importing the appropriate data 
clear; clc;
file_name = "cavity_500nm_optimization_trial2.csv";
data = importdata(file_name);
Vmode = data(:,8);
Qwvg = data(:,9);
Qsc = data(:,10);
Q = data(:,15);
detuning_wavelength = data(:,17);
%% account for unrealistically small mode volumes 
for i = 1:length(Vmode)
    if Vmode(i) < 0.4
        Vmode(i) = 1000000;
    end
end 
%% plotting the fitness 
fitness = (Q.*Qsc)./(Vmode.*Vmode).*exp(-(detuning_wavelength/(1e-9)).^2);
figure; plot(fitness);
%% plotting the Q ratios
ratio = Qsc./Qwvg;
figure; plot(ratio);
%% plotting Qsc vs Qwvg
figure; hold on; 
plot(Qsc); plot(Qwvg);


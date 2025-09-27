%% Figure2a_Simulation.m
% This script generates **Figure 2a** from the paper:
% "Functional stability and recurrent STDP in rhythmogenesis"
%
% Description:
% Computes and plots the phase diagram of the neural network in the
% [J_EI, J_IE] plane, showing the oscillation frequency as a contour plot.
% Uses parfor loop for computing the frequency contour
%
% Dependencies:
% - FrequencyContourJeiJie.m 
%
% Authors: Gabi Socolovsky & Maoz Shamir
% Date: 2025-07-09

%% Parameters
% Units
Tunits  = 5e-3;      % Time unit = Tau_m = 5 ms

% Time and delay
dt      = 0.01;      % Time step 
D       = 0.4;       % Synaptic delay
tf      = 200;       % Simulation Duration

% Order parameters of Synaptic weights
Jee     = 0.6; % Excitatory to excitatory synaptic weight (order parameter)
Jii     = 0.4; % Inhibitory to inhibitory synaptic weight (order parameter)
Jhat    = (Jee + Jii) / 2;

%% Find bifurcation frequency and J̄_D value
syms wD JbarD

if Jee >= Jii
    range = [0.1 5; 0.01 pi / (2*D)];
    Y = vpasolve([
        sqrt(JbarD^2 - Jee*Jii) == 1 / cos(wD*D - acos((Jee - Jii) / (2*sqrt(JbarD^2 - Jee*Jii)))),
        wD == -tan(wD*D - acos((Jee - Jii) / (2*sqrt(JbarD^2 - Jee*Jii))))
    ], [JbarD, wD], range);
else
    error('This script does not support Jee < Jii configuration.');
end

% Extract numerical values
wD    = double(Y.wD);       % Frequency on bifurcation line
JbarD = double(Y.JbarD);    % J-bar on the bifurcation line

%% Compute frequency contour in [Jei, Jie] plane
Jei_arr = 0:0.01:(1 + Jii);
Jie_arr = 0:0.02:20;
f       = nan(length(Jei_arr), length(Jie_arr));

parfor i = 1:length(Jei_arr)
    Jei = Jei_arr(i);
    f(i, :) = FrequencyContourJeiJie(JbarD, Jee, Jei, Jie_arr, Jii, dt, tf);
end

f = f / Tunits; % Convert frequency to Hz

%% Plot phase diagram
[JEI_arr, JIE_arr] = meshgrid(Jei_arr, Jie_arr);
[C, h] = contourf(JEI_arr, JIE_arr, f', 100);
set(h, 'LineColor', 'none');
grid on; hold on;

% Overlay bifurcation line: J_IE = J̄_D² / J_EI
Jei_arr = 0.01:0.01:(1 + Jii);
plot(Jei_arr, JbarD^2 ./ Jei_arr, 'k', 'LineWidth', 3.5);

% Vertical asymptote at J_EI = 1 + J_II
plot((1 + Jii) * ones(1, length(0:20)), 0:20, 'k', 'LineWidth', 3.5);

% Example parameter marker
plot(0.6, 10, '*', 'Color', 'red', 'MarkerSize', 10);

% Axis formatting
xlim([0 2]);
ylim([0 20]);
xlabel('$J_{EI}$', 'Interpreter', 'latex', 'FontSize', 18);
ylabel('$J_{IE}$', 'Interpreter', 'latex', 'FontSize', 18);
set(gca, 'FontSize', 14, 'TickLabelInterpreter', 'latex');

% Region labels
text(0.3,  4,  '$\mathrm{FP}$', 'Interpreter', 'latex', 'FontSize', 26);
text(0.95, 13, '$\mathrm{R}$',  'Interpreter', 'latex', 'FontSize', 26);
text(1.65, 10, '$\mathrm{OI}$', 'Interpreter', 'latex', 'FontSize', 26);

% Annotation for bifurcation line
annotation(figure(1), 'textarrow', ...
    [0.328571428571428 0.240357142857143], ...
    [0.732380952380956 0.660476190476192], ...
    'String', '$J_{IE} = \frac{\bar{J}_d^2}{J_{EI}}$', ...
    'Interpreter', 'latex', ...
    'FontSize', 14);

%% Colorbar
colormap parula;
c = colorbar;

% Remove default colorbar label
c.Label.String = '';
c.Label.Interpreter = 'latex';

% Get colorbar position
cb_pos = c.Position;  % [left bottom width height]

% Create LaTeX-styled title above the colorbar
x = cb_pos(1) + cb_pos(3) / 2;
y = cb_pos(2) + cb_pos(4) + 0.02;

annotation('textbox', [x - 0.05, y, 0.1, 0.03], ...
    'String', '$\mathit{f}\ \mathrm{[Hz]}$', ...
    'Interpreter', 'latex', ...
    'FontSize', 14, ...
    'EdgeColor', 'none', ...
    'HorizontalAlignment', 'center');

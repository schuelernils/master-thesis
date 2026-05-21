clc; clear; close all;

%% === Einstellungen ===
% Cantera (MATLAB Toolbox)
addpath('/Users/nils/Applications/Cantera/matlab/toolbox')

% Eingangszusammensetzung (H2, N2, CO, CO2, CH4, C2H6,  H2O)  [mol]
N_in = [9.33; 0.49; 6.93; 1.03; 22.13; 0.79; 59.29];              % -> Dry Reforming Start: CH4 + CO2
N_tot = sum(N_in);
x_in  = N_in / N_tot;                   % Molanteile am Einlass


% Druck & Temperaturbereich
p     = 51e5;                           % Pa
Tin   = 273.15;                         % K
Tend  = 1600;                           % K
dT    = 10;                             
Trange = Tin:dT:Tend;                   % K

% Species-Namen (Reihenfolge passend zu x_in)
spNames = {'H2','N2','CO','CO2','CH4','C2H6','H2O'};

% Speicher
nT     = numel(Trange);
Xconv  = zeros(1, nT);
x_H2    = zeros(1, nT);
x_N2    = zeros(1, nT);
x_CO    = zeros(1, nT);
x_CO2   = zeros(1, nT);
x_CH4   = zeros(1, nT);
x_C2H6  = zeros(1, nT);
x_H2O   = zeros(1, nT);
Keq_DRM = nan(1, nT);                   % Gleichgewichtskonstante für DRM: CH4 + CO2 = 2CO + 2H2

%% === Cantera: Gas anlegen (GRI30) ===
% Modern: g = Solution('gri30.yaml');
% Kompatibel zu deinem Code:
g = GRI30; 

% Startzustand setzen (nur einmal, dann T im Loop ändern)
set(g, 'T', Trange(1), 'P', p, 'X', sprintf(['H2:%.12g, N2:%.12g, CO:%.12g, CO2:%.12g, ' ...
             'CH4:%.12g, C2H6:%.12g, H2O:%.12g'], ...
             x_in(1), x_in(2), x_in(3), x_in(4), x_in(5), x_in(6), x_in(7)));

w_CH4_in = massFraction(g,'CH4');

%% === Schleife über Temperatur ===
for i = 1:nT
    Tabs = Trange(i);
    set(g, 'T', Tabs, 'P', p);   % nur T ändern, spart etwas Overhead
    equilibrate(g, 'TP');

    % Molanteile im Gleichgewicht
    x_out = moleFractions(g); % Vektor über alle Species
    % Wir picken die relevanten Species:
    x_H2(i)   = moleFraction(g,'H2');
    x_N2(i)   = moleFraction(g,'N2');
    x_CO(i)   = moleFraction(g,'CO');
    x_CO2(i)  = moleFraction(g,'CO2');
    x_CH4(i)  = moleFraction(g,'CH4');
    x_C2H6(i) = moleFraction(g,'C2H6');
    x_H2O(i)  = moleFraction(g,'H2O');

    w_CH4(i) = massFraction(g,'CH4');

    % Umsatz CH4
    if w_CH4_in(1) > 0
        Xconv(i) = (w_CH4_in(1) - w_CH4(i)) / w_CH4_in(1);
    else
        Xconv(i) = 0;
    end
    % Numerisch sauber einklemmen
    Xconv(i) = max(0, min(1, Xconv(i)));
end

%% === Mol-Verhältnis rein/raus ===
T_target_C = 1100;
T_target_K = T_target_C + 273.15;

% Index der Zieltemperatur finden
%[~, idx] = min(abs(Trange - T_target_K));
set(g, 'T', T_target_K, 'P', p, 'X', ...
    sprintf(['H2:%.12g, N2:%.12g, CO:%.12g, CO2:%.12g, ' ...
             'CH4:%.12g, C2H6:%.12g, H2O:%.12g'], ...
             x_in(1), x_in(2), x_in(3), x_in(4), x_in(5), x_in(6), x_in(7)));

equilibrate(g, 'TP');

% Eingang
N_tot_in = sum(N_in);

% Ausgang: normiert → Summe = 1
x_out = [moleFraction(g,'H2');
         moleFraction(g,'N2');
         moleFraction(g,'CO');
         moleFraction(g,'CO2');
         moleFraction(g,'CH4');
         moleFraction(g,'C2H6');
         moleFraction(g,'H2O')];

% Elementmatrix (C, H, O)
% Reihenfolge: H2, N2, CO, CO2, CH4, C2H6, H2O
A = [
    0 0 1 1 1 2 0;   % C
    2 0 0 0 4 6 2;   % H
    0 0 1 2 0 0 1;   % O
];

% Elementmengen inlet
b_in = A * N_in;

% Elementmengen outlet (normiert auf 1 mol)
b_out = A * x_out;

% Skalierungsfaktor (mol_out)
scale = b_in ./ b_out;

% robust: Mittelwert über nicht-NaN
mol_out = mean(scale(~isnan(scale) & isfinite(scale)));

ratio = mol_out / N_tot_in;

fprintf('\nMol-Verhältnis:\n');
fprintf('Mol rein  = %.6f\n', N_tot_in);
fprintf('Mol raus  = %.6f\n', mol_out);
fprintf('Verhältnis (out/in) = %.6f\n', ratio);

%% === Ausgabe der Gleichgewichtszusammensetzung bei 1100 °C ===

fprintf('\nGleichgewichtszusammensetzung bei %.1f °C (%.2f K):\n', T_target_C, T_target_K);
fprintf('H2   = %.6f\n', moleFraction(g,'H2'));
fprintf('N2   = %.6f\n', moleFraction(g,'N2'));
fprintf('CO   = %.6f\n', moleFraction(g,'CO'));
fprintf('CO2  = %.6f\n', moleFraction(g,'CO2'));
fprintf('CH4  = %.6f\n', moleFraction(g,'CH4'));
fprintf('C2H6 = %.6f\n', moleFraction(g,'C2H6'));
fprintf('H2O  = %.6f\n', moleFraction(g,'H2O'));


%% === Plots ===
% 2) Umsatz von CH4
fig = figure(14); clf(fig);
fig.Units = 'centimeters';
fig.Position = [10 10 13 10];

plot(Trange - 273.15, Xconv, 'LineWidth', 1.5);
grid on;
xl = xlabel('Temperature $T$ [$^\circ$C]');
set(xl,'Interpreter','latex');

yl = ylabel('Conversion $X_{CH_4}$ [-]');
set(yl,'Interpreter','latex');

ylim([0 1]);

print(fig, 'CH4_conversion.emf', '-dpdf');   % Für PowerPoint

% 3) Molanteile ausgewählter Komponenten
fig = figure(15); clf(fig);
fig.Units = 'centimeters';
fig.Position = [10 10 13 10];
hold on;

p1 = plot(Trange - 273.15, x_H2,   'LineWidth', 1.8, 'LineStyle', '-');
p3 = plot(Trange - 273.15, x_CO,   'LineWidth', 1.8, 'LineStyle', '-.');
p4 = plot(Trange - 273.15, x_CO2,  'LineWidth', 1.8, 'LineStyle', ':');
p5 = plot(Trange - 273.15, x_CH4,  'LineWidth', 1.8, 'LineStyle', '-');
p7 = plot(Trange - 273.15, x_H2O,  'LineWidth', 1.8, 'LineStyle', '-.');
grid on; 
xl = xlabel('Temperature $T$ [$^\circ$C]');
set(xl,'Interpreter','latex');
yl = ylabel('Mole Fraction [-]');
set(yl,'Interpreter','latex');
axis([Tin-273.15, Tend-273.15, 0, 0.6]);
legend([p1 p3 p4 p5 p7], {'H_2','CO','CO_2','CH_4','H_2O'}, 'Location','best');
hold off;

print(fig, 'mole_fraction.emf', '-dpdf');   % Für PowerPoint

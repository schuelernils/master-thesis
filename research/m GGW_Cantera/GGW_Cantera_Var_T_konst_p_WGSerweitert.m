clc; clear; close all;

%% === Einstellungen ===
% Cantera (MATLAB Toolbox)
addpath('/Users/nils/Applications/Cantera/matlab/toolbox')

% Eingangszusammensetzung (H2, N2, CO, CO2, CH4, C2H6,  H2O)  [mol]
N_in = [12.96; 0.96; 10.09; 1.39; 20.10; 0.71; 53.79];              % -> Dry Reforming Start: CH4 + CO2
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



%% =========================================================================
%%   WGS-REAKTOR  –  Einzelreaktion:  CO + H2O  ⇌  CO2 + H2
%%
%%   Dieser Block schließt direkt an GGW_Cantera_Var_T_konst_p.m an.
%%   Voraussetzung: g, x_out, mol_out, x_in, p sind noch im Workspace.
%%
%%   Methode:
%%     Keq(T) = exp(–ΔG°(T)/RT)  aus Cantera-Standardzustandsdaten
%%     Gleichgewicht: analytische Lösung der quadratischen Gleichung
%%     (kein Solver, keine Vollgleichgewicht-Berechnung mit equilibrate)
%% =========================================================================

%% === 1) Feed des WGS-Reaktors ============================================
%  Einlass = Reformer-Gleichgewichtsgas bei 1100 °C  (aus vorherigem Block)
%  x_out  : Molanteilvektor [H2, N2, CO, CO2, CH4, C2H6, H2O]
%  mol_out: absolute Gesamtmolzahl (skaliert auf Einlass-Elementbilanz)

N_WGS_in = x_out * mol_out;   % mol jeder Spezies

n_H2_in   = N_WGS_in(1);
n_N2_in   = N_WGS_in(2);
n_CO_in   = N_WGS_in(3);
n_CO2_in  = N_WGS_in(4);
n_CH4_in  = N_WGS_in(5);
n_C2H6_in = N_WGS_in(6);
n_H2O_in  = N_WGS_in(7);
N_WGS_tot = sum(N_WGS_in);    % mol gesamt (konstant, da Δn_WGS = 0)

fprintf('\n=== WGS-Reaktor  (Einzelreaktion: CO + H2O ⇌ CO2 + H2) ===\n');

%% === 2) Temperaturbereich WGS ============================================
%  Typischer Betriebsbereich: 200–600 °C
%  (Hochtemperatur-WGS ~300–500 °C, Niedertemperatur-WGS ~200–300 °C)

T_WGS_range_C = 200:5:600;             % °C
T_WGS_range   = T_WGS_range_C + 273.15; % K
nT_WGS        = numel(T_WGS_range);

% Speicherfelder
Keq_WGS  = zeros(1, nT_WGS);
xi_eq    = zeros(1, nT_WGS);
X_CO     = zeros(1, nT_WGS);
x_H2_W   = zeros(1, nT_WGS);
x_CO_W   = zeros(1, nT_WGS);
x_CO2_W  = zeros(1, nT_WGS);
x_H2O_W  = zeros(1, nT_WGS);

%% === 3) Schleife: Keq + Gleichgewichtsextent =============================
%
%  Standardzustandseigenschaften (reine T-Funktionen für ideale Gase):
%    enthalpies_RT(g) → h°_i / (R·T)   [dimensionslos, Vektor über alle Spezies]
%    entropies_R(g)   → s°_i / R        [dimensionslos, Vektor über alle Spezies]
%
%  Gibbs-Energie:  g°_i/(RT) = h°_i/(RT) – s°_i/R
%  Reaktion:       ΔG°/(RT)  = Σ νᵢ · g°_i/(RT)
%  Gleichgewicht:  Keq       = exp(–ΔG°/(RT))
%
%  Gleichgewicht (analytisch):
%    Keq · (n_CO – ξ)(n_H2O – ξ) = (n_CO2 + ξ)(n_H2 + ξ)
%    ⟹ (Keq–1)ξ² – [Keq(n_CO+n_H2O) + (n_CO2+n_H2)]ξ
%       + [Keq·n_CO·n_H2O – n_CO2·n_H2] = 0

% Spezies-Indizes einmalig holen
iH2  = speciesIndex(g, 'H2');
iCO  = speciesIndex(g, 'CO');
iCO2 = speciesIndex(g, 'CO2');
iH2O = speciesIndex(g, 'H2O');

xi_max = min(n_CO_in, n_H2O_in);   % physikalische Obergrenze für ξ

for i = 1:nT_WGS
    T_i = T_WGS_range(i);

    % --- Cantera: nur T und P setzen (X irrelevant für Standardzustand) ---
    set(g, 'T', T_i, 'P', p);

    % --- Dimensionslose Standardzustands-Gibbs-Energien ---
    h_RT = enthalpies_RT(g);     % h°_i / (R·T)
    s_R  = entropies_R(g);       % s°_i / R
    g_RT = h_RT - s_R;           % g°_i / (R·T)

    % ΔG°/(RT) = g°_CO2 + g°_H2 – g°_CO – g°_H2O   (alle durch RT)
    dG_RT = g_RT(iCO2) + g_RT(iH2) - g_RT(iCO) - g_RT(iH2O);
    Keq_i = exp(-dG_RT);
    Keq_WGS(i) = Keq_i;

    % --- Quadratische Gleichung für ξ ---
    %   a·ξ² + b·ξ + c = 0
    a_q =  Keq_i - 1;
    b_q = -(Keq_i*(n_CO_in + n_H2O_in) + (n_CO2_in + n_H2_in));
    c_q =  Keq_i*n_CO_in*n_H2O_in - n_CO2_in*n_H2_in;

    if abs(a_q) < 1e-12        % Sonderfall Keq ≈ 1 → lineare Gleichung
        xi_sol = -c_q / b_q;
    else
        disc = b_q^2 - 4*a_q*c_q;
        if disc < 0
            disc = 0;          % numerische Toleranz
        end
        xi1 = (-b_q + sqrt(disc)) / (2*a_q);
        xi2 = (-b_q - sqrt(disc)) / (2*a_q);

        % Physikalisch sinnvolle Lösung wählen: 0 ≤ ξ ≤ min(n_CO, n_H2O)
        candidates = [xi1, xi2];
        mask = candidates >= -1e-10 & candidates <= xi_max + 1e-10;
        if any(mask)
            xi_sol = min(candidates(mask));   % kleinere (konservativere) Wurzel
        else
            xi_sol = 0;
        end
    end

    % Numerisch einklemmen
    xi_sol   = max(0, min(xi_sol, xi_max));
    xi_eq(i) = xi_sol;

    % --- Gleichgewichts-Molmengen ---
    n_H2_eq  = n_H2_in  + xi_sol;
    n_CO_eq  = n_CO_in  - xi_sol;
    n_CO2_eq = n_CO2_in + xi_sol;
    n_H2O_eq = n_H2O_in - xi_sol;
    % N_tot konstant (Δn = 0 für WGS)

    x_H2_W(i)  = n_H2_eq  / N_WGS_tot;
    x_CO_W(i)  = n_CO_eq  / N_WGS_tot;
    x_CO2_W(i) = n_CO2_eq / N_WGS_tot;
    x_H2O_W(i) = n_H2O_eq / N_WGS_tot;

    if n_CO_in > 0
        X_CO(i) = xi_sol / n_CO_in;
    end
end

%% === 4) Textausgabe bei Zieltemperatur ===================================
T_WGS_out_C = 440;
[~, idx_out] = min(abs(T_WGS_range_C - T_WGS_out_C));

xi_out = xi_eq(idx_out);

% Absolute Molmengen aller Spezies am Auslass
n_out_H2   = n_H2_in   + xi_out;
n_out_N2   = n_N2_in;
n_out_CO   = n_CO_in   - xi_out;
n_out_CO2  = n_CO2_in  + xi_out;
n_out_CH4  = n_CH4_in;
n_out_C2H6 = n_C2H6_in;
n_out_H2O  = n_H2O_in  - xi_out;
% N_tot unveraendert (delta_n = 0)

% Molanteile aller Spezies
x_out_all = [n_out_H2; n_out_N2; n_out_CO; n_out_CO2; ...
             n_out_CH4; n_out_C2H6; n_out_H2O] / N_WGS_tot;

fprintf('\n--- WGS-Gleichgewicht bei %.0f °C (%.2f K) ---\n', ...
        T_WGS_out_C, T_WGS_range(idx_out));
fprintf('  Keq             = %10.4f\n',        Keq_WGS(idx_out));
fprintf('  xi_eq           = %10.6f mol\n',    xi_out);
fprintf('  X_CO            = %10.4f  (%.2f %%)\n', X_CO(idx_out), X_CO(idx_out)*100);

fprintf('\nGleichgewichtszusammensetzung:\n');
fprintf('H2   = %.6f\n', x_out_all(1));
fprintf('N2   = %.6f\n', x_out_all(2));
fprintf('CO   = %.6f\n', x_out_all(3));
fprintf('CO2  = %.6f\n', x_out_all(4));
fprintf('CH4  = %.6f\n', x_out_all(5));
fprintf('C2H6 = %.6f\n', x_out_all(6));
fprintf('H2O  = %.6f\n', x_out_all(7));


%% === 5) Plots ============================================================

% --- Plot B: CO-Umsatz über T ---
fig_wgs2 = figure(22); clf(fig_wgs2);
fig_wgs2.Units    = 'centimeters';
fig_wgs2.Position = [10 10 13 10];
plot(T_WGS_range_C, X_CO, 'LineWidth', 1.8);
grid on;
xl = xlabel('Temperature $T$ [$^\circ$C]'); set(xl,'Interpreter','latex');
yl = ylabel('Conversion $X_\mathrm{CO}$ [-]'); set(yl,'Interpreter','latex');
ylim([0 1]);
title('WGS: CO-Umsatz im thermodynamischen GGW');
print(fig_wgs2, 'WGS_CO_conversion.pdf', '-dpdf');

% --- Plot C: Molanteile über T ---
fig_wgs3 = figure(23); clf(fig_wgs3);
fig_wgs3.Units    = 'centimeters';
fig_wgs3.Position = [10 10 13 10];
hold on;
pw1 = plot(T_WGS_range_C, x_H2_W,  'LineWidth', 1.8, 'LineStyle', '-');
pw2 = plot(T_WGS_range_C, x_CO_W,  'LineWidth', 1.8, 'LineStyle', '-.');
pw3 = plot(T_WGS_range_C, x_CO2_W, 'LineWidth', 1.8, 'LineStyle', ':');
pw4 = plot(T_WGS_range_C, x_H2O_W, 'LineWidth', 1.8, 'LineStyle', '--');
grid on;
xl = xlabel('Temperature $T$ [$^\circ$C]'); set(xl,'Interpreter','latex');
yl = ylabel('Mole Fraction [-]');            set(yl,'Interpreter','latex');
legend([pw1 pw2 pw3 pw4], {'H_2','CO','CO_2','H_2O'}, 'Location','best');
title('WGS: Molanteile im thermodynamischen GGW');
hold off;
print(fig_wgs3, 'WGS_mole_fractions.pdf', '-dpdf');
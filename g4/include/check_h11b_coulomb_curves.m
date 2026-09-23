% check_h11b_coulomb_curves.m
% MATLAB check script for the latest H11B GEANT4 Coulomb-penetrability version.
%
% Put this file in the GEANT4 project root directory, i.e. the directory
% containing:
%   include/CoulombPenetrabilityData.hh
%
% It reads the precomputed log(P_L) tables from CoulombPenetrabilityData.hh,
% uses log(P) vs log(E) interpolation, and plots:
%   1) 165-keV p+11B resonance cross section
%   2) 675-keV p+11B resonance cross section
%   3) 8Be(2+) excitation-energy sampling weight
%
% Units:
%   Energies in keV unless otherwise stated.
%   Cross sections are plotted in barn and include the same 1e5 sampling
%   scale factor used in H11BCrossSection.cc.

clear; clc; close all;

%% ------------------------------------------------------------------------
%  Locate and read the C++ penetrability table
% -------------------------------------------------------------------------
tableFile = fullfile('include', 'CoulombPenetrabilityData.hh');

if ~isfile(tableFile)
    error(['Cannot find %s. Please run this MATLAB script from the GEANT4 ', ...
           'project root directory.'], tableFile);
end

txt = fileread(tableFile);

P11B_E        = extract_cpp_array(txt, 'P11B_energy_keV');
P11B_L0_logP  = extract_cpp_array(txt, 'P11B_L0_logP');
P11B_L1_logP  = extract_cpp_array(txt, 'P11B_L1_logP');

AA_E          = extract_cpp_array(txt, 'AlphaAlpha_energy_keV');
AA_L2_logP    = extract_cpp_array(txt, 'AlphaAlpha_L2_logP');

%% ------------------------------------------------------------------------
%  Constants, consistent with H11BCrossSection.cc
% -------------------------------------------------------------------------
amu_c2_keV    = 931494.10242;
hbarc_keV_fm  = 197326.9804;

projectA = 1.0;
targetA  = 11.0;
mu_p11B_keV = projectA * targetA / (projectA + targetA) * amu_c2_keV;

% Same artificial sampling scale still present in the GEANT4 code
CrossSectionSamplingScale = 1.0;

% If you have changed these in Constants.hh, change them here too.
H11B165SpinStatFactor = 5.0/8.0;
H11B675SpinStatFactor = 5.0/8.0;

% To check the strictly physical Breit-Wigner normalization instead, use:
% H11B165SpinStatFactor = 5.0/8.0;
% H11B675SpinStatFactor = 5.0/8.0;

%% ------------------------------------------------------------------------
%  p + 11B resonance parameters from Constants.hh
% -------------------------------------------------------------------------
% 165-keV resonance, Ecm value
Er165        = 148.3;
GammaTot165_R = 5.3;
GammaP165_R = 0.0215;
GammaA165   = 0.26 + 5.0;  % alpha0 + alpha1, gamma channel not included

% 675-keV resonance, Ecm value
Er675        = 618.75;
GammaP675_R = 150.0;
GammaA675   = 150.0;       % alpha0 = 0, alpha1 = 150 in default code

%% ------------------------------------------------------------------------
%  8Be(2+) parameters from Constants.hh
% -------------------------------------------------------------------------
Ex8Be2Plus = 3030.0;       % keV, relative to 8Be ground state
Gamma8Be2Plus_R = 1513.0; % keV
Ex8BeGroundAbove2Alpha = 91.84; % keV
Eaa8Be2Plus = Ex8Be2Plus + Ex8BeGroundAbove2Alpha; % keV, relative to 2-alpha threshold

%% ------------------------------------------------------------------------
%  Energy grids
% -------------------------------------------------------------------------
E_p11B = linspace(1.0, 1200.0, 4000);     % Ecm for p+11B
Ex8Be  = linspace(0.0, 8000.0, 4000);     % 8Be excitation relative to 8Be ground
Eaa    = Ex8Be + Ex8BeGroundAbove2Alpha;  % alpha-alpha relative energy

%% ------------------------------------------------------------------------
%  Compute 165 and 675 cross sections using full Coulomb penetrability ratio
% -------------------------------------------------------------------------
% 165: p + 11B entrance channel, L = 1
P_ratio_165 = penetrability_ratio(E_p11B, Er165, P11B_E, P11B_L1_logP);
GammaP165_E = GammaP165_R .* P_ratio_165;

% Preserve the total width at the resonance energy by carrying any missing
% width as an energy-independent remainder.
GammaMissing165 = max(0.0, GammaTot165_R - GammaP165_R - GammaA165);
GammaTot165_E = GammaP165_E + GammaA165 + GammaMissing165;

sigma165 = sigma_breit_wigner(E_p11B, Er165, GammaP165_E, GammaA165, ...
                              GammaTot165_E, H11B165SpinStatFactor, ...
                              mu_p11B_keV, hbarc_keV_fm, CrossSectionSamplingScale);

% Original code window for 165 was ecm < 1 or ecm > 400 => 0
sigma165(E_p11B < 1.0 | E_p11B > 400.0) = 0.0;

% 675: p + 11B entrance channel, L = 0
P_ratio_675 = penetrability_ratio(E_p11B, Er675, P11B_E, P11B_L0_logP);
GammaP675_E = GammaP675_R .* P_ratio_675;
GammaTot675_E = GammaP675_E + GammaA675;  % no gamma channel, no missing width

sigma675 = sigma_breit_wigner(E_p11B, Er675, GammaP675_E, GammaA675, ...
                              GammaTot675_E, H11B675SpinStatFactor, ...
                              mu_p11B_keV, hbarc_keV_fm, CrossSectionSamplingScale);

% Latest modified version removes the old 300-keV lower cutoff.
% Keep only the broad safety range used in the modified code.
sigma675(E_p11B < 1.0 | E_p11B > 3500.0) = 0.0;

%% ------------------------------------------------------------------------
%  Compute 8Be(2+) sampling weight using alpha-alpha L = 2 penetrability
% -------------------------------------------------------------------------
P_ratio_aa_L2 = penetrability_ratio(Eaa, Eaa8Be2Plus, AA_E, AA_L2_logP);
Gamma8Be_E = Gamma8Be2Plus_R .* P_ratio_aa_L2;

delta8Be = Eaa - Eaa8Be2Plus;
Weight8Be = Gamma8Be_E ./ (delta8Be.^2 + Gamma8Be_E.^2 ./ 4.0);

% Same protection as the C++ code
Weight8Be(Eaa <= 0) = 0.0;

%% ------------------------------------------------------------------------
%  Plot: three curves
% -------------------------------------------------------------------------
figure('Color','w','Position',[100,100,950,780]);

subplot(3,1,1);
plot(E_p11B, sigma165, 'LineWidth', 1.5);
xlim([100, 400]);
xlabel('E_{cm} (keV)');
ylabel('\sigma_{165} (barn)');
title('165-keV resonance: p + ^{11}B entrance, L = 1, full Coulomb P_L');
grid on;

subplot(3,1,2);
plot(E_p11B, sigma675, 'LineWidth', 1.5);
xlim([100, 1200]);
xlabel('E_{cm} (keV)');
ylabel('\sigma_{675} (barn)');
title('675-keV resonance: p + ^{11}B entrance, L = 0, full Coulomb P_L');
grid on;

subplot(3,1,3);
plot(Ex8Be, Weight8Be ./ max(Weight8Be), 'LineWidth', 1.5);
xlim([0, 8000]);
xlabel('E_x(^{8}Be) relative to ^{8}Be ground state (keV)');
ylabel('Normalized weight');
title('^{8}Be(2^+) sampling weight: \alpha+\alpha, L = 2, full Coulomb P_L');
grid on;

exportgraphics(gcf, 'h11b_three_coulomb_curves.png', 'Resolution', 200);

%% ------------------------------------------------------------------------
%  Extra plot: 165 and 675 together on log scale
% -------------------------------------------------------------------------
figure('Color','w','Position',[120,120,900,520]);
semilogy(E_p11B, sigma165, 'LineWidth', 1.5); hold on;
semilogy(E_p11B, sigma675, 'LineWidth', 1.5);
xlim([100, 1200]);
ylim([1e-3, max([sigma165, sigma675])*2]);
xlabel('E_{cm} (keV)');
ylabel('Cross section (barn, including code \times10^5 scale)');
title('p + ^{11}B cross sections with full Coulomb penetrability');
legend('165 keV resonance, L=1', '675 keV resonance, L=0', 'Location','best');
grid on;
exportgraphics(gcf, 'h11b_165_675_coulomb_log.png', 'Resolution', 200);

%% ------------------------------------------------------------------------
%  Print useful diagnostic values
% -------------------------------------------------------------------------
fprintf('\nDiagnostic values, including code sampling scale x1e5:\n');

print_value('sigma165 at 148.3 keV', interp1(E_p11B, sigma165, 148.3), 'barn');
print_value('sigma675 at 148.3 keV', interp1(E_p11B, sigma675, 148.3), 'barn');
print_value('sigma675 at 618.75 keV', interp1(E_p11B, sigma675, 618.75), 'barn');

[~, iMax8] = max(Weight8Be);
print_value('8Be weight maximum Ex', Ex8Be(iMax8), 'keV');

fprintf('\nOutput figures:\n');
fprintf('  h11b_three_coulomb_curves.png\n');
fprintf('  h11b_165_675_coulomb_log.png\n');

%% ========================================================================
%  Local functions
% ========================================================================

function arr = extract_cpp_array(txt, name)
    % Extracts:
    %   static const double name[] = { ... };
    pattern = ['static\s+const\s+double\s+', name, '\s*\[\]\s*=\s*\{(.*?)\};'];
    token = regexp(txt, pattern, 'tokens', 'once');
    if isempty(token)
        error('Could not find array "%s" in CoulombPenetrabilityData.hh.', name);
    end

    body = token{1};
    body = regexprep(body, '//[^\n\r]*', ' ');
    nums = regexp(body, '[-+]?\d*\.?\d+(?:[eE][-+]?\d+)?', 'match');
    arr = str2double(nums);
end

function logP = interp_logP(E, E_tab, logP_tab)
    % Interpolate log(P) as a function of log(E), matching the C++ table logic.
    E_safe = max(E, min(E_tab));
    E_safe = min(E_safe, max(E_tab));
    logP = interp1(log(E_tab), logP_tab, log(E_safe), 'linear', 'extrap');
end

function ratio = penetrability_ratio(E, Er, E_tab, logP_tab)
    logP_E  = interp_logP(E,  E_tab, logP_tab);
    logP_Er = interp_logP(Er, E_tab, logP_tab);
    ratio = exp(logP_E - logP_Er);
end

function sigma_barn = sigma_breit_wigner(E, Er, GammaIn_E, GammaOut, GammaTot_E, ...
                                         spinFactor, mu_keV, hbarc_keV_fm, scaleFactor)
    pi_over_k2_fm2 = pi .* hbarc_keV_fm.^2 ./ (2.0 .* mu_keV .* E);
    denominator = (E - Er).^2 + GammaTot_E.^2 ./ 4.0;
    bw = spinFactor .* GammaIn_E .* GammaOut ./ denominator;

    sigma_cm2 = pi_over_k2_fm2 .* bw .* 1e-26; % fm^2 -> cm^2
    sigma_cm2 = sigma_cm2 .* scaleFactor;
    sigma_barn = sigma_cm2 ./ 1e-24;           % cm^2 -> barn
end

function print_value(label, value, unit)
    fprintf('  %-28s = %.6g %s\n', label, value, unit);
end

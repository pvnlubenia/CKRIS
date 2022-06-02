% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
%      Power Law Approximation of Hill-Type Kinetics                      %
%                                                                         %
% OUTPUT: Returns the parameter values of the power law approximation of  %
%    Hill-type kinetic reactions (v29 and v33) of Brännmark et al's       %
%    insulin resistance signaling model [1]. We follow the procedure      %
%    detailed in Appendix B of Fortun et al [2] to approximate a non-     %
%    power law kinetic equation with power law kinetics.                  %
%                                                                         %
% References:                                                             %
%    [1] Brännmark C, Nyman E, Fagerholm S, Bergenholm L, Ekstrand E,     %
%           Cedersund G, Stralfors P (2013) Insulin signaling in type 2   %
%           diabetes: experimental and modeling analyses reveal           %
%           mechanisms of insulin resistance in human adipocytes. J Biol  %
%           Chem 288(14):9867--9880.                                      %
%           https://doi.org/10.1074/jbc.m112.432062                       %
%    [2] Fortun N, Mendoza E, Razon L, Lao A (2018) A deficiency-one      %
%           algorithm for power-law kinetic systems with reactant-        %
%           determined interactions. J Math Chem 56:2929--2962.           %
%           https://doi.org/10.1007/s10910-018-0925-2                     %
%                                                                         %
% Created: 29 May 2022                                                    %
% Last Modified: 2 June 2022                                              %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

% Variables as symbols
syms k29 x28 n6 km6 x34
syms k33 x36 n9 km9 x31

% Original reactions
v29 = k29 * ((x28^n6)/(km6^n6 + x28^n6)) * x34;
v33 = k33 * x36 * ((x31^n9)/(km9^n9 + x31^n9));

% Partial derivatives
pv29x28 = diff(v29, x28);
pv29x34 = diff(v29, x34);

pv33x36 = diff(v33, x36);
pv33x31 = diff(v33, x31);

% Constants
%    - Based on Matlab file from BioModels:
% ﻿        https://www.ebi.ac.uk/biomodels/BIOMD0000000448)
k29 = 36.93;
n6 = 2.137;
km6 = 30.54;
k33 = 0.1298;
n9 = 0.9855;
km9 = 5873.0;

% Equilibria
%    - Run the Matlab file from BioModels for tspan = [0:0.01:1000]
x28 = 37.9227;
x34 = 29.576;
x36 = 96.0473;
x31 = 78.7914;

% Power law approximations:
%    v29 = k29p x28^p29 x34^q29
%    v33 = k33p x36^p33 x31^q33

% Solve for the parameter values
p29 = double(subs(pv29x28 * (x28 / v29)))
q29 = double(subs(pv29x34 * (x34 / v29)))
k29p = double(subs(v29 / ((x28^p29) * (x34^q29))))
p33 = double(subs(pv33x36 * (x36 / v33)))
q33 = double(subs(pv33x31 * (x31 / v33)))
k33p = double(subs(v33 / ((x36^p33) * (x31^q33))))
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                         %
%      Braennmark et al Insulin Resistance Signaling (INRES) Model        %
%                                                                         %
% OUTPUT: Returns the final values of the state variables of Braennmark   %
%    et al's insulin resistance signaling model over a period of time. It %
%    also returns the graph of 4 species (AS160, AS160p, S6K, S6Kp) and   %
%    their approximation when their Hill-type reactions are approximated  %
%    using power law kinetics.                                            %
%                                                                         %
% Reference: ﻿Braennmark C, Nyman E, Fagerholm S, Bergenholm L,            %
%    Ekstrand E, Cedersund G, Stralfors P (2013) Insulin signaling in     %
%    type 2 diabetes: experimental and modeling analyses reveal           %
%    mechanisms of insulin resistance in human adipocytes. J Biol Chem    %
%    288(14):9867--9880. https://doi.org/10.1074/jbc.m112.432062          %
%                                                                         %
% Created: 2 June 2022                                                    %
% Last Modified: 2 June 2022                                              %
%                                                                         %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %

%
% Main function
%   - Contains the initial values, time series data, and graph
%

function inres()
	
	% Initial values
	%    - Variables in the comments are the ones used in CKRIS research
	x0 = zeros(31,1);             % Originally just 27 species
	x0(1) = 99.8737104842408;     % X2  = IR
	x0(2) = 0.00186253217635894;  % X4  = IRp
	x0(3) = 0.0;                  % X3  = IRins
	x0(4) = 0.0188430465801578;   % X7  = IRip
	x0(5) = 0.105583925473107;    % X6  = IRi
	x0(6) = 82.9671997523599;     % X9  = IRS1
	x0(7) = 0.00119481841136737;  % X10 = IRS1p
	x0(8) = 0.327454355438396;    % X22 = IRS1p307
	x0(9) = 16.7041510257561;     % X23 = IRS1307
	x0(10) = 99.9983336594667;    % X24 = X
	x0(11) = 0.00166634053318549; % X25 = Xp
	x0(12) = 68.1806649661901;    % X26 = PKB
	x0(13) = 13.2964849666951;    % X27 = PKB308p
	x0(14) = 16.8171941560617;    % X28 = PKB473p
	x0(15) = 1.70566051030056;    % X29 = PKB308p473p
	x0(16) = 86.5002472240273;    % X30 = mTORC1
	x0(17) = 13.4997527759726;    % X31 = mTORC1a
	x0(18) = 99.8478148461591;    % X32 = mTORC2
	x0(19) = 0.152185153840861;   % X33 = mTORC2a
	x0(20) = 83.8141018591099;    % X34 = AS160
	x0(21) = 16.1858981408903;    % X35 = AS160p
	x0(22) = 26.523878746229;     % X21 = GLUT4m
	x0(23) = 73.476121253771;     % X20 = GLUT4
	x0(24) = 99.2731987219547;    % X36 = S6K
	x0(25) = 0.72680127804522;    % X37 = S6Kp
	x0(26) = 92.7596420796038;    % X38 = S6
	x0(27) = 7.24035792039603;    % X39 = S6p
	% For approximation of Hill-type kinetics with power law
	x0(28) = 83.8141018591099;    % X34 ~ AS160 approx
	x0(29) = 16.1858981408903;    % X35 ~ AS160p approx
	x0(30) = 99.2731987219547;    % X36 ~ S6K approx
	x0(31) = 0.72680127804522;    % X37 ~ S6Kp approx
	
	% Time series
	tspan = [0:0.01:100]; % Use up to 1000 for equilibrium
	opts = odeset('AbsTol',1e-3);
	[t,x] = ode23tb(@f,tspan,x0,opts);    
	x(size(x,1),:) % Outputs final values of the time series
	
	% Plot (all)
	%    - Uncomment to show plot of all species
	%  	plot(t,x);
	
	% Plot (approximated species only)
	%   - Comment the lines below when uncommenting "Plot (all)"
	plot(t,x(:,20),'-', t,x(:,28),'--', t,x(:,21),'-', t,x(:,29),'--', t,x(:,24),'-', t,x(:,30),'--', t,x(:,25),'-', t,x(:,31),'--');
	legend('X34', 'X34 approx', 'X35', 'X35 approx', 'X36', 'X36 approx', 'X37', 'X37 approx', 'location', 'east');
	legend boxoff
	xlabel('time (min)');
	ylabel('mol/L');
	set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 0.75 0.75]);
end





%
% Function for ODEs
%

function xdot = f(t,x)
	
	% Parameters
	%    - Comments refer to variables used in Braennmark et al
	k1 = 0.6331*10;    % k1a
	k2 = 0.03683;      % k1basal
	k3 = 0.8768;       % k1c
	k4 = 31.01;        % k1d
	k5 = 1840.0;       % k1f
	k6 = 1944.0;       % k1g
	k7 = 0.5471;       % k1r
	k8 = 3.227;        % k2a
	k9 = 3424.0;       % k2b
	k10 = 5759.0;      % k2c
	k11 = 280.8;       % k2d
	k12 = 2.913;       % k2f
	k13 = 0.04228;     % k2basal
	k14 = 0.2671;      % k2g
	k15 = 0.001377;    % k3a
	k16 = 0.09876;     % k3b
	k17 = 5790.0;      % k4a
	k18 = 34.8;        % k4b
	k19 = 4.456;       % k4c
	k20 = 42.84;       % k4e
	k21 = 143.6;       % k4f
	k22 = 0.5361;      % k4h
	k23 = 1.842;       % k5a1
	k24 = 0.05506;     % k5a2
	k25 = 24.83;       % k5b
	k26 = 0.08575;     % k5c
	k27 = 1.06;        % k5d
	k28 = 2.652;       % k6f1
	k29 = 36.93;       % k6f2
	k30 = 65.18;       % k6b
	k31 = 50.98;       % k7f
	k32 = 2286.0;      % k7b
	k33 = 0.1298;      % k9f1
	k34 = 0.04441;     % k9b1
	k35 = 31.0;        % k9b2
	k36 = 3.329;       % k9f2
	km6 = 30.54;       % km6
	km9 = 5873.0;      % km9
	n6 = 2.137;        % n6
	n9 = 0.9855;       % n9
	% For approximation of Hill-type kinetics with power law
	k29p = 1.1265;     % k9'
	k33p = 2.6265e-05; % k33'
	p29 = 0.8256;      % p29
	q33 = 0.9716;      % q33
	
	% Reactions
	%    - Comments refer to variables used in Braennmark et al
	v1   = k1 * x(1);                                    % v1a
	v2   = k2 * x(1);                                    % v1basal
	v3   = k3 * x(3);                                    % v1c
	v4   = k4 * x(2);                                    % v1d
	v5   = k5 * x(4) * x(11);                            % v1e
	v6   = k6 * x(2);                                    % v1g
	v7   = k7 * x(5);                                    % v1r
	v8   = k8 * x(4) * x(6);                             % v2a
	v9   = k9 * x(7);                                    % v2b
	v10  = k10 * x(7) * x(17);                           % v2c
	v11  = k11 * x(8);                                   % v2d
	v12  = k12 * x(8);                                   % v2f
	v13  = k13 * x(6);                                   % v2basal
	v14  = k14 * x(9);                                   % v2g
	v15  = k15 * x(7) * x(10);                           % v3a
	v16  = k16 * x(11);                                  % v3b
	v17  = k17 * x(7) * x(12);                           % v4a
	v18  = k18 * x(13);                                  % v4b
	v19  = k19 * x(13) * x(19);                          % v4c
	v20  = k20 * x(8) * x(14);                           % v4e
	v21  = k21 * x(15);                                  % v4f
	v22  = k22 * x(14);                                  % v4h
	v23  = k23 * x(15) * x(16);                          % v5a
	v24  = k24 * x(13) * x(16);                          % v5a
	v25  = k25 * x(17);                                  % v5b
	v26  = k26 * x(4) * x(18);                           % v5c
	v27  = k27 * x(19);                                  % v5d
	v28  = k28 * x(15) * x(20);                          % v6f1
	v29  = k29 * (x(14)^n6/(km6^n6 + x(14)^n6)) * x(20); % v6f1
	v30  = k30 * x(21);                                  % v6b1
	v31  = k31 * x(21) * x(23);                          % v7f
	v32  = k32 * x(22);                                  % v7b
	v33  = k33 * x(24) * (x(17)^n9/(km9^n9 + x(17)^n9)); % v9f1
	v34  = k34 * x(25);                                  % v9b1
	v35  = k35 * x(27);                                  % v9b2
	v36  = k36 * x(25) * x(26);                          % v9f2
	% For approximation of Hill-type kinetics with power law
	v29p = k29p * x(14)^p29 * x(28);                     % v29p
	v33p = k33p * x(30) * x(17)^q33;                     % v33p
	
	% Actual ODEs
	xdot = zeros(31,1); % Originally just 27 species
	xdot(1)  = -v1 - v2 + v6 + v7;
	xdot(2)  = v2 + v3 - v4 - v6;
	xdot(3)  = v1 - v3;
	xdot(4)  = v4 - v5;
	xdot(5)  = v5 - v7;
	xdot(6)  = -v8 + v9 - v13 + v14;
	xdot(7)  = v8 - v9 - v10 + v11;
	xdot(8)  = v10 - v11 - v12;
	xdot(9)  = v12 + v13 - v14;
	xdot(10) = -v15 + v16;
	xdot(11) = v15 - v16;
	xdot(12) = -v17 + v18 + v22;
	xdot(13) = v17 - v18 - v19;
	xdot(14) = -v20 + v21 - v22;
	xdot(15) = v19 + v20 - v21;
	xdot(16) = -v23 - v24 + v25;
	xdot(17) = v23 + v24 - v25;
	xdot(18) = -v26 + v27;
	xdot(19) = v26 - v27;
	xdot(20) = -v28 - v29 + v30;
	xdot(21) = v28 + v29 - v30;
	xdot(22) = v31 - v32;
	xdot(23) = -v31 + v32;
	xdot(24) = -v33 + v34;
	xdot(25) = v33 - v34;
	xdot(26) = -v36 + v35;
	xdot(27) = v36 - v35;
	% For approximation of Hill-type kinetics with power law
	xdot(28) = -v28 - v29p + v30;
	xdot(29) = v28 + v29p - v30;
	xdot(30) = -v33p + v34;
	xdot(31) = v33p - v34;
end
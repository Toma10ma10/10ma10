%% beam_cc_viv_forced.m
clear; clc;

E   = 210e9; I = 1.0e-8; rho = 7850; A = 1.0e-4; L = 1.0;
nNodes = 100;

[M, K] = assembleBeamEB(E,I,rho,A,L,nNodes);

alpha = 1e-3; beta = 1e-6;
C = alpha*M + beta*K;

fixed = [1 2, 2*nNodes-1, 2*nNodes];
[Mr, Cr, Kr, free] = reduceByBC(M,C,K,fixed);

tF = linspace(0, 20, 4001).';
fviv_hz = 0.5;
Famp = 10.0;
Fviv_w = zeros(length(tF), nNodes);
Fviv_w(:, :) = Famp * sin(2*pi*fviv_hz*tF);

Ffull = expandNodalWForceToFullDOF(Fviv_w);
Fr = Ffull(:, free);

tspan = [tF(1) tF(end)];
x0 = zeros(length(free),1);
v0 = zeros(length(free),1);
z0 = [x0; v0];

Fr_fun = @(t) interp1(tF, Fr, t, 'linear', 'extrap').';

ode = @(t,z) stateSpaceRHS(t, z, Mr, Cr, Kr, Fr_fun);
opts = odeset('RelTol',1e-6,'AbsTol',1e-9);
[tSol, zSol] = ode45(ode, tspan, z0, opts);

x_r = zSol(:, 1:length(free));
ndof = 2*nNodes;
x_full = zeros(length(tSol), ndof);
x_full(:, free) = x_r;

w = x_full(:, 1:2:end);
xpos = linspace(0, L, nNodes);

figure; plot(tSol, w(:, floor(nNodes/2)));
grid on; xlabel('t'); ylabel('w(mid)');
title('Midspan response under VIV forcing (input-driven)');

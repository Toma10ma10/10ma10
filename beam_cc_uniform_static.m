%% beam_cc_uniform_static.m
clear; clc;

% ---- Parameters (SI) ----
E   = 210e9;     % Pa
I   = 1.0e-8;    % m^4
rho = 7850;      % kg/m^3 (unused in static)
A   = 1.0e-4;    % m^2   (unused in static)
L   = 1.0;       % m

nNodes = 100;            % "100点くらい"
q = -100.0;              % N/m  (downward)

% ---- Assemble ----
[M, K] = assembleBeamEB(E, I, rho, A, L, nNodes);
f = assembleUniformLoadEB(q, L, nNodes);

% ---- Clamped-Clamped BC ----
fixed = [1 2, 2*nNodes-1, 2*nNodes];   % w1,th1,wN,thN
[Kr, fr, free] = applyDirichletBC(K, f, fixed);

% ---- Static solve: K x = f ----
xr = Kr \ fr;

% expand to full dof vector
ndof = 2*nNodes;
x = zeros(ndof,1);
x(free) = xr;

% extract deflection and rotation
w = x(1:2:end);
theta = x(2:2:end);

fprintf("w(0)=%.3e, w(L)=%.3e, w(mid)=%.3e\n", w(1), w(end), w(floor(nNodes/2)));

xpos = linspace(0, L, nNodes).';
figure; plot(xpos, w, '-o'); grid on;
xlabel('x [m]'); ylabel('w [m]');
title('Clamped-Clamped Beam: Static deflection under uniform load');

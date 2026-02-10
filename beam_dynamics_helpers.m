function [Mr, Cr, Kr, free] = reduceByBC(M,C,K,fixed)
ndof = size(M,1);
mask = true(ndof,1);
mask(fixed) = false;
free = find(mask);

Mr = M(free, free);
Cr = C(free, free);
Kr = K(free, free);
end

function Ffull = expandNodalWForceToFullDOF(Fw)
[nT, nNodes] = size(Fw);
ndof = 2*nNodes;
Ffull = zeros(nT, ndof);
Ffull(:, 1:2:end) = Fw;
end

function dz = stateSpaceRHS(t, z, M, C, K, Ffun)
n = size(M,1);
x = z(1:n);
v = z(n+1:2*n);
f = Ffun(t);
a = M \ (f - C*v - K*x);
dz = [v; a];
end

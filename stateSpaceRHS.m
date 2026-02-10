function dz = stateSpaceRHS(t, z, M, C, K, Ffun)
n = size(M,1);
x = z(1:n);
v = z(n+1:2*n);
f = Ffun(t);
a = M \ (f - C*v - K*x);
dz = [v; a];
end

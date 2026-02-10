function [M, K] = assembleBeamEB(E, I, rho, A, L, nNodes)
nElem = nNodes - 1;
Le = L / nElem;
EI = E * I;
rhoA = rho * A;

ndof = 2*nNodes;
M = zeros(ndof, ndof);
K = zeros(ndof, ndof);

[me, ke] = beamElementEB(EI, rhoA, Le);

for e = 1:nElem
    n1 = e;
    n2 = e+1;
    dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2];

    M(dofs, dofs) = M(dofs, dofs) + me;
    K(dofs, dofs) = K(dofs, dofs) + ke;
end
end

function [me, ke] = beamElementEB(EI, rhoA, Le)
ke = (EI / Le^3) * ...
    [ 12,    6*Le,  -12,    6*Le;
      6*Le, 4*Le^2, -6*Le, 2*Le^2;
     -12,   -6*Le,   12,   -6*Le;
      6*Le, 2*Le^2, -6*Le, 4*Le^2];

me = (rhoA * Le / 420) * ...
    [156,     22*Le,   54,    -13*Le;
     22*Le,  4*Le^2,  13*Le,  -3*Le^2;
     54,      13*Le,  156,    -22*Le;
    -13*Le, -3*Le^2, -22*Le,   4*Le^2];
end

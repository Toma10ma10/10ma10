function f = assembleUniformLoadEB(q, L, nNodes)
nElem = nNodes - 1;
Le = L / nElem;

ndof = 2*nNodes;
f = zeros(ndof, 1);

fe = (q * Le / 12) * [6; Le; 6; -Le];

for e = 1:nElem
    n1 = e;
    n2 = e+1;
    dofs = [2*n1-1, 2*n1, 2*n2-1, 2*n2];
    f(dofs) = f(dofs) + fe;
end
end

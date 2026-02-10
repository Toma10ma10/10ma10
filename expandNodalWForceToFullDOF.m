function Ffull = expandNodalWForceToFullDOF(Fw)
[nT, nNodes] = size(Fw);
ndof = 2*nNodes;
Ffull = zeros(nT, ndof);
Ffull(:, 1:2:end) = Fw;
end

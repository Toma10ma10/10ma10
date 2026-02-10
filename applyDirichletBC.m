function [Kr, fr, free] = applyDirichletBC(K, f, fixed)
ndof = size(K,1);
mask = true(ndof,1);
mask(fixed) = false;
free = find(mask);

Kr = K(free, free);
fr = f(free);
end

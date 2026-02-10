function [Mr, Cr, Kr, free] = reduceByBC(M,C,K,fixed)
ndof = size(M,1);
mask = true(ndof,1);
mask(fixed) = false;
free = find(mask);

Mr = M(free, free);
Cr = C(free, free);
Kr = K(free, free);
end

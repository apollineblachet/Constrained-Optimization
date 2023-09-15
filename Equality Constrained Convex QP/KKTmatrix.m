function [K,d] = KKTmatrix(H,g,A,b)
[n,m] = size(A);
K = [H -A; -A' zeros(m,m)];
d = [-g; b];
end
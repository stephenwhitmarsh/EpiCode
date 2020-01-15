function d = dtwdist(Xi, Xj, varargin)

[m,n] = size(Xj);

d = zeros(m,1);

for j = 1 : m
    d(j) = dtw(Xi,Xj(j,:), varargin{:});
end

function h = subplottight(n,m,i)
[c,r] = ind2sub([m n], i);
ax = subplot('Position', [(c-1)/m + c*0.01 + 0.03, 1-(r)/n - r *0.01, 1/m*0.75, 1/n]);
if (nargout>0)
    h = ax;
end

function N = weighted_norm(a)
S = eye(size(a,1));
N = a'*S*a;
end

function Sigma = solve_ricatti(F,Q)

% Find Sigma that satisfies:
% Sigma = F * Sigma F' + Q;

r = size(F,1);
% Use VEC formula in Hamilton (1994)
Sigma_vec = (eye(r^2)-kron(F,F))\Q(:);

% reconstruct the matrix

Sigma = reshape(Sigma_vec,[r,r]);

% check it solves the equation
%check = Sigma - F*Sigma*F'-Q;
% it works!

end
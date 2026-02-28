
function J = QPO_Add_Harmonics(w,w_Q,max_int_2,max_int_1,J_QPO)
% QPO_Add_Harmonics  Identify new integer frequency combinations via CVX.
%
%   Given a set of frequencies w extracted from NAFF on the Jacobi
%   constant residual, solves an integer least-squares problem to find
%   which combinations [j1, j2] of the fundamental frequencies w_Q best
%   match them. New combinations not already in J_QPO are returned.
%
%   Uses CVX with the Gurobi solver for the mixed-integer L1 program:
%     minimize  ||J*w_Q - w||_1
%     subject to  ||J(:,1)||_inf <= max_int_1
%                 ||J(:,2)||_inf <= max_int_2
%
%   INPUTS:
%     w         - [M x 1] frequencies from NAFF [rad/TU]
%     w_Q       - [2 x 1] fundamental frequencies [w0; wP]
%     max_int_2 - maximum integer order for wP combinations
%     max_int_1 - maximum integer order for w0 combinations
%     J_QPO     - [K x 2] existing integer combination matrix
%
%   OUTPUT:
%     J         - [P x 2] new integer combinations to add (not in J_QPO)

cvx_clear
cvx_solver Gurobi
cvx_precision high
cvx_begin quiet
    variable J(length(w), length(w_Q)) integer
    minimize( norm((J*w_Q - w.'), 1) )
    subject to
        norm(J(:,2), inf) <= max_int_2
        norm(J(:,1), inf) <= max_int_1
cvx_end

J = round(full(J));

% Keep only entries that closely match a target frequency (within 5e-2 rad/TU)
index1 = abs((J*w_Q - w.')) <= 5e-2;
% Remove combinations already present in the basis
index2 = ismember(J, J_QPO, 'rows');
% Flip sign of any negative-frequency combinations
index3 = J*w_Q < 0;
J(index3,:) = -J(index3,:);

J = J(index1 & ~index2, :);
J = unique(J, 'rows');
end


% SOLVEMS

function [V, J, cpu_time] = solveMS(A, B, Q, R, N, nx, nu, xbar0, xbarf)

    % Define optimization variables
    V = sdpvar(nu * N, 1);
    
    % Useful matrices
    EN = [zeros(nx, N * nx) eye(nx)];
    
    % System Dynamics
    Xbar = A * xbar0 + B * V;

    % Define constraints
    Constraints = [];
    
    % Final mean state constraint
    Constraints = [Constraints, EN * Xbar == xbarf];      
    
    % Objective Function
    Objective = Xbar' * Q * Xbar + V' * R * V;

    % Solve the Problem
    options = sdpsettings('solver','mosek');
    sol = optimize(Constraints,Objective,options);
    
    % Get performance measure values
    cpu_time = sol.solvertime;
    J = value(Objective);

    % Return optimal control
    V = value(V);
end

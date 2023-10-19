% SOLVECS

function [K, F, J, cpu_time, problem] = solveCS(B, Q, R, N, nx, nu, PZ, ...
    PhatPrior0, PhatPriorf)

    % Define optimization variables
    F = [];
    for i = 1:N            
        % F = [F, zeros(nu * (i - 1), nx); sdpvar(nu, nx * i)];  % Full lower-triangular,  optimal but slow for large N            
        % F = blkdiag(F, sdpvar(nu, nx));                        % F restricted to block-diagonal (suboptimal, but way faster)
        F = blkdiag(F, [diag(sdpvar(nu, 1)), diag(sdpvar(nu, 1))]);   % nx = 2*nu, for double integrator each dimension maybe treated independently  
    end
    F = [F,zeros(nu * N, nx)];
    
    % Useful matrices
    IplusBF = eye((N + 1) * nx) + B * F;
    EN = [zeros(nx, N * nx) eye(nx)];
    PZhalf = computeS(PZ)';
    PtargetHalfInv = inv(computeS(PhatPriorf));

    % Define constraints
    Constraints = [];  
    
    % Maximum final covariance constraint
    Constraints = [Constraints, ...
        norm(PZhalf * IplusBF' * EN' * PtargetHalfInv) <= 1];
    
    % Objective Function
%     tic
%     Objective = trace((IplusBF' * Q * IplusBF + F' * R * F) * PZ);
%     toc
    Objective = trace(Q * PZ) + 2 * trace(PZ * Q * B * F ) + trace(F * PZ * F' * (B' * Q * B + R));
    
    % Solve the Problem
    options = sdpsettings('solver','mosek');
    sol = optimize(Constraints,Objective,options);
    problem = sol.problem;
    
    % Get performance measure values
    cpu_time = sol.solvertime;
    J = value(Objective);

    % Return optimal control
    F = value(F);
    K = F/(eye((N + 1) * nx) + B * F); 
end

function PZhalf = computeS(SigmaY)
    [L,D] = ldl(SigmaY);
    PZhalf = L * sqrtm(D);
end

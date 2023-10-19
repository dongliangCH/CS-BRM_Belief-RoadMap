% SOLVECS

function [K, J, cpu_time, problem] = solveCS(B, Q, R, N, nx, nu, PZ, ...
    PhatPrior0, PhatPriorf, useHistory)

    % Define optimization variables
    K = [];
    for i = 1:N
        if useHistory
            % With history: optimal, but slow for large N
            K = [K, zeros(nu * (i - 1), nx); sdpvar(nu, nx * i)];
        else
            % Without feedback history (suboptimal, but way faster)
            K = blkdiag(K, sdpvar(nu, nx));
        end
    end
    K = [K,zeros(nu * N, nx)];
    
    % Useful matrices
    IplusBK = eye((N + 1) * nx) + B * K;
    EN = [zeros(nx, N * nx) eye(nx)];
    PZhalf = computeS(PZ)';
    PtargetHalfInv = inv(computeS(PhatPriorf));

    % Define constraints
    Constraints = [];  
    
%     % Intermediate covariance constraint
%     k = ceil( N/2 );
%     Ek = [zeros(nx, k * nx) eye(nx) zeros(nx,(N - k) * nx)];
%     PhatPriorinterm = 3 * (PhatPrior0 + PhatPriorf) / 2;
%     PintermediateHalfInv = inv(computeS(PhatPriorinterm ));
%     Constraints = [Constraints, ...
%         norm(PZhalf * IplusBK' * Ek' * PintermediateHalfInv) <= 1];
    
    % Maximum final covariance constraint
    Constraints = [Constraints, ...
        norm(PZhalf * IplusBK' * EN' * PtargetHalfInv) <= 1];
    
    % Objective Function
    Objective = trace((IplusBK' * Q * IplusBK + 1 * K' * R * K) * PZ);

    % Solve the Problem
    options = sdpsettings('solver','mosek');
    sol = optimize(Constraints,Objective,options);
    problem = sol.problem;
    
    if problem
        temp = 1;
    end
    
    % Get performance measure values
    cpu_time = sol.solvertime;
    J = value(Objective);

    % Return optimal control
    K = value(K);
end

function PZhalf = computeS(SigmaY)
    [L,D] = ldl(SigmaY);
    PZhalf = L * sqrtm(D);
end

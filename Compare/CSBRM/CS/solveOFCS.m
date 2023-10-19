% SOLVEOFCS

function [K, V, J, cpu_time] = solveOFCS(A, B, Q, R, N, nx, nu, PZ, ...
    xbar0, xbarf, Pf, PPt, useHistory)

    % Define optimization variables
    V = sdpvar(nu * N, 1);
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
    
    % Filter error covariance
    PtildeN = PPt(:, :, end);
    PPtHalf = zeros(size(PPt));
    for k = 1:(N + 1)
        PPtHalf(:, :, k) = computeS(PPt(:, :, k))';
    end
    
    % Useful matrices
    IplusBK = eye((N + 1) * nx) + B * K;
    EN = [zeros(nx, N * nx) eye(nx)];
    PZhalf = computeS(PZ)';
    PtargetHalfInv = inv(computeS(Pf - PtildeN));

%     SigmaY = ScriptA*PS.Sigma0*ScriptA'+ScriptD*ScriptD';
%     SY = computeS(SigmaY);
    
    % System Dynamics
    Xbar = A * xbar0 + B * V;
%     Z = EN * IplusBK * PZhalf;

    % Define constraints
    Constraints = [];
    
    % Final mean state constraint
    Constraints = [Constraints, EN * Xbar == xbarf];      
    
    % Maximum final covariance constraint
%     Constraints = [Constraints, [(PS.Sigmaf - PS.PtildeN) Z; Z' eye((N+1)*nx)] >= 0];
    Constraints = [Constraints, ...
        norm(PZhalf * IplusBK' * EN' * PtargetHalfInv) <= 1];
    
    % Objective Function
    Objective = Xbar' * Q * Xbar + V' * R * V ...
        + trace((IplusBK' * Q * IplusBK + K' * R * K) * PZ);

    % Solve the Problem
    options = sdpsettings('solver','mosek');
    sol = optimize(Constraints,Objective,options);
    
    % Get performance measure values
    cpu_time = sol.solvertime;
    J = value(Objective);

    % Return optimal control
    V = value(V);
    K = value(K);
end

function PZhalf = computeS(SigmaY)
    [L,D] = ldl(SigmaY);
    PZhalf = L * sqrtm(D);
end

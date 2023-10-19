function [CovCost, CollisonCost, K, problem] = PRMCovC(PhatPrior0, PtildePrior0, PhatPriorf, PtildePriorf, N, param, Xbar, V, world)

% Notation: Ak is single step of matrix, AA is 3d array with each Ak
% stacked in the 3rd dimension, A is block matrix

%% Problem setup

% Number of steps N
% Choose time step dt = 0.2

% Dynamics
dt = param.dt;
nx = 4;
nu = 2;
ny = 4;
nw = 4;
Ak = [1 0 dt 0; 0 1 0 dt; 0 0 1 0; 0 0 0 1];
Bk = [dt ^ 2 / 2 0; 0 dt ^ 2 / 2; dt 0; 0 dt];
Gk = sqrt(dt) * diag([0.01 0.01 0.01 0.01]);

% Sequences of system matricies
AA = repmat(Ak, [1, 1, N]);
BB = repmat(Bk, [1, 1, N]);
GG = repmat(Gk, [1, 1, N]);

% Make block matricies for dynamics
A = makeStateMatrix(AA, N);
B = makeInputMatrix(AA, BB, N);

% Cost
Qk = blkdiag(2, 2, 4, 4);
Rk = blkdiag(1, 0.1);

% Make block matricies for cost
Q = [];
R = [];
for i = 1:N
    Q = blkdiag(Q, Qk);
    R = blkdiag(R, Rk);
end
Q = blkdiag(Q, zeros(nx));


%% Kalman filter
% Observation model

Ck = zeros(4,4);
CC = repmat(Ck, [1, 1, N + 1]);
Dk = 0.01 * eye(ny);
DD = repmat(Dk, [1, 1, N + 1]);

for i = 1 : N + 1
    
    if 3 < Xbar(1 + 4*(i-1)) && Xbar(1 + 4*(i-1)) < 5 && 0 < Xbar(2 + 4*(i-1)) && Xbar(2 + 4*(i-1)) < 2
        Ck = eye(4,4);
    else
        if 3 < Xbar(1 + 4*(i-1)) && Xbar(1 + 4*(i-1)) < 5 && 3.5 < Xbar(2 + 4*(i-1)) && Xbar(2 + 4*(i-1)) < 5
            Ck = eye(4,4);
        else
            Ck = eye(4); 
        end
    end
    CC(:, :, i) = Ck;
end


% Initial state covariance
P0 = PhatPrior0 + PtildePrior0;

% Initial estimation error covariance  PtildePrior0
% Initial estimated state covariance  PhatPrior0 

% Target maximum final covariance
Pf = PhatPriorf + PtildePriorf;

% Solve for sequence of kalmain gains LL, error covariances PPtilde, and prior
% error covariances PPtm
[LL, PPtilde, PPtm] = solveKF(AA, GG, CC, DD, N, PtildePrior0);

CollisonCost = 10000;
CovCost = 10000;
problem = 1;
K = 0;

if ~isMatrixPD(PtildePriorf  * 1.2 + 0.01 * eye(4) - PPtm(:, :, end))
    
    
L = makeInputMatrixType2(AA, LL, N);

%% Covariance steering
% Innovation process covariance
PYtilde = zeros((N + 1) * ny, (N + 1) * ny);
for i = 0:N
    ndx = (i * ny) + (1:ny);
    PYtilde(ndx, ndx) = CC(:, :, i + 1) * PPtm(:, :, i + 1) ...
        * CC(:, :, i + 1)' + DD(:, :, i + 1) * DD(:, :, i + 1)';
end

% Feedback process covariance
PZ = A * PhatPrior0 * A' + L * PYtilde * L';

% Flag if should use history of disturbances in feedback
useHistory = false;

% Solve covariance steering with yalmip
[K, CovCost, cpu_time, problem] = solveCS(B, Q, R, N, nx, nu, PZ, PhatPrior0, PhatPriorf + PPtm(:, :, end) - PPtilde(:, :, end), useHistory);
CovCost = dt * CovCost;

if ~problem

    % Filtered state covariance (block matrix)
    IplusBK = eye((N + 1) * nx) + B * K;
    Phat = IplusBK * PZ * IplusBK';

    % Mean and covariance as normal sequences
    PPhat = zeros(nx, nx, N + 1);
    PP = zeros(nx, nx, N + 1);
    xbar = zeros(nx, N + 1);
    for k = 0:N
        Ek = [zeros(nx, k * nx) eye(nx) zeros(nx,(N - k) * nx)];
        PPhat(:, :, k + 1) = Ek * Phat * Ek';
        PP(:, :, k + 1) = PPhat(:, :, k + 1) + PPtilde(:, :, k + 1);
        xbar(:, k + 1) = Ek * Xbar;
    end

    %% Monte Carlo
    % Number of monte carlo trials to display (can be zero)
    MCnum = 50;

    figure(10); 
    hold on
    CollisonCost = montecarlo_brm(MCnum, AA, BB, GG, CC, DD, LL, xbar(:,1), ...
        N, PhatPrior0, PtildePrior0, V, K, useHistory, world);

    %% Plot result

    % Step numbers to draw covariance ellipses at
    ellipse_steps = 0:N:N;

    figure(1);
    drawResult_brm(MCnum, AA, BB, GG, CC, DD, LL, xbar(:,1), xbar(:,end), xbar, ...
        Pf, N, PhatPrior0, PtildePrior0, V, K, PP, PPhat, ...
        PPtilde, ellipse_steps, useHistory)

end

end

end


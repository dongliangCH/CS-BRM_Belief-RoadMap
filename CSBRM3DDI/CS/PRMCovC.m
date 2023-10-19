function [CovCost, CollisonCost, K, problem] = PRMCovC(PhatPrior0, PtildePrior0, PhatPriorf, PtildePriorf, N, param, Xbar, V, world)

% Notation: Ak is single step of matrix, AA is 3d array with each Ak
% stacked in the 3rd dimension, A is block matrix

%% Problem setup

% Number of steps N
% Choose time step dt = 0.2

% Dynamics
dt = param.dt;
nx = 6;
nu = 3;
ny = 6;
nw = 6;
Ak = [1 0 0 dt 0 0; 0 1 0 0 dt 0; 0 0 1 0 0 dt; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
Bk = [dt ^ 2 / 2 0 0; 0 dt ^ 2 / 2 0; 0  0 dt ^ 2 / 2; dt 0 0 ; 0 dt 0; 0 0 dt];
Gk = sqrt(dt) * diag([0.001 0.001 0.001 0.04 0.04 0.04]);

% Sequences of system matricies
AA = repmat(Ak, [1, 1, N]);
BB = repmat(Bk, [1, 1, N]);
GG = repmat(Gk, [1, 1, N]);

% Make block matricies for dynamics
A = makeStateMatrix(AA, N);
B = makeInputMatrix(AA, BB, N);

% Cost
Qk = blkdiag(5, 5, 5, 1, 1, 1);
Rk = blkdiag(0.1, 0.1, 0.1);

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
Ck = eye(6);
CC = repmat(Ck, [1, 1, N + 1]);
% Dk = 0.06 * eye(ny);
% DD = repmat(Dk, [1, 1, N + 1]);

for i = 1 : N + 1    
    obs = [0.75, 2; 1.25, 2; 1.25, 0.75; 2.25, 1.25; 2.75, 0.75; 3, 2];
    dist = [Xbar(1 + 6*(i-1)), Xbar(2 + 6*(i-1))] - obs;
    min_dist = min(sqrt(sum(dist.^2, 2)));
    if min_dist < 0.6
       Dk = 0.06 * eye(ny);
    else
       Dk = 0.03 * eye(ny);
    end
    DD(:, :, i) = Dk;
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

if ~isMatrixPD(PtildePriorf - PPtm(:, :, end))
    return;
end
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



% Solve covariance steering with yalmip
[K, F, CovCost, cpu_time, problem] = solveCS(B, Q, R, N, nx, nu, PZ, PhatPrior0, PhatPriorf + PPtm(:, :, end) - PPtilde(:, :, end));

CovCost = dt * CovCost;

if ~problem

    % Filtered state covariance (block matrix)
    IplusBF = eye((N + 1) * nx) + B * F;
    Phat = IplusBF * PZ * IplusBF';

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
    MCnum = 100;

    figure(10); 
    hold on
    useHistory = true;
    CollisonCost = montecarlo_brm(MCnum, AA, BB, GG, CC, DD, LL, xbar, ...
        N, PhatPrior0, PtildePrior0, V, K, useHistory, world);

end

end


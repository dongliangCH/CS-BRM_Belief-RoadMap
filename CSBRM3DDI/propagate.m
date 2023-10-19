function [endP0, endPtilde, CovCost, CollisonProb, K] = propagate(P0, PtildePrior0, N, param, Xbar, world)

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
Gk = sqrt(dt) * diag([0.002 0.002 0.002 0.05 0.05 0.05]);

% Sequences of system matricies
AA = repmat(Ak, [1, 1, N]);
BB = repmat(Bk, [1, 1, N]);
GG = repmat(Gk, [1, 1, N]);

% Make block matricies for dynamics
A = makeStateMatrix(AA, N);
B = makeInputMatrix(AA, BB, N);

% Cost
Qk = 2 * blkdiag(2, 2, 2, 1, 1, 1);
Rk = 0.5 * blkdiag(1, 1, 1);

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
Dk = 0.08 * eye(ny);
DD = repmat(Dk, [1, 1, N + 1]);

% Initial state covariance
PhatPrior0 = P0 - PtildePrior0;

% Initial estimation error covariance  PtildePrior0
% Initial estimated state covariance  PhatPrior0 

% Solve for sequence of kalmain gains LL, error covariances PPtilde, and prior
% error covariances PPtm
[LL, PPtilde, PPtm] = solveKF(AA, GG, CC, DD, N, PtildePrior0);

%% Time-varying LQR tracking

S(:,:,N+1) = Qk;
for k = N:-1:1
    
    K(:,:,k) = (BB(:,:,k)' * S(:,:,k+1) * BB(:,:,k) + Rk) \ BB(:,:,k)' * S(:,:,k+1) * AA(:,:,k);
    S(:,:,k) = Qk + AA(:,:,k)' * S(:,:,k+1) * AA(:,:,k) - AA(:,:,k)' * S(:,:,k+1) * BB(:,:,k) * K(:,:,k);

end

%% State covariance and estimated state covariance

PPhat(:,:,1) = P0 - PPtilde(:,:,1);
PP(:,:,1) = P0(1:3,1:3);
Jc = 0;
for k = 1:N
    PPhat(:,:,k+1) = (Ak - Bk * K(:,:,k)) * PPhat(:,:,k) * (Ak - Bk * K(:,:,k))' + LL(:,:,k+1) * Ck * PPtm(:,:,k+1);
    PP(:,:,k+1) = PPhat(1:3,1:3,k+1) + PPtilde(1:3,1:3,k+1);
    PP_all = PPhat(:,:,k+1) + PPtilde(:,:,k+1);
    Jc = Jc + trace( PP_all * (Qk + K(:,:,k)' * Rk * K(:,:,k)) );
end

endP0 = PPhat (:, :, end) + PPtilde(:, :,end);
endPtilde = PPtilde(:, :,end);

CovCost = dt * Jc;

xbar = zeros(nx, N + 1);
for k = 0:N
    Ek = [zeros(nx, k * nx) eye(nx) zeros(nx,(N - k) * nx)];
    xbar(:, k + 1) = Ek * Xbar;
end

%% Monte Carlo
% Number of monte carlo trials to display (can be zero)
MCnum = 100;

% figure(10); hold on
CollisonProb = montecarlo(MCnum, AA, BB, GG, CC, DD, LL, xbar, ...
    N, PhatPrior0, PtildePrior0, K, world);

%% Plot result

% % Step numbers to draw covariance ellipses at
% ellipse_steps = 0:1:N;
% 
% figure(1);
% drawResult(xbar, PP, PPtilde, ellipse_steps)

end


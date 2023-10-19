function  [x_MC_end, xhatPrior_MC_end] = MC_plan(xbar0, x0_MC, xhatPrior0_MC, Xbar, PtildePrior0, V, K, useHistory, param, world)

% Number of steps N
N = size(K,1)/2; 

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

%% Kalman filter
% Observation model

% Ck = zeros(6,6);
% CC = repmat(Ck, [1, 1, N + 1]);
Dk = 0.01 * eye(ny);
DD = repmat(Dk, [1, 1, N + 1]);

for i = 1 : N + 1
    
%     if 3 < Xbar(1 + 4*(i-1)) && Xbar(1 + 4*(i-1)) < 5 && 0 < Xbar(2 + 4*(i-1)) && Xbar(2 + 4*(i-1)) < 2
%         Ck = eye(4,4);
%     else
%         if 3 < Xbar(1 + 4*(i-1)) && Xbar(1 + 4*(i-1)) < 5 && 3.5 < Xbar(2 + 4*(i-1)) && Xbar(2 + 4*(i-1)) < 5
%             Ck = eye(4,4);
%         else
            Ck = eye(4); 
%         end
%     end
    CC(:, :, i) = Ck;
end

% Solve for sequence of kalmain gains LL, error covariances PPtilde, and prior
% error covariances PPtm
[LL, ~, ~] = solveKF(AA, GG, CC, DD, N, PtildePrior0);

%%
% Sizes of inputs
nv = ny;

hold on
% rng(0);

    zPrior0 = xhatPrior0_MC - xbar0;

    U = zeros(nu,N);
    x_MC = zeros(nx,N+1);
    xhat_MC = zeros(nx,N+1);
    xhatPrior_MC = zeros(nx,N+1);
    z_MC = zeros(nx,N+1);
    y_MC = zeros(ny,N+1);
    x_MC(:,1) = x0_MC;

    for k = 1:(N + 1)

        v = randn(nv, 1);
        w = randn(nw, 1);

        % Measurement
        y_MC(:, k) = CC(:, :, k) * x_MC(:, k) + DD(:, :, k) * v;

        % Prior estimate
        if k == 1
            xhatPrior_MC(:,k) = xhatPrior0_MC;
        else
            xhatPrior_MC(:, k) = AA(:, :, k - 1) * xhat_MC(:, k - 1) ...
                + BB(:, :, k - 1) * U(:, k - 1);
        end

        % Innovation process
        ytilde_MC(:, k) = y_MC(:, k) - CC(:, :, k) * xhatPrior_MC(:, k);

        % Feedback process
        if k == 1
            z_MC(:,k) = zPrior0 + LL(:,:,k) * ytilde_MC(:,k);
        else
            z_MC(:,k) = AA(:, :, k - 1) * z_MC(:,k-1) ...
                + LL(:,:,k) * ytilde_MC(:,k);
        end

        % Control
        if k <= N
            if useHistory
                U(:, k) = V((k - 1) * nu + 1:k * nu) ...
                    + K((k - 1) * nu + 1:k * nu, 1:k * nx) * vec(z_MC(:, 1:k));
            else
                U(:, k) = V((k - 1) * nu + 1:k * nu) ...
                    + K((k - 1) * nu + 1:k * nu, ...
                    (k - 1) * nx + 1:k * nx) * z_MC(:,k);
            end
        end

        % Dynamics
        if k <= N
            x_MC(:, k + 1) = AA(:, :, k) * x_MC(:, k) ...
                + BB(:, :, k) * U(:, k) + GG(:, :, k) * w;
        end

        % Filtered process
        if k == 1
            xhat_MC(:,k) = xhatPrior0_MC + LL(:,:,k) * ytilde_MC(:,k);
        else
            xhat_MC(:,k) = AA(:, :, k - 1) * xhat_MC(:,k-1) ...
                + BB(:, :, k - 1)*U(:,k-1)...
                + LL(:,:,k) * ytilde_MC(:,k);
        end

    end
    
    x_MC_end = x_MC(:, end);
    xhatPrior_MC_end = xhatPrior_MC(:,end);
    hold on;
    plot(x_MC(1,:), x_MC(2,:), 'color', 0.5 * ones(3, 1));

end

function b = vec(a)

b = a(:);

end
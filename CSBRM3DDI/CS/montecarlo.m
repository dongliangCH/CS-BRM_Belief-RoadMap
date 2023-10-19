

function collisionProb = montecarlo(MCnum, AA, BB, GG, CC, DD, LL, xbar, ...
    N, PhatPrior0, PtildePrior0, K, world)


% Sizes of inputs
nx = size(AA, 1);
nu = size(BB, 2);
ny = size(CC, 1);
nw = size(GG, 2);
nv = ny;

% rng(0);
collision = 0;

for mc = 1:MCnum

    xhatPrior0_MC = mvnrnd(xbar(:,1) - xbar(:,1), PhatPrior0 + 0.00001 * eye(6), 1)';
    xtildePrior0 = mvnrnd(zeros(nx,1), PtildePrior0, 1)';
    x0_MC = xhatPrior0_MC + xtildePrior0;

    U = zeros(nu,N);
    x_MC = zeros(nx,N+1);
    xhat_MC = zeros(nx,N+1);
    xhatPrior_MC = zeros(nx,N+1);
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

        xhat_MC(:,k) =  xhatPrior_MC(:,k) + LL(:,:,k) * (y_MC(:,k) - CC(:,:,k) * xhatPrior_MC(:,k));
        
        % Control
        if k <= N
            U(:, k) = - K(:,:,k) * xhat_MC(:,k);
        end

        % Dynamics
        if k <= N
            x_MC(:, k + 1) = AA(:, :, k) * x_MC(:, k) ...
                + BB(:, :, k) * U(:, k) + GG(:, :, k) * w;
        end

    end
    
% %     hold on;
% %     plot3( x_MC(1,:) + xbar(1,:), x_MC(2,:) + xbar(2,:), x_MC(3,:) + xbar(3,:), 'color', 0.6 * ones(3, 1));
    
    if MeanCollisionCheck(x_MC(1:3,:) + xbar(1:3,:), world, 3)
        collision = collision + 1;
    end
    
end

collisionProb = collision / MCnum;

end


function collisionCost = montecarlo(MCnum, AA, BB, GG, CC, DD, LL, xbar0, ...
    N, PhatPrior0, PtildePrior0, V, K, useHistory)

global Cspace

% Sizes of inputs
nx = size(AA, 1);
nu = size(BB, 2);
ny = size(CC, 1);
nw = size(GG, 2);
nv = ny;

hold on
rng(0);
collision = 0;

for mc = 1:MCnum

    xhatPrior0_MC = mvnrnd(xbar0, PhatPrior0, 1)';
    xtildePrior0 = mvnrnd(zeros(nx,1), PtildePrior0, 1)';
    x0_MC = xhatPrior0_MC + xtildePrior0;

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
    hold on;
    plot(x_MC(2,:).*20, x_MC(1,:).*20, 'color', 0.6 * ones(3, 1));
    
    if ~MeanCollisionCheck(20 * [x_MC(1,:); x_MC(2,:)], Cspace)
        collision = collision + 1;
    end
    
end

collisionP = collision / MCnum;
collisionCost = 1000 * collisionP;

x1 = xlabel('$x_1$', 'interpreter', 'latex');
y1 = ylabel('$x_2$', 'interpreter', 'latex');
set(x1,'FontSize',18);
set(y1,'FontSize',18);
set(gca,'FontSize',16,'FontName','Times');
xlim([-10 510])
ylim([-10 510])    
end

function b = vec(a)

b = a(:);

end
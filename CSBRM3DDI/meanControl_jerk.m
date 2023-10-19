function [xbar, V, Mcost, N, Xbar] = meanControl_jerk(startp, endp, param)

% Initial mean
xbar0 = startp;

% Target mean
xbarf = endp;

velavg = param.velavg;  % Average velocity
dim = 3;
tf = norm (xbar0(1:dim) - xbarf(1:dim))/velavg; % Final time

% Number of steps
dt = param.dt;
N = ceil (tf / dt);

xbar0 = [xbar0; 0; 0; 0];
xbarf = [xbarf; 0; 0; 0];

% Dynamics
nx = 9;
nu = 3;
Ak = [ 1 0 0 dt 0  0  dt^2/2 0      0;
       0 1 0 0  dt 0  0      dt^2/2 0;
       0 0 1 0  0  dt 0      0      dt^2/2;
       0 0 0 1  0  0  dt     0      0;
       0 0 0 0  1  0  0      dt     0;
       0 0 0 0  0  1  0      0      dt;
       0 0 0 0  0  0  1      0      0;
       0 0 0 0  0  0  0      1      0;
       0 0 0 0  0  0  0      0      1];
Bk = [dt^3/4 0      0;
      0      dt^3/4 0; 
      0      0      dt^3/4;
      dt^2/2 0      0;
      0      dt^2/2 0; 
      0      0      dt^2/2; 
      dt     0      0; 
      0      dt     0; 
      0      0      dt];
% Ak = [1 0 0 dt 0 0; 0 1 0 0 dt 0; 0 0 1 0 0 dt; 0 0 0 1 0 0; 0 0 0 0 1 0; 0 0 0 0 0 1];
% Bk = [dt ^ 2 / 2 0 0; 0 dt ^ 2 / 2 0; 0  0 dt ^ 2 / 2; dt 0 0 ; 0 dt 0; 0 0 dt];

% Cost
Qk = 0 * blkdiag(0, 0, 0, 0, 0, 0, 2, 2, 2);
Rk = 1 * blkdiag(2, 2, 2);

% Sequences of system matricies
AA = repmat(Ak, [1, 1, N]);
BB = repmat(Bk, [1, 1, N]);

% Make block matricies for cost
Q = [];
R = [];
for i = 1:N
    Q = blkdiag(Q, Qk);
    R = blkdiag(R, Rk);
end
Q = blkdiag(Q, zeros(nx));

% Make block matricies for dynamics
A = makeStateMatrix(AA, N);
B = makeInputMatrix(AA, BB, N);

t = linspace(0, tf, N+1);
Vr = xbar0 + (xbarf - xbar0) * t / tf;
Xr = [];
for i=1 : N + 1
    Xr((i-1) * nx + 1 : i * nx, 1) = Vr(:,i);
end

% Mean control
M = B' * Q * B + R;
Minv = inv(M);
AN = A((nx * N + 1) : nx * (N + 1), :);
BN = B((nx * N + 1) : nx * (N + 1), :);
TempM = B' * Q * (A * xbar0 - Xr);
V = Minv * ( - TempM + BN' * inv(BN * Minv * BN') * (xbarf - AN * xbar0 + BN * Minv * TempM));

% Mean trajectory
Xbar = A * xbar0 + B * V;

% Cost for mean control
Mcost = dt * (Xbar' * Q * Xbar + V' * R * V);

for k = 0:N
    Ek = [zeros(nx, k * nx) eye(nx) zeros(nx,(N - k) * nx)];
    xbar(:, k + 1) = Ek * Xbar;
end

V = [];
Xbar = [];
for i = 1:size(xbar,2)
    Xbar = [Xbar; xbar(1:6,i)];
end
for i = 2:size(xbar,2)
    V = [V; xbar(7:9,i)];
end

xbar = xbar(1:6, :);

end
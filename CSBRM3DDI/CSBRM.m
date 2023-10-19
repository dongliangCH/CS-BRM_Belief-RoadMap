clear
addpath(genpath('C:\Program Files\Mosek\10.1'))
addpath(genpath('C:\Users\Dongliang\Documents\MATLAB\YALMIP-master'))
addpath('CS');
rng(1000);
dim = 3; 
state_dim = 6;

% Initial state distribution (mean/covariance)
start_cord = [0.75, 1.25, 1.5, 0, 0, 0];
P0 = 0.4 * diag([0.01, 0.01, 0.01, 0.01, 0.01, 0.01]);

% Initial estimation error covariance
PtildePrior0 = 0.4 * P0;

% Initial estimated state covariance
PhatPrior0 = P0 - PtildePrior0;

landingPoints = [2, 0.5, 0.3, 0, 0, 0; 
                 2, 2.5, 0.3, 0, 0, 0;
                 3.5, 1.75, 0.3, 0, 0, 0];

% Create random world
world = createKnownWorld([4,3,2.2],[0,0,0.2],dim);

param.dt = 0.05;
param.velavg = 1;

% Maximum steplength
segmentLength = 1;
% Neighbor distance
neb = 1;

samples = 30;
points = start_cord;
Covs(1,1) = {P0};
Covs(1,2) = {PtildePrior0};

figure(1); hold on
plot3(start_cord(1), start_cord(2), start_cord(3), 'Marker','s','MarkerSize',8,'MarkerEdgeColor','[0.8500 0.3250 0.0980]','MarkerFaceColor','[0.8500 0.3250 0.0980]');
for i = 1:size(landingPoints,1)
    plot3(landingPoints(i,1),landingPoints(i,2),landingPoints(i,3),'Marker','.','MarkerSize',12,'MarkerEdgeColor','[0.8500 0.3250 0.0980]');
end
plotCovariance(start_cord, P0);
plotWorld(world, dim); 
 
%% Sample nodes
i = 1;
while i < samples        
    randomPoint = zeros(1, state_dim);
    min_dist = 0;
    while min_dist <= 0.5            
        for j = 1:dim
            randomPoint(1, j) = (world.endcorner(j) - world.origincorner(j) - 0.7) * rand + world.origincorner(j) + 0.4;
        end
        % find the vertix that is closest to randomPoint (Eucl. dist. between positions)
        tmp = points(:, 1 : dim) - randomPoint(1 : dim);
        sqr_dist = sum(tmp.^2,2);
        [min_dist, idx] = min(sqr_dist);             
    end 
    if min_dist > segmentLength^2
        % generate a new point that is closest to randomPoint, segmentLength away from points(idx,1:dim)
        Vect = randomPoint(1:dim)-points(idx,1:dim);
        Vect = Vect/norm(Vect);
        new_point(1 : dim) = points(idx, 1 : dim) + Vect * segmentLength;
    else
        new_point(1 : dim) = randomPoint(1 : dim);
    end
    new_point(dim + 1 : state_dim) = 0.5 * [-1 + (1 - (-1)) * rand, -1 + (1 - (-1)) * rand, -1 + (1 - (-1)) * rand];
    
    % check if the new_point is in collision
    if collision_point(new_point, world) == 0
%         [meanTraj, V, MCost, N, Xbar]  = meanControl(points(idx, :)', new_point', param);  
        [meanTraj, V, ~, N, Xbar]  = meanControl_jerk(points(idx, :)', new_point', param); 
        if ~MeanCollisionCheck(meanTraj(1:dim, :), world, dim)
            [endP0, endPtildePrior0, CovCost, CollisionProb, K] = propagate( Covs{idx, 1}, Covs{idx, 2}, N, param, Xbar, world ); 
            if CollisionProb < 0.1
                plot3(meanTraj(1,:), meanTraj(2,:), meanTraj(3,:), 'color', 'g', 'LineWidth', 1);
                i = i + 1;                        
                points = [points; new_point];
                idx_newpoint = size(points,1);
                Covs(idx_newpoint,1) = {endP0};
                Covs(idx_newpoint,2) = {endPtildePrior0};
                plotCovariance(new_point, Covs{idx_newpoint,1});
            end                           
        end
    end
end

for i = 1:size(landingPoints,1)             
    new_point = landingPoints(i, :);
    % find the vertix that is closest to randomPoint (Eucl. dist. between positions)
    tmp = points(:, 1 : dim) - new_point(1 : dim);
    sqr_dist = sum(tmp.^2,2);
    [min_dist, idx] = min(sqr_dist);          
    % check if the new_point is in collision
    if collision_point(new_point, world) == 0
%         [meanTraj, V, MCost, N, Xbar]  = meanControl(points(idx, :)', new_point', param);  
        [meanTraj, V, ~, N, Xbar]  = meanControl_jerk(points(idx, :)', new_point', param); 
        if ~MeanCollisionCheck(meanTraj(1:dim, :), world, dim)
            [endP0, endPtildePrior0, CovCost, CollisionProb, K] = propagate( Covs{idx, 1}, Covs{idx, 2}, N, param, Xbar, world ); 
            plot3(meanTraj(1,:), meanTraj(2,:), meanTraj(3,:), 'color', 'g', 'LineWidth', 1);                       
            points = [points; new_point];
            idx_newpoint = size(points,1);
            Covs(idx_newpoint,1) = {endP0};
            Covs(idx_newpoint,2) = {endPtildePrior0};
            plotCovariance(new_point, Covs{idx_newpoint,1});                        
        end
    end
end
[points, Covs, DistM] = SampleMultipleVel(points, Covs, dim, size(landingPoints,1));

%% Graph construction 
% The ith row of ChildM contains the Children of node i 
ChildM = [];
% The EdgesCost(i,j) is the cost to go from node i to node j
EdgesCost = [];

% Adding edges
for i = 1:size(points, 1)   
    n = 0;    
    vec = points(i,1:dim) - points(:,1:dim);  % Euclidean distances
    dist = sqrt(sum(vec.^2, 2));     
    for k = 1:length(dist)                    % not connect nodes that are too close
        if dist(k) < 0.25
            dist(k) = 20;
        end
    end    
    near_idx = find(dist <= DistM(i));
    for j = 1:size(near_idx, 1)            
        vel = points(near_idx(j),1:dim) - points(i,1:dim);
        vel = vel / norm(vel);
        if norm(points(near_idx(j),4:6))<0.00001 || norm(points(i,4:6))<0.00001 || sum(abs(vel - points(i,4:6)/norm(points(i,4:6)))) < 0.00001 
            if norm(points(near_idx(j),4:6))>0.00002 && sum(abs(vel + points(near_idx(j),4:6)/norm(points(near_idx(j),4:6)))) < 0.00001 
                continue
            end
            [~, ~, MCost, N, ~]  = meanControl( points(i,1:6)', points(near_idx(j),1:6)', param);
            [meanTraj, V, ~, ~, Xbar]  = meanControl_jerk( points(i,1:6)', points(near_idx(j),1:6)', param);
            if ~MeanCollisionCheck(meanTraj(1:3, :), world, dim)
                [CovCost, CollisionCost, K, problem] = PRMCovC( Covs{i,1}-Covs{i,2}, Covs{i,2}, Covs{near_idx(j), 1}-Covs{near_idx(j), 2}, Covs{near_idx(j), 2}, N, param, Xbar, V, world );   
                if ~problem
                    plot3(meanTraj(1,:), meanTraj(2,:), meanTraj(3,:), 'color', [0.8 0.8 0.8], 'LineWidth', 1);
                    plot3(meanTraj(1,end), meanTraj(2,end), meanTraj(3,end), 'Marker','.','MarkerSize',12,'MarkerEdgeColor','[0.8500 0.3250 0.0980]');
                    n = n + 1;                    
                    ChildM(i, n) = near_idx(j);     % near_idx(j): (the index of) the nth child of node i      
                    EdgeTraj(i, n) = {meanTraj};
                    EdgeControlK(i, n) = {K};
                    EdgeControlV(i, n) = {V};
                    EdgesCostMean(i, n) = MCost;
                    EdgesCostCov(i, n) = CovCost;
                    EdgesCostColli(i, n) = CollisionCost;               
                end
            end
        end
    end
    points(i,7) = n;
end
plotWorld(world, dim); 

EdgesCost = EdgesCostMean + 1 * EdgesCostCov + 1 * EdgesCostColli;
% Each node has a state and a ChildNum
Nodes = points;

%% Finding Path
start_idx = 1;
goal_idx = 60;
weight = 1;
[count_vertex, Path, cost] = Astar(Nodes, ChildM, EdgesCost, start_idx, goal_idx, weight);

%% Plot Path
figure(2);
plotWorld(world, dim); 
pathCost = 0;
MC_num = 20;
PhatPrior0 = Covs{Path(1),1}-Covs{Path(1),2};
PtildePrior0 = Covs{Path(1),2};
rng(0); 

for mc = 1 : MC_num
    xbar0 = Nodes(Path(1), 1:6)';
    xhatPrior0_MC = mvnrnd(xbar0, PhatPrior0, 1)';
    xtildePrior0 = mvnrnd(zeros(6,1), PtildePrior0, 1)';
    x0_MC = xhatPrior0_MC + xtildePrior0;
    xhatPrior0_MC = x0_MC;      
    for i = 1:size(Path, 2)-1       
        xbar0 = Nodes(Path(i), 1:6)';
        PtildePrior0 = Covs{Path(i),2};
        idx = ChildM(Path(i), :) == Path(i+1);
        [x0_MC, xhatPrior0_MC] = MC_plan(x0_MC, xhatPrior0_MC, EdgeTraj{Path(i),idx}, PtildePrior0, EdgeControlV{Path(i),idx}, EdgeControlK{Path(i),idx}, true, param, world);
    end
end
figure(2)
for i = 1:size(Path, 2)-1
    Traj = EdgeTraj{Path(i),ChildM(Path(i), :) == Path(i+1)};
    plot3( Traj(1, :), Traj(2, :),  Traj(3, :), 'color', [0 0.4470 0.7410],'Linewidth', 1.5); hold on
    plot3(Nodes(Path(i), 1), Nodes(Path(i), 2), Nodes(Path(i), 3), 'Marker','.','MarkerSize',12,'MarkerEdgeColor','[0.8500 0.3250 0.0980]');
    pathCost = pathCost + EdgesCost(Path(i),ChildM(Path(i), :) == Path(i+1));   
end
plot3(Nodes(Path(end), 1), Nodes(Path(end), 2), Nodes(Path(end), 3), 'Marker','.','MarkerSize',12,'MarkerEdgeColor','[0.8500 0.3250 0.0980]')
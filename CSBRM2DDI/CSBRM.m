clear
addpath('CS');
addpath(genpath('C:\Program Files\Mosek\9.3'))
addpath(genpath('C:\Users\Dongliang\Documents\MATLAB\YALMIP-master'))
rng(0);
dim = 2; 
state_dim = 4;

% Initial state distribution (mean/covariance)
start_cord = [1, 1, 0, 0];
P0 = 0.8 * diag([0.01, 0.01, 0.01, 0.01]);

% Initial estimation error covariance
PtildePrior0 = 0.2 * P0;

% Initial estimated state covariance
PhatPrior0 = P0 - PtildePrior0;

goal_cord = [4, 3, 0, 0];

% Create random world
world = createKnownWorld([5,4],[0,0],dim);

param.dt = 0.2;
param.velavg = 1;

% Neighbor distance
r = 1;
% Maximum steplength
segmentLength = 1;

samples = 10;

points = start_cord;
Covs(1,1) = {P0};
Covs(1,2) = {PtildePrior0};

% figure(1); hold on
plot(start_cord(1), start_cord(2), 'Marker','s','MarkerSize',8,'MarkerEdgeColor','[0.8500 0.3250 0.0980]','MarkerFaceColor','[0.8500 0.3250 0.0980]'); hold on;
plotCovariance(start_cord, P0);
plotWorld(world, dim); 
 
%% Sample nodes
i = 1;
while i <= samples
        
        randomPoint = zeros(1, state_dim);
        min_dist = 0;

        while min_dist < 0.3            
            for j = 1:dim
                randomPoint(1, j) = (world.endcorner(j) - world.origincorner(j)-0.4) * rand + world.origincorner(j)+0.2;
            end
            % find the vertix that is closest to randomPoint (Eucl. dist. between positions)
            tmp = points(:, 1 : dim) - randomPoint(1 : dim);
            sqr_dist = sqr_eucl_dist(tmp, dim);
            [min_dist, idx] = min(sqr_dist);             
        end
        min_parent_idx = idx;
        
        Vect = randomPoint(1:dim)-points(idx,1:dim);
        Vect = Vect/norm(Vect);
        % find new_point that is within the range of idx
        if min_dist > segmentLength^2
            % generate a new point that is closest to randomPoint, segmentLength away from points(idx,1:dim)
            new_point(1 : dim) = points(idx, 1 : dim) + Vect * segmentLength;
        else
            new_point(1 : dim) = randomPoint(1 : dim);
        end

        if dim == 2
            new_point(dim + 1 : state_dim) = 0.2 * [-1 + (1 - (-1)) * rand, -1 + (1 - (-1)) * rand];     
        else
            if dim == 3
                new_point(dim + 1 : state_dim) = 0.2 * [-1 + (1 - (-1)) * rand, -1 + (1 - (-1)) * rand, -1 + (1 - (-1)) * rand];
            end
        end
        
        % check if the new_point is in collision
        if collision_point(new_point, world) == 0
            [meanTraj, V, MCost, N, Xbar]  = meanControl(points(idx, :)', new_point', param);  
            if ~MeanCollisionCheck(meanTraj(1:dim, :), world, dim)
                    [endP0, endPtildePrior0, CovCost, CollisionProb, K] = propagate( Covs{idx, 1}, Covs{idx, 2}, N, param, Xbar, V, world );                   
                    if CollisionProb < 0.1
                        i = i + 1;                         
    
                        points = [points; new_point];
                        idx_newpoint = size(points,1);
                        
                        Covs(idx_newpoint,1) = {endP0};
                        Covs(idx_newpoint,2) = {endPtildePrior0};
                    end                           
            end
        end
end

[points, Covs] = SampleMultipleVel(points, Covs, dim);
points = [points; goal_cord];
Covs(size(points,1), 1) = {0.015 * eye(4)};
Covs(size(points,1), 2) = {0.003 * eye(4)};

%% Graph construction 
% Neighbor distance
neb = 1.2;
% The ith row of ChildM contains the Children of node i 
ChildM = [];
% The EdgesCost(i,j) is the cost to go from node i to node j
EdgesCost = [];

% Adding edges
for i = 1:size(points, 1)   
    n = 0;
    
    vec = points(i,1:dim) - points(:,1:dim);  % Euclidean distances
    dist = sqrt(sum(vec.*vec, 2));     
%     dist(i) = dist(i) + 20;                 % Exclude the node itself as its' neighbor
    for k = 1:length(dist)                    % not connect nodes that are too close
        if dist(k) < 0.3
            dist(k) = 20;
        end
    end
    
    near_idx = find(dist < neb);
    for j = 1:size(near_idx, 1)
%         if collisioncheck(points(i,1:2), points(near_idx(j),1:2), Cspace)
        vel = points(near_idx(j),1:dim) - points(i,1:dim);
        vel = vel / norm(vel);
        if norm(points(near_idx(j),3:4))<0.00001 || norm(points(i,3:4))<0.00001 || sum(abs(vel - points(i,3:4)/norm(points(i,3:4)))) < 0.00001 
            [meanTraj, V, MCost, N, Xbar]  = meanControl( points(i,1:4)', points(near_idx(j),1:4)', param);
            if ~MeanCollisionCheck(meanTraj(1:2, :), world, dim)
                [CovCost, CollisionCost, K, problem] = PRMCovC( Covs{i,1}-Covs{i,2}, Covs{i,2}, Covs{near_idx(j), 1}-Covs{near_idx(j), 2}, Covs{near_idx(j), 2}, N, param, Xbar, V, world );   
                if ~problem
                    plot(meanTraj(1,:), meanTraj(2,:), 'color', [0 0.4470 0.7410], 'LineWidth', 1.5);
                    plot(meanTraj(1,end), meanTraj(2,end), 'Marker','.','MarkerSize',12,'MarkerEdgeColor','[0.8500 0.3250 0.0980]');
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
%         end
    end
    points(i,7) = n;
end
plotWorld(world, dim); 

EdgesCost = 1 * EdgesCostMean + EdgesCostCov + 0.5 * EdgesCostColli;

samp = size(points,1);

% Each node has a state and a ChildNum
Nodes = points;

%% Value iteration
% Utility function / Value function 
CostValue = 1000 * ones(samp, 1);

% Set starting point
StartIdx = 1;
CostValue(StartIdx) = 0;

% Convergence threshhold
epsn = 0.01;
delta = 0;

for i = 1:samp
    % The Value at the initial node is zero and not updated
    tempValue = [];
    for j = 1:Nodes(i, 7)        
        tempValue(j) = CostValue(i) + EdgesCost(i, j);
        if tempValue(j) < CostValue(ChildM(i, j))            
            if abs(CostValue(ChildM(i, j)) - tempValue(j)) > delta 
                delta = abs(CostValue(ChildM(i, j)) - tempValue(j));
            end
            CostValue(ChildM(i, j)) = tempValue(j);            
        end       
    end
end

iter = 1;

while delta > epsn
    iter = iter+1;
    delta = 0;
    for i = 1:samp
        tempValue = [];
        for j = 1:Nodes(i, 7)        
            tempValue(j) = CostValue(i) + EdgesCost(i, j);
            if tempValue(j) < CostValue(ChildM(i, j))                
                if abs(CostValue(ChildM(i, j)) - tempValue(j)) > delta 
                    delta = abs(CostValue(ChildM(i, j)) - tempValue(j));
                end
                CostValue(ChildM(i, j)) = tempValue(j);                
            end        
        end
    end
end

%% Finding Path

% Setting Goal point
GoalIdx = 36;
NextP = GoalIdx;
Path = GoalIdx;
tempValue = [];
while NextP ~= StartIdx
    Parents = sum( ChildM == NextP, 2 );
    parent_idx = find(Parents>=1);
    for i = 1:length(parent_idx)
        tempValue(i) = CostValue(parent_idx(i)) + EdgesCost(parent_idx(i), ChildM(parent_idx(i),:) == NextP);
    end
   
    [~, indx] = min( tempValue );
    NextP = parent_idx(indx);
    Path = [NextP, Path];
end


%% Plot Path
figure(10);
plotWorld(world, dim); 
pathCost = 0;

MC_num = 20;
PhatPrior0 = Covs{Path(1),1}-Covs{Path(1),2};
PtildePrior0 = Covs{Path(1),2};
rng(0); 

for mc = 1 : MC_num
    xbar0 = Nodes(Path(1), 1:4)';
    xhatPrior0_MC = mvnrnd(xbar0, PhatPrior0, 1)';
    xtildePrior0 = mvnrnd(zeros(4,1), PtildePrior0, 1)';
    x0_MC = xhatPrior0_MC + xtildePrior0;
    
    for i = 1:size(Path, 2)-1       
        xbar0 = Nodes(Path(i), 1:4)';
        PtildePrior0 = Covs{Path(i),2};
        [x0_MC, xhatPrior0_MC] = MC_plan(xbar0, x0_MC, xhatPrior0_MC, EdgeTraj{Path(i),ChildM(Path(i), :) == Path(i+1)}(1:3, :), PtildePrior0, EdgeControlV{Path(i),ChildM(Path(i), :) == Path(i+1)}, EdgeControlK{Path(i),ChildM(Path(i), :) == Path(i+1)}, false, param, world);
    end
end

for i = 1:size(Path, 2)-1
    plot3( EdgeTraj{Path(i),ChildM(Path(i), :) == Path(i+1)}(1, :), EdgeTraj{Path(i),ChildM(Path(i), :) == Path(i+1)}(2, :),  EdgeTraj{Path(i),ChildM(Path(i), :) == Path(i+1)}(3, :), 'color', [0 0.4470 0.7410],'Linewidth', 1.5); hold on
    plot3(Nodes(Path(i), 1), Nodes(Path(i), 2), Nodes(Path(i), 3), 'Marker','.','MarkerSize',12,'MarkerEdgeColor','[0.8500 0.3250 0.0980]');
    pathCost = pathCost + EdgesCost(Path(i),ChildM(Path(i), :) == Path(i+1));   
end
plot3(Nodes(Path(end), 1), Nodes(Path(end), 2), Nodes(Path(end), 3), 'Marker','.','MarkerSize',12,'MarkerEdgeColor','[0.8500 0.3250 0.0980]')
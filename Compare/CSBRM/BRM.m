% %% Generate Graph
% addpath('CS');
% rng(0);
% global Cspace
% 
% figure(1); hold on
% Cspace = im2bw(imread('E5.bmp')); % input map read from a bmp file. for new maps write the file name here
% start = [200,40]; % start position in Y, X format
% goal = [340,463]; % goal position in Y, X format
% k = 16; % number of points in the PRM
% display = true; % display processing of nodes
% 
% if ~checkpoint(start,Cspace), error('source lies on an obstacle or outside map'); end
% if ~checkpoint(goal,Cspace), error('goal lies on an obstacle or outside map'); end
% 
% imshow(Cspace);
% rectangle('position',[1 1 size(Cspace)-1],'edgecolor','k')
% hold on
% points = [start;goal]; % start and goal taken as additional vertices in the path planning to ease planning of the robot
% if display, rectangle('Position',[points(1,2)-6,points(1,1)-6,12,12],'Curvature',[1,1],'FaceColor','r'); end
% if display, rectangle('Position',[points(2,2)-6,points(2,1)-6,12,12],'FaceColor','blue'); end
% 
% % Signal positions
% Signal = [400, 10; 477, 380; 260, 12; 250, 460; 489, 110; 480, 250];
% for i = 1:size(Signal,1)
%     plot(Signal(i,2), Signal(i,1), 'kp', 'MarkerSize', 5, 'LineWidth', 1.5);
% end
% 
% param.dt = 0.2;
% param.vv = 0.2;
% param.vp = 0.1;
% param.vCons = 0.01;
% param.velavg = 4;
% param.scale = 20;
% param.neb = 200;    % Neighbor distance
% 
% [points, velRand, Covs] = sampleNodes(Cspace, k, points, Signal, display, param);
% 
% figure(20); 
% imshow(Cspace);
% rectangle('position',[1 1 size(Cspace)-1],'edgecolor','k');
% hold on
% for i = 1:size(Signal,1)
%     plot(Signal(i,2), Signal(i,1), 'kp', 'MarkerSize', 5, 'LineWidth', 1.5); hold on
% end
% for i = 1:size(points,1)
%     rectangle('Position',[points(i,2)-5,points(i,1)-5,10,10],'Curvature',[1,1],'FaceColor','[0.5 0.5 0.5]'); hold on
% end
% if display, rectangle('Position',[points(3,2)-6,points(3,1)-6,12,12],'Curvature',[1,1],'FaceColor','r'); end
% if display, rectangle('Position',[points(6,2)-6,points(6,1)-6,12,12],'FaceColor','blue'); end
% 
% % The ith row of Edges contains the Children of node i 
% ChildM = [];
% % The EdgesCost(i,j) is the cost to go from node i to node j
% EdgesCost = [];
% 
% % Adding edges
% for i = 1:size(points) 
%     n = 0;
%     % Euclidean distances
%     vec = points(i,1:2) - points(:,1:2);
%     dist = sqrt(sum(vec.*vec, 2));
%     for k = 1:length(dist)     % not connect nodes that are too close
%         if dist(k) < 40
%             dist(k) = size(Cspace,1) * 2;
%         end
%     end
%     near_idx = find(dist < param.neb);
%     for j = 1:size(near_idx, 1)
%         if collisioncheck(points(i,1:2), points(near_idx(j),1:2), Cspace)
%             vel = points(near_idx(j),1:2) - points(i,1:2);
%             vel = vel / norm(vel);
%             if sum(abs(vel - velRand(i,:)/norm(velRand(i,:)))) < 0.00001                
%                 [meanTraj, V, MCost, N, Xbar]  = meanControl( [points(i,1:2), velRand(i,:)]', [points(near_idx(j),1:2), velRand(near_idx(j),:)]', param);
%                 if MeanCollisionCheck(param.scale * meanTraj(1:2, :), Cspace)
%                     [CovCost, CollisionCost, K, problem] = PRMCovC( Covs{i,1}, Covs{i,2}, Covs{near_idx(j), 1}, Covs{near_idx(j), 2}, N, Signal, param, Xbar, V, Cspace );   
%                     if ~problem
%                         figure(1)
%                         plot(param.scale * meanTraj(2,:), param.scale * meanTraj(1,:), 'color', 'g'); hold on
%                         n = n + 1;
%                         % near_idx(j): (the index of) the nth child of node i
%                         ChildM(i, n) = near_idx(j); 
%                         EdgeTraj(i, n) = {meanTraj};
%                         EdgeControlK(i, n) = {K};
%                         EdgeControlV(i, n) = {V};
%                         EdgesCostMean(i, n) = MCost;
%                         EdgesCostCov(i, n) = CovCost; 
%                         EdgesCostColli(i, n) = CollisionCost; 
%                     end
%                 end                  
%             end
%         end
%     end
%     points(i,3) = n;
% end
EdgesCost = EdgesCostMean + EdgesCostCov + 1 * EdgesCostColli;

% Sample points
samp = size(points, 1);
% Each node has three elements: x, y, ChildNum
Nodes = points;

% Utility function / Value function 
CostValue = zeros(samp, 1);

% Plot the graph
figure(10); 
imshow(Cspace); hold on
rectangle('position',[1 1 size(Cspace)-1],'edgecolor','k')
for i = 1:size(Signal,1)
    plot(Signal(i,2), Signal(i,1), 'kp', 'MarkerSize', 5, 'LineWidth', 1.5); hold on
end
for i = 1:size(points,1)
    rectangle('Position',[points(i,2)-5,points(i,1)-5,10,10],'Curvature',[1,1],'FaceColor','[0.5 0.5 0.5]'); hold on
end
if display, rectangle('Position',[points(3,2)-6,points(3,1)-6,12,12],'Curvature',[1,1],'FaceColor','r'); end
if display, rectangle('Position',[points(6,2)-6,points(6,1)-6,12,12],'FaceColor','blue'); end

for i = 1:samp
    for j = 1:Nodes(i, 3) 
        plot([Nodes(i, 2); Nodes(ChildM(i, j), 2)], [Nodes(i, 1); Nodes(ChildM(i, j), 1)], 'r', 'marker','o'); hold on
    end
end

%% Value Iteration

% Set goal point
GoalIdx = 5;

% Convergence threshhold
epsn = 0.001;
delta = 0;

for i = 1:samp
    % The Value at the goal node is zero and not updated
    if i ~= GoalIdx
        tempValue = [];
        for j = 1:Nodes(i, 3)
            tempValue(j) = EdgesCost(i, j) + CostValue(ChildM(i, j));
        end
        NewValue = min( tempValue );
        if abs(CostValue(i) - NewValue) > delta 
            delta = abs(CostValue(i) - NewValue);
        end
        CostValue(i) = NewValue;  
    end
end

iter = 1;

while delta > epsn
    iter = iter+1;
    delta = 0;
    for i = 1:samp
        if i ~= GoalIdx
            tempValue = [];
            for j = 1:Nodes(i, 3)
                tempValue(j) = EdgesCost(i, j) + CostValue(ChildM(i, j));
            end
            NewValue = min( tempValue );
            if abs(CostValue(i) - NewValue) > delta 
                delta = abs(CostValue(i) - NewValue);
            end
            CostValue(i) = NewValue;  
        end
    end
end

%% Finding Path

% Starting point, the first node is the starting point
NextP = 2;
Path = 2;
tempValue = [];
while NextP ~= GoalIdx    
    for i = 1:Nodes(NextP, 3)
        tempValue(i) = EdgesCost(NextP, i) + CostValue(ChildM(NextP, i));
    end
    [~,Idx] = min(tempValue);
    NextP = ChildM(NextP, Idx);
    Path = [Path, NextP];
end

% Plot Path
figure(10);
pathCost = 0;
for i = 1:size(Path, 2)-1
    [meanTraj, V, MCost, N, Xbar]  = meanControl( [points(Path(i),1:2), velRand(Path(i),:)]', [points(Path(i+1),1:2), velRand(Path(i+1),:)]', param);
    [CovCost, CollisionCost, K, problem] = PRMCovC( Covs{Path(i),1}, Covs{Path(i),2}, Covs{Path(i+1), 1}, Covs{Path(i+1), 2}, N, Signal, param, Xbar, V, Cspace );   
end
for i = 1:size(Path, 2)-1
    plot(param.scale * EdgeTraj{Path(i),ChildM(Path(i), :) == Path(i+1)}(2, :), param.scale * EdgeTraj{Path(i),ChildM(Path(i), :) == Path(i+1)}(1, :), 'color', 'g','Linewidth', 2); hold on
    pathCost = pathCost + EdgesCost(Path(i),ChildM(Path(i), :) == Path(i+1));   
end




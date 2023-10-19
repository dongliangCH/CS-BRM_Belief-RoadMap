%% Value iteration
% Utility function / Value function 
CostValue = zeros(samp, 1);

% Set goal point
GoalIdx = 10;

% Convergence threshhold
epsn = 0.01;
delta = 0;

for i = 1:samp
    % The Value at the goal node is zero and not updated
    if i ~= GoalIdx
        tempValue = [];
        if Nodes(i, 5) == 0
            CostValue(i) = 1000;
        else    
            for j = 1:Nodes(i, 5)
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

iter = 1;

while delta > epsn
    iter = iter+1;
    delta = 0;
    for i = 1:samp
        if i ~= GoalIdx
            tempValue = [];
            if Nodes(i, 5) == 0
                CostValue(i) = 1000;
            else
                for j = 1:Nodes(i, 5)
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
end

%% Finding Path

% Starting point, the first node is the starting point
NextP = 1;
Path = 1;
tempValue = [];
while NextP ~= GoalIdx    
    for i = 1:Nodes(NextP, 5)
        tempValue(i) = EdgesCost(NextP, i) + CostValue(ChildM(NextP, i));
    end
    [~,Idx] = min(tempValue);
    NextP = ChildM(NextP, Idx);
    Path = [Path, NextP];
end





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
GoalIdx = 10;
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
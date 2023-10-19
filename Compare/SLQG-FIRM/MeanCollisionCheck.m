function feasible = MeanCollisionCheck(posi,map)
feasible=true;
for i = 1:2:size(posi,2)
    posCheck=posi(:, i);
    if ~(checkpoint(ceil(posCheck),map) && checkpoint(floor(posCheck),map) && ... 
            checkpoint([ceil(posCheck(1)) floor(posCheck(2))],map) && checkpoint([floor(posCheck(1)) ceil(posCheck(2))],map))
        feasible=false;break;
    end
end
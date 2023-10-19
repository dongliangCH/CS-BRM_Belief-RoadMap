function [world, NumObstacles] = createKnownWorld(endcorner, origincorner, dim)

  if dim == 2
      
    NumObstacles = 3;
  % check to make sure that the region is nonempty
    if (endcorner(1) <= origincorner(1)) || (endcorner(2) <= origincorner(2))
      disp('Not valid corner specifications!')
      world=[];
  % create world data structure
    else
    world.NumObstacles = NumObstacles;
    world.endcorner = endcorner;
    world.origincorner = origincorner;

    % create NumObstacles
    
        x = 1;
        y = 2;
        a = 1;
        b = 1;
        world.ox(1) = x;
        world.oy(1) = y;
        world.oa(1) = a;
        world.ob(1) = b;

        x = 3;
        y = 0.5;
        a = 1.5;
        b = 1.5;
        world.ox(2) = x;
        world.oy(2) = y;
        world.oa(2) = a;
        world.ob(2) = b;
        
        x = 1.5;
        y = 2.5;
        a = 1;
        b = 1;
        world.ox(3) = x;
        world.oy(3) = y;
        world.oa(3) = a;
        world.ob(3) = b;
        
    end

  elseif dim == 3

    NumObstacles = 6;
    % check to make sure that the region is nonempty
    if (endcorner(1) <= origincorner(1)) || (endcorner(2) <= origincorner(2)) || (endcorner(3) <= origincorner(3))
      disp('Not valid corner specifications!')
      world=[];

    % create world data structure
    else
    world.NumObstacles = NumObstacles;
    world.endcorner = endcorner;
    world.origincorner = origincorner;

    % create NumObstacles
        
        x = 0;
        y = 1.5;
        z = 1;
        a = 1.5;
        b = 1.5;
        c = 2;
        world.ox(1) = x;
        world.oy(1) = y;
        world.oz(1) = z;
        world.oa(1) = a;
        world.ob(1) = b;
        world.oc(1) = c;
        
        x = 2.5;
        y = 0.4;
        z = 1; 
        a = 1;
        b = 1.2;
        c = 2;
        world.ox(3) = x;
        world.oy(3) = y;
        world.oz(3) = z;
        world.oa(3) = a;
        world.ob(3) = b;
        world.oc(3) = c;

        x = 4;
        y = 2.2;
        z = 1;
        a = 1;
        b = 0.8;
        c = 2;
        world.ox(2) = x;
        world.oy(2) = y;
        world.oz(2) = z;
        world.oa(2) = a;
        world.ob(2) = b;
        world.oc(2) = c;
        
        x = 4.6;
        y = 0;
        z = 1;
        a = 1.5;
        b = 1.5;
        c = 2;
        world.ox(6) = x;
        world.oy(6) = y;
        world.oz(6) = z;
        world.oa(6) = a;
        world.ob(6) = b;
        world.oc(6) = c;  
        
        x = 6.5;
        y = 2;
        z = 1;
        a = 1.5;
        b = 1;
        c = 2;
        world.ox(5) = x;
        world.oy(5) = y;
        world.oz(5) = z;
        world.oa(5) = a;
        world.ob(5) = b;
        world.oc(5) = c;
        
%         x = 6;
%         y = 3.5;
%         z = 1;
%         a = 1.5;
%         b = 1;
%         c = 1.5;
%         world.ox(4) = x;
%         world.oy(4) = y;
%         world.oz(4) = z;
%         world.oa(4) = a;
%         world.ob(4) = b;
%         world.oc(4) = c;
%         x = 8;
%         y = 2;
%         z = 0;
%         a = 1;
%         b = 1;
%         c = 2;
%         world.ox(7) = x;
%         world.oy(7) = y;
%         world.oz(7) = z;
%         world.oa(7) = a;
%         world.ob(7) = b;
%         world.oc(7) = c;
        
     end
   end
end
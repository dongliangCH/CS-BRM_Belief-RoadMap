function plotWorld(world,dim)
  % the first element is the north coordinate
  % the second element is the south coordinate
  if dim ==2

  axis([world.origincorner(1),world.endcorner(1),...
      world.origincorner(2), world.endcorner(2)]);
  hold on

  for i=1:world.NumObstacles
      
      X = [world.ox(i)    world.ox(i)+world.oa(i)    world.ox(i)+world.oa(i)  world.ox(i)  world.ox(i)];
      Y = [world.oy(i)    world.oy(i)                world.oy(i)+world.ob(i)  world.oy(i)+world.ob(i)   world.oy(i)];
      fill(X, Y, [0.5 0.5 0.5]);

  end
  
  plot([world.origincorner(1) world.endcorner(1)  world.endcorner(1) world.origincorner(1) world.origincorner(1)], [world.origincorner(2) world.origincorner(2) world.endcorner(2) world.endcorner(2) world.origincorner(2)], 'color', 'k', 'linewidth', 2)
  axis off
%   plot(2, 5, 'Marker','^','MarkerSize',12,'MarkerEdgeColor','[0.8500 0.3250 0.0980]','MarkerFaceColor','[0.8500 0.3250 0.0980]');
%   plot(8, 5, 'Marker','^','MarkerSize',12,'MarkerEdgeColor','[0.8500 0.3250 0.0980]','MarkerFaceColor','[0.8500 0.3250 0.0980]');
%   plot(4.5, 0, 'Marker','v','MarkerSize',12,'MarkerEdgeColor','[0.8500 0.3250 0.0980]','MarkerFaceColor','[0.8500 0.3250 0.0980]');
  
%   X = [3  5  5  3  3];
%   Y = [0  0  2  2  0];
%   fill(X, Y, 'blue');
% 
%   X = [3  5  5  3  3];
%   Y = [3.5  3.5  5  5  3.5];
%   fill(X, Y, 'blue');
   
  elseif dim ==3
  axis([world.origincorner(1),world.endcorner(1),...
      world.origincorner(2), world.endcorner(2),...
      world.origincorner(3), world.endcorner(3)]);
  hold on

  for i=1:world.NumObstacles
      
      vert = [world.ox(i)               world.oy(i)              world.oz(i);
              world.ox(i)+world.oa(i)   world.oy(i)              world.oz(i);
              world.ox(i)+world.oa(i)   world.oy(i)+world.ob(i)  world.oz(i);
              world.ox(i)               world.oy(i)+world.ob(i)  world.oz(i);
              world.ox(i)               world.oy(i)              world.oz(i)+world.oc(i);
              world.ox(i)+world.oa(i)   world.oy(i)              world.oz(i)+world.oc(i);
              world.ox(i)+world.oa(i)   world.oy(i)+world.ob(i)  world.oz(i)+world.oc(i);
              world.ox(i)               world.oy(i)+world.ob(i)  world.oz(i)+world.oc(i)];
      face = [1 2 3 4;
              1 2 6 5;
              2 3 7 6;
              1 4 8 5;
              5 8 7 6;
              3 4 8 7];
      patch('Faces',face,'Vertices',vert,'FaceColor','[0.4 0.4 0.4]','EdgeColor','[0.4 0.4 0.4]','FaceAlpha',.8)
%       material shiny;
%       alpha('color');
%       alphamap('rampdown');
      
  end
%   axis off
%   plot3(2, 5, 2.8, 'Marker','^','MarkerSize',12,'MarkerEdgeColor','[0.8500 0.3250 0.0980]','MarkerFaceColor','[0.8500 0.3250 0.0980]');
%   plot3(8, 5, 2.8, 'Marker','^','MarkerSize',12,'MarkerEdgeColor','[0.8500 0.3250 0.0980]','MarkerFaceColor','[0.8500 0.3250 0.0980]');
%   plot3(4.5, 0, 2.8, 'Marker','^','MarkerSize',12,'MarkerEdgeColor','[0.8500 0.3250 0.0980]','MarkerFaceColor','[0.8500 0.3250 0.0980]');
  
  end
%   axis equal
%   ylim([0,3]);
%   zlim([1.5,3]);  
%   xl = xlabel('$x  (m)$','Interpreter','LaTeX');
%   yl = ylabel('$y  (m)$','Interpreter','LaTeX');
%   zl = zlabel('$z  (m)$','Interpreter','LaTeX');
%   set(xl,'FontSize',18);
%   set(yl,'FontSize',18);
%   set(zl,'FontSize',18);
%   set(gca,'FontSize',16,'FontName','Times');
  
end

% goplotpso4demo.m
% simple graphing prog for use with PSO demo script

% Brian Birge
% Rev 1.0
% 1/1/3
set(gcf,'Position',[693    28   586   500]); % this is computer dependent
set(gcf,'Doublebuffer','on');

subplot(1,2,1)
 ptsx=30;
 ptsy=30;
 rngx=max(pbest(:,1))-min(pbest(:,1));
 rngy=max(pbest(:,2))-min(pbest(:,2));
 clear x y f6val
 x=min(pbest(:,1)):rngx/(ptsx-1):max(pbest(:,1));
 y=min(pbest(:,2)):rngy/(ptsy-1):max(pbest(:,2));
 
 x=-20:.5:20;
 y=-20:.5:20;
 ptsx=length(x);
 ptsy=length(y);
 
 for i=1:ptsx
     for j=1:ptsy
         f6val(i,j)=f6([x(i),y(j)]);
      %   f6val(i,j)=f6mod([x(i),y(j)]);                  
      %   f6val(i,j)=f6_linear_dyn([x(i),y(j)]);         
      %   f6val(i,j)=f6_spiral_dyn([x(i),y(j)]);                  
      %   f6val(i,j)=f6_bubbles_dyn([x(i),y(j)]);         
      %   f6val(i,j)=ackley([x(i),y(j)]);
      %   f6val(i,j)=tripod([x(i),y(j)]);         
      %   f6val(i,j)=alpine([x(i),y(j)]);         
      %   f6val(i,j)=NDparabola([x(i),y(j)]);         
      %   f6val(i,j)=Griewank([x(i),y(j)]);         
      %   f6val(i,j)=Rastrigin([x(i),y(j)]);         
      %   f6val(i,j)=Rosenbrock([x(i),y(j)]);         
      %   f6val(i,j)=DeJong_f2([x(i),y(j)]);         
      %   f6val(i,j)=DeJong_f4([x(i),y(j)]);         
      %   f6val(i,j)=DeJong_f4_noise([x(i),y(j)]);         
      %   f6val(i,j)=DeJong_f3([x(i),y(j)]);         
      %   f6val(i,j)=foxhole([x(i),y(j)]);         
     end
 end
 surf(x,y,f6val')
% contour(x,y,f6val',1,'k');
 %axis([min(x) max(x) min(y) max(y) 0 1]);
 colormap(copper),shading interp
 %grid on
 view([-22.5,76]) % good 3d view
% view([-25.5,-20]) % good 3d view 
% view([1.2650e+002  7.4000e+001]) % good 3d view 
% view([3.8500e+001  4.2000e+001]) % good 3d view 
% view(2)
%view([60,16]);
%view(3)
% view([60.5,-16]) % uncomment this for alternate 3D view
% view([0,-90])   % underside 2D view
 xlabel('x')
 ylabel('y')
 zlabel('cost')
 title([num2str(gbestval),' = ',functname,'([',num2str(gbest(1)),', ',num2str(gbest(2)),'])'])
 hold on     
 plot3(pos(:,1),pos(:,2),out,'b.','Markersize',7)
 plot3(pbest(:,1),pbest(:,2),pbestval,'g.','Markersize',7);
 plot3(gbest(1),gbest(2),gbestval,'r.','Markersize',25);
 %plot3(gbest(1),gbest(2),0,'r.','Markersize',10);
% try
%   axis([min(pbest(:,1)) max(pbest(:,1)) ...
%         min(pbest(:,2)) max(pbest(:,2)) ...
%         min(min(f6val)) max(max(f6val))]);
% catch
%   axis tight  
% end
axis([-20 20 -20 20 -1 1])
%axis([-1 20 -1 20 0 1])
 %axis equal
 hold off
 drawnow
 
subplot(1,2,2)
 semilogy(tr(find(~isnan(tr)))),xlabel('epoch');,ylabel('gbest');
 if trelea==0
       title(['ps=',num2str(ps),', Inertia wt=',num2str(iwt(length(iwt))),...
               ', Common PSO, ',num2str(D),'D']);
     elseif trelea==2 | trelea==1
       title(['ps=',num2str(ps),', Trelea mod=',num2str(trelea),...
          ', ',num2str(D),'D']);
     elseif trelea==3
       title(['ps=',num2str(ps),', Clerc Type 1", \chi = ',num2str(chi),...
          ', ',num2str(D),'D']);
     end
 grid on
 axis tight
 drawnow
 
% evalin('base','tempmvframe=getframe(gcf);');
% evalin('base','aviobj=addframe(aviobj,tempmvframe);');
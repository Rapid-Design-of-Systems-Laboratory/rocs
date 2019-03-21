% goplotpso.m
% simple graphing prog called by pso_Trelea_vectorized.m
% Blue dots are current positions, green dots are personal bests, and red
% is global best

% Brian Birge
% Rev 1.0
% 1/1/3

set(gcf,'Doublebuffer','on');
%set(gcf,'Position',[3 30 683 497]);
set(gcf,'Position',[694    32   584   476]);
    subplot(2,1,1)
     semilogy(tr(1:length(tr))),xlabel('epoch');,ylabel('gbest');
     if trelea==0
       title(['ps=',num2str(ps),', Inertia wt=',num2str(iwt(i)),...
               ', Common PSO, ',num2str(D),'D, cost=',num2str(gbestval)]);
     elseif trelea==2 | trelea==1
       title(['ps=',num2str(ps),', Trelea mod=',num2str(trelea),...
          ', ',num2str(D),'D, cost=',num2str(gbestval)]);
     elseif trelea==3
       title(['ps=',num2str(ps),', Clerc Type 1", \chi = ',num2str(chi),...
          ', ',num2str(D),'D, cost=',num2str(gbestval)]);
     end
     hold on
     semilogy(1:te,ones(size(1:te))*errgoal,'r-.');
     grid on
     hold off
     drawnow
     
    subplot(2,1,2)
%     beststrg=[];
%     for plotloopD=1:D
%         beststrg=[beststrg,num2str(gbest(plotloopD)),', '];
%     end
     plot(pos(:,1),pos(:,D),'b.')     
     hold on
%    % for asd=1:ps
%    %   plot([pos(asd,1),pbest(asd,1)],[pos(asd,D),pbest(asd,D)],'g-')
%    % end
     plot(pbest(:,1),pbest(:,D),'g.');
%     
%%    % interconnectedness plot
%%     conmatrix=ones(ps,ps);
%%     conmatrix=triu(conmatrix,1);  % all ones above main diagonal
%%     conindx=find(conmatrix~=0);
%%     [i_idx,j_idx]=ind2sub([ps,ps],conindx);  
%%     for asd2=1:length(conindx)
%%         plot([pos(i_idx,1),pos(j_idx,1)],[pos(i_idx,D),pos(j_idx,D)],'b');
%%     end
%     
     xlabel('pos dim 1');
     ylabel(['pos dim ',num2str(D)]);
%     title(['Epoch: ',num2str(i),', ',num2str(gbestval),' = ',functname,...
%           '( ',beststrg(1:length(beststrg)-2),' )']);
     grid on
     plot(gbest(1),gbest(D),'r*');    
     hold off
%     if rem(cnt,10)==0 | cnt==1
%         %tmpax=[min(pos(:,1)) max(pos(:,1)) min(pos(:,2)) max(pos(:,2))];
%         tmpax=[min(pbest(:,1)) max(pbest(:,1)) min(pbest(:,2)) max(pbest(:,2))];         
%     end
     
    % try 
    %    axis(tmpax);
    % catch
  %      axis tight;
    % end
  %   drawnow     
     
 %   subplot(2,2,4)
 %    plot(vel(:,1),vel(:,D),'b.');
 %    xlabel('vel dim 1');
  %   ylabel(['vel dim ',num2str(D)]);
  %   grid on
  %   hold on
  %   plot(vel(idx1,1),vel(idx1,D),'r*');
  %   hold off
% %    axis([-10 10 -10 10]);
  %   axis tight
  %   axis equal
  %   drawnow
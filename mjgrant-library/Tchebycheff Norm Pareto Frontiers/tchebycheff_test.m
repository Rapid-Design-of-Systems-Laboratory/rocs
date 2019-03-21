%TCHEBYCHEFF_TEST provides a proof of concept example of how the simple 
%additive weighting (L1 norm) scheme for a multi-attribute problem is
%flawed compared to the Tchebycheff (infinity) norm, which will identify
%all Pareto-optimal points under varying weights -- even for concave Pareto
%frontiers.
%
%Jarret Lafleur, NASA JSC/EG5, May 2010

%Initialize rand seed (for consistency)
rand('state',0);

y = rand(5000,2);   %Create 5000 possible designs, rated in 2 attributes
exclval = 0.75;     %Radius of circle to exclude

%Remove designs to create a concave Pareto frontier of radius "exclval"
ym1 = y-1;
Rs  = ym1.^2;
sums = sum(Rs,2);
goodinds = find(sums>exclval.^2);
ygood = y(goodinds,:);

%figure(1); plot(y(:,1),y(:,2),'+'); title('Before Filter');
%axis([0 1 0 1]); axis square
figure(2); plot(ygood(:,1),ygood(:,2),'+'); title('Full Data Set');
axis([0 1 0 1]); axis square;
xlabel('Objective #1'); ylabel('Objective #2');

weightlist = linspace(0,1,1000);
for w = 1:length(weightlist)
    w1 = weightlist(w); w2=1-w1;
    wts = [w1 w2]';
    
    %Compute simple additive weighting scores.
    %SAW_scores = w1*ygood(:,1) + w2*ygood(:,2);
    SAW_scores = ygood*wts;     %More general
    [SAW_opts(w),SAW_windex(w)] = max(SAW_scores);
    
    %Compute Tchebycheff norm scores.
    %Assumes the ideal criterion vector is (1,1), which is true here.
    %Tcheb_scores = -max(w1*abs(1-(ygood(:,1))),w2*abs(1-(ygood(:,2))));
    Tcheb_scores = -max(abs(1-ygood).*repmat(wts',size(ygood,1),1),[],2);   %More general
    [Tcheb_opts(w),Tcheb_windex(w)] = max(Tcheb_scores);
end

SAW_paretos = ygood(SAW_windex,:);
Tcheb_paretos = ygood(Tcheb_windex,:);
figure(3); plot(SAW_paretos(:,1),SAW_paretos(:,2),'+');
axis([0 1 0 1]); axis square; title('Simple Additive Weighting Pareto-Optimal Points');
xlabel('Objective #1'); ylabel('Objective #2');

figure(4); plot(Tcheb_paretos(:,1),Tcheb_paretos(:,2),'+');
axis([0 1 0 1]); axis square; title('Weighted Tchebycheff Norm Pareto-Optimal Points');
xlabel('Objective #1'); ylabel('Objective #2');

%x=0:0.01:1;
%y=x;
%[xx,yy]=meshgrid(x,y);
%Tscores = -max(w1*abs(1-(xx)),w2*abs(1-(yy)));
%figure(209)
%pcolor(xx,yy,Tscores); axis square
big=max(fit); little=6000;%min(fit);
range=big-little;

for n=1:length(fit)
    if isnan(fit(n))~=1 & fit(n)>6000
        col=[(fit(n)-little)/range 0 1-(fit(n)-little)/range];
        plot3(popb1(n),popgamma(n),fit(n),'.','Color',col)
        hold on
    end
end
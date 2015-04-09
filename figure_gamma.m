% Load all states of a chain
f1 = load('results/odemcmc-jakstat-log-1-2000000-50000-1-1.mat');
t0 = 5e4;
nn = 2e6;
fx = f1.x((t0+1):(t0+nn),:);
etamax = 4;
etas = ceil(logspace(0,etamax,20));
[n,m] = size(fx);
Vf = (1/n)*(sum(fx.^2)) - ((1/n)*sum(fx)).^2;
fm = mean(fx);
for j=1:m
    for i=1:length(etas)
        rho(i,j) = sum((fx(1:(n-etas(i)),j)-fm(j)).*(fx((etas(i)+1):n,j)-fm(j)))/(n-etas(i));
        gammas(i,j) = 1-(abs(rho(i,j))/Vf(j)).^(1/etas(i));
    end
end
[gamma,gammas2,ming,eta2] = getGammaIter(fx);

markers = {'x','^','v','o','s'};
figure; hold on;
for j=1:m
    plot(etas,gammas(:,j),['k-' markers{j}]);
end

markersize = 8;
plot([eta2(end),eta2(end)],[1e-4,1e-1],'k--','markers',markersize);

set(gca,'xscale','log');
set(gca,'yscale','log');
legend({'k_1','k_2','k_3','k_4'});
xlim([1 10^etamax]);
xlabel('$\eta$','Interpreter','LaTex','fontsize',13);
ylabel('$\hat{\gamma}^* := 1- (|\hat{\rho}_{\eta}| / V_{\pi}{f})^{1/\eta}$','Interpreter','LaTex','fontsize',13);
annotation('textarrow',[0.6 0.69],[0.3 0.35],'String','Final \eta')
set(findall(gcf,'type','text'),'FontSize',16)
set(gca,'FontSize',14)
saveas(gcf,'figs/figure_gamma.eps','eps');

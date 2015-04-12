% Figure for distibution of error rate and bound for fixed sample size test

clear variables;

reload = false;
if reload
	nRuns = 1000;
	N = 2e6;
	% Load parameters for each chain
	load('results/gammas_2e6_from2e5.mat');

	epsilon = 0.01;
	mtrue = 0.88746;
	Ns = ceil(logspace(3,6,15));
	deltas = [0.02,0.05,0.1];
	fn = zeros(length(deltas),length(Ns));
	fp = zeros(length(deltas),length(Ns));
	for k=1:nRuns
		fname = sprintf('/home/beni/data/phd/mcmchyp/odemcmc-jakstat-1-2000000-50000-1-%d.mat',k);
		fprintf('Reading %s\n',fname);
		f = load(fname);
		sfx = cumsum(f.f);
		gamma = gammas(k);
		for j=1:length(deltas)
			delta = deltas(j);
			r = mtrue-delta;
			for i=1:length(Ns)
				ns = Ns(i);
				if sfx(ns) <= ns*r
					fn(j,i) = fn(j,i) + 1;
				end
			end
		end
	end
	fn = fn/nRuns;
else
	load('results/figure_fix_dat.mat');
end


%% Plotting
fh = figure; hold on;
fn(fn==0) = 0.5*(1/nRuns); % Set a small value below ylimit to avoid gap on log scale plot

markersize = 8;
plot(Ns,fn(1,:),'x','markers',markersize);
plot(Ns,fn(2,:),'o','markers',markersize);
plot(Ns,fn(3,:),'^','markers',markersize);
plot([1],nan,'--');
plot([1],nan,'-');
plot(Ns,fn(1,:),'-','markers',markersize);
plot(Ns,fn(2,:),'-','markers',markersize);
plot(Ns,fn(3,:),'-','markers',markersize);

for j=1:length(deltas)
	for i=1:length(Ns)
		Nfix(j,i) = mean(exp(-Ns(i)*gammas*deltas(j)^2));
	end
end
plot(Ns,Nfix(1,:),'x--','markers',markersize);
plot(Ns,Nfix(2,:),'o--','markers',markersize);
plot(Ns,Nfix(3,:),'^--','markers',markersize);

legends = {'\delta=0.02','\delta=0.05','\delta=0.1','\epsilon_n','E_n'};

xlabel('Number of samples (n)','FontSize',12);
ylabel('Error rate and bound (E_n and \epsilon_n)','FontSize',12);
ylim([1/nRuns,1]);
set(findall(gcf,'type','text'),'FontSize',16)
set(gca,'FontSize',14)
set(gca,'Xscale','log');
set(gca,'Yscale','log');
legend(legends,'Location','SouthWest');
saveas(fh, 'figs/figure_fix','eps');

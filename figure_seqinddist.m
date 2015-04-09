% Figure for distibution of stopping times for sequential test with
% indifference region
clear variables

reload = false;
if reload
	nRuns = 1000;
	N = 2e6;
	load('results/gammas_2e6_from2e5.mat');

	rs = [0.1, 0.5, 0.7, 0.8, 0.95];
	epsilon = 0.01;
	delta = 0.05;

	T = N*ones(length(rs),nRuns);
	for k=1:nRuns
		gamma = gammas(k);
		fname = sprintf('results/odemcmc-jakstat-1-2000000-50000-1-%d.mat',k);
		fprintf('Reading %s\n',fname);
		f = load(fname);
		sfx = cumsum(f.f);
		N = length(f.f);
		Nfix(k) = log(1/epsilon)/(gamma*delta^2);
		for i=1:length(rs)
			r = rs(i);
			M = log(epsilon*gamma*delta^2/2)/ (-2*gamma*delta - (gamma*delta^2)/(1-r));
			L = (1:N)'*r - M;
			U = (1:N)'*r + M;
			tj = find(sfx<L | sfx>U,1,'first');
			if ~isempty(tj)
				T(i,k) = tj;
			else
				T(i,k) = nan;
			end
		end
	end
else
	load('results/figure_seqinddist_dat.mat');
end

%% Plotting
fh = figure; hold on;
colors = {'x','^','o','v','s'};

markersize = 8;
for i=1:length(rs)
	[f,x]=hist(T(i,:),200);
	f = f./nRuns;
	cmf = cumsum(f);
	xend = find(cmf>0.9,1,'first');
	xtick = x(1:round(xend/10):xend);
	ytick = cmf(1:round(xend/10):xend);
	plot(xtick,ytick,colors{i},'markers',markersize);
	legends{i} = ['r=' num2str(rs(i))];
end

[f,x]=hist(Nfix,200);
cmf = cumsum(f)/nRuns;
plot(x,cmf,'--');

for i=1:length(rs)
	[f,x]=hist(T(i,:),200);
	f = f./nRuns;
	cmf = cumsum(f);
	plot([0 x],[0 cmf]);
end

legend([legends, 'fixed'],'location','SouthEast');
xlim([0,1e6]);
ylim([0 1])
grid off;
xlabel('Number of samples','FontSize',12);
ylabel('Ratio of chains stopped','FontSize',12);
set(findall(gcf,'type','text'),'FontSize',16);
set(gca,'FontSize',14)
saveas(fh,'figs/figure_seqinddist','eps');

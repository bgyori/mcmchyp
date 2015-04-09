clear variables;
reload = false;
if reload
	nRuns = 1000;
	N = 2e6;
	load('results/gammas_2e6_from2e5.mat');

	rs = [0.1, 0.5, 0.7, 0.8, 0.95];
	epsilon = 0.01;

	T = N*ones(length(rs),nRuns);
	for k=1:nRuns
		gamma = gammas(k);
		fname = sprintf('results/odemcmc-jakstat-1-2000000-50000-1-%d.mat',k);
		fprintf('Reading %s\n',fname);
		load(fname);
		sfx = cumsum(f);
		N = length(f);
		for i=1:length(rs)
			r = rs(i);
			g = sqrt((1:N)/gamma * log(1/epsilon) + 1 + 2*log(1:N))';
			L = (1:N)'*r - g;
			U = (1:N)'*r + g;
			tj = find(sfx<L | sfx>U,1,'first');
			if ~isempty(tj)
				T(i,k) = tj;
			else
				T(i,k) = nan;
			end
		end
	end
else
	load('results/fig_seqnoinddist_dat.mat');
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
for i=1:length(rs)
	[f,x]=hist(T(i,:),200);
	f = f./nRuns;
	cmf = cumsum(f);
	plot([0 x],[0 cmf]);
end

legend(legends,'location','SouthEast');
xlim([0,1e6]);
ylim([0 1])
grid off;
xlabel('Number of samples','FontSize',12);
ylabel('Ratio of chains stopped','FontSize',12);
set(findall(gcf,'type','text'),'FontSize',16);
set(gca,'FontSize',14)
saveas(fh,'figs/figure_seqnoinddist','eps');

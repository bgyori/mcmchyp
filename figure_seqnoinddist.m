clear variables;
reload = false;
if reload
	nRuns = 1000;
	load('results/gammas_2e6_from2e5.mat');

	rs = [0.1, 0.5, 0.7, 0.8, 0.95];
	epsilon = 0.01;
    xi = 0.3;

	T = ones(length(rs),nRuns);
	for k=1:nRuns
		gamma = gammas(k);
		fname = sprintf('~/data/phd/mcmchyp/odemcmc-jakstat-1-2000000-50000-1-%d.mat',k);
		fprintf('Reading %s\n',fname);
		load(fname);
		sfx = cumsum(f);
		N = length(f);
		n0 = floor(100/gamma);
		imax = floor(log(N/n0)/log(1+xi));
		ntest = floor(n0*(1+xi).^(1:imax));
            
		g = sqrt((ntest/gamma) .* (log(1/epsilon) + 1 + 2*log(1:imax)))';
		for i=1:length(rs)
			r = rs(i);
			
			L = ntest'*r - g;
			U = ntest'*r + g;
			tj = find(sfx(ntest)<L | sfx(ntest)>U,1,'first');
			if ~isempty(tj)
				T(i,k) = ntest(tj);
			else
				T(i,k) = nan;
			end
		end
	end
	save('results/figure_seqnoinddist_dat.mat',...
		'T','xi','rs','nRuns','epsilon')
else
	load('results/figure_seqnoinddist_dat.mat');
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
xlim([0,2e6]);
set(gca,'xscale','log')
ylim([0 1])
grid off;
xlabel('Number of samples','FontSize',12);
ylabel('Ratio of chains stopped','FontSize',12);
set(findall(gcf,'type','text'),'FontSize',16);
set(gca,'FontSize',14)
saveas(fh,'figs/figure_seqnoinddist','eps');

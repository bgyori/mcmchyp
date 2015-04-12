% Figure for distibution of stopping times for sequential test with
% indifference region
clear variables

reload = false;
if reload
	nRuns = 1000;
	load('results/gammas_2e6_from2e5.mat');

	rs = [0.1, 0.5, 0.7, 0.8, 0.95];
	rms = min([1./(1-rs);1./rs]);
	epsilon = 0.01;
	delta = 0.05;
    
    xi = 0.3;

	T = N*ones(length(rs),nRuns);
	for k=1:nRuns
		gamma = gammas(k);
		fname = sprintf('~/data/phd/mcmchyp/odemcmc-jakstat-1-2000000-50000-1-%d.mat',k);
		fprintf('Reading %s\n',fname);
		f = load(fname);
		sfx = cumsum(f.f);
		N = length(f.f);
		
		Nfix(k) = log(1/epsilon)/(gamma*delta^2);
		M = log(2/sqrt(epsilon*xi))/(2*gamma*delta);
		for i=1:length(rs)
			r = rs(i);
			rm = rms(i);
            n0 = floor(M*rm);
            imax = floor(log(N/n0)/log(1+xi));
            ntest = floor(n0*(1+xi).^(1:imax));
            fprintf('M: %.1f, n0: %d,imax: %d\n',M,n0,imax);
            L = ntest'*r - M;
			U = ntest'*r + M;
			tj = find(sfx(ntest)<L | sfx(ntest)>U,1,'first');
			if ~isempty(tj)
				T(i,k) = ntest(tj);
			else
				T(i,k) = nan;
			end
		end
	end
	save('results/figure_seqinddist_dat.mat',...
		'T','xi','rs','nRuns','Nfix','epsilon','delta')
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
xlim([0,2e6]);
set(gca,'xscale','log');
ylim([0 1])
grid off;
xlabel('Number of samples','FontSize',12);
ylabel('Ratio of chains stopped','FontSize',12);
set(findall(gcf,'type','text'),'FontSize',16);
set(gca,'FontSize',14)
saveas(fh,'figs/figure_seqinddist','eps');

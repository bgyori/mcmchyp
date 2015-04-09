clear variables;

reload = false;
if reload
	mtrue = 0.88746;
	delta = 0.05;
	nRuns = 1000;
	load('results/gammas_2e6_from2e5.mat');

	rs = [linspace(0.01,0.70,8), ...
		  linspace(0.75,1,20)];

	epsilons = [0.01,0.02,0.05];

	errs = zeros(length(rs),length(epsilons));
	ET = 1e8*ones(length(rs),length(epsilons));
	T = nan(nRuns,length(rs),length(epsilons));
	Nfix = zeros(nRuns,length(epsilons));

	for i=1:nRuns
		fname = sprintf('results/odemcmc-jakstat-1-2000000-50000-1-%d.mat',i);
		fprintf('Reading %s\n',fname);
		load(fname);
		sfx = cumsum(f);
		N = length(f);	
		gamma = gammas(i);

		for j=1:length(rs)
			r = rs(j);
			sfxz = sfx-(1:N)'*r;
			for k=1:length(epsilons)
				epsilon = epsilons(k);
				g = sqrt((1:N)/gamma * log(1/epsilon) + 1 + 2*log(1:N))';
				tj = find(abs(sfxz)>g,1,'first');
				Nfix(i,k) = log(1/epsilon)/(gamma*delta^2);
				if ~isempty(tj)
					T(i,j,k) = tj;
					if (r < mtrue) && (sfxz(tj) < 0)
						errs(j,k) = errs(j,k) + 1;
					elseif (r > mtrue) && (sfxz(tj) > 0)
						errs(j,k) = errs(j,k) + 1;
					end
				end
			end
		end
	end

	for j=1:length(rs)
		for k=1:length(epsilons)
			if ~any(isnan(T(:,j,k)))
				ET(j,k) = mean(T(:,j,k));
			else
				ET(j,k) = nan;
			end
		end
	end
else
	load('results/figure_seqnoindeteps_dat.mat');
end

%% Plotting (e)
fh = figure; hold on;
markersize = 8;
plot(rs,ET(:,1),'x','markers',markersize);
plot(rs,ET(:,2),'o','markers',markersize);
plot(rs,ET(:,3),'^','markers',markersize);
plot(1,nan,'-');
plot(1,nan,'--');
plot(rs,ET(:,1),'-','markers',markersize);
plot(rs,ET(:,2),'-','markers',markersize);
plot(rs,ET(:,3),'-','markers',markersize);
legend({'\epsilon=0.01','\epsilon=0.02','\epsilon=0.05'});
ylim([0 1e6]);
xlabel('Probability threshold (r)','FontSize',12);
ylabel('Average stopping time (number of samples)','FontSize',12);
set(findall(gcf,'type','text'),'FontSize',16)
set(gca,'FontSize',14)
saveas(fh,'figs/figure_seqnoindeteps','eps');
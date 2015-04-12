clear variables;

reload = false;
if reload
    mtrue = 0.88746;
    epsilon = 0.01;
    nRuns = 1000;
    load('results/gammas_2e6_from2e5.mat');

    rs = [linspace(0.01,0.70,8), ...
          linspace(0.75,1,20)];
    rms = min([1./(1-rs);1./rs]);
    deltas = [0.05,0.08,0.1];

    errs = zeros(length(rs),length(deltas));
    ET = 1e8*ones(length(rs),length(deltas));
    T = nan(nRuns,length(rs),length(deltas));
    Nfix = zeros(nRuns,length(deltas));
    
    xi = 0.3;

    for i=1:nRuns
        fname = sprintf('~/data/phd/mcmchyp/odemcmc-jakstat-1-2000000-50000-1-%d.mat',i);
        fprintf('Reading %s\n',fname);
        load(fname);
        sfx = cumsum(f);
        N = length(f);	
        gamma = gammas(i);

        for j=1:length(rs)
            r = rs(j);
			rm = rms(j);
            sfxz = sfx-(1:N)'*r;
            for k=1:length(deltas)
                delta = deltas(k);
                M = log(2/sqrt(epsilon*xi))/(2*gamma*delta);
				n0 = floor(M*rm);
				imax = floor(log(N/n0)/log(1+xi));
				ntest = floor(n0*(1+xi).^(1:imax));
				tj = find(abs(sfxz(ntest))>M,1,'first');
                Nfix(i,k) = log(1/epsilon)/(gamma*delta^2);
                if ~isempty(tj)
                    T(i,j,k) = ntest(tj);
                    if (r < mtrue) && (sfxz(ntest(tj)) < 0)
                        errs(j,k) = errs(j,k) + 1;
                    elseif (r > mtrue) && (sfxz(ntest(tj)) > 0)
                        errs(j,k) = errs(j,k) + 1;
                    end
                end
            end
        end
    end

    for j=1:length(rs)
        for k=1:length(deltas)
            if ~any(isnan(T(:,j,k)))
                ET(j,k) = mean(T(:,j,k));
            else
                ET(j,k) = nan;
            end

        end
	end
	save('results/figure_seqindetdelta_dat.mat',...
		'ET','xi','rs','nRuns','Nfix','epsilon','deltas','mtrue','errs')
else
    load('results/figure_seqindetdelta_dat.mat');
end


%% Plotting
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
Nfixm = mean(Nfix);
plot(0:0.2:1,Nfixm(1)*ones(1,6),'x--','markers',markersize);
plot(0:0.2:1,Nfixm(2)*ones(1,6),'o--','markers',markersize);
plot(0:0.2:1,Nfixm(3)*ones(1,6),'^--','markers',markersize);
legend({'\delta=0.05','\delta=0.08','\delta=0.1','seq.','fixed'});
ylim([0 1e6]);
set(gca,'yscale','linear')
xlabel('Probability threshold (r)','FontSize',12);
ylabel('Average stopping time (number of samples)','FontSize',12);
set(findall(gcf,'type','text'),'FontSize',16)
set(gca,'FontSize',14)
saveas(fh, 'figs/figure_seqindetdelta','eps');

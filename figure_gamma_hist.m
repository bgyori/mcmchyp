load results/gammas_2e6_from2e5.mat;
figure;
hist(gammas,50);
xlabel('$\hat{\gamma}^*$','interpreter','latex','fontsize',16);
set(findall(gcf,'type','text'),'FontSize',16)
set(gca,'FontSize',14)
h = findobj(gca,'Type','patch');
set(h,'FaceColor',[0 0 0],'EdgeColor','k')
saveas(gcf,'figs/figure_gamma_hist.eps','eps');
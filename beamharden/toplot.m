
mincost=min(min(tosave(1:1200,2)),min(tosave(1:1200,3)));
figure; semilogy(tosave(1:1200,1),tosave(1:1200,2)-mincost); hold on; 
semilogy(tosave(1:1200,1),tosave(1:1200,3)-mincost,'r--');
legend('FISTA-ADMM','NCG-PR');
ylabel('f(\alpha)-f*'); xlabel('# of iterations');
saveas(gcf,'costVsItr.eps','psc2');

figure; semilogy(tosave(1:1200,1),tosave(1:1200,4)); hold on; 
semilogy(tosave(1:1200,1),tosave(1:1200,5),'r--');
legend('FISTA-ADMM','NCG-PR');
ylabel('Relative Square Error'); xlabel('# of iterations');
saveas(gcf,'rmseVsItr.eps','psc2');

figure; semilogy(tosave(1:1200,6),tosave(1:1200,2)-mincost); hold on; 
semilogy(tosave(1:1200,7),tosave(1:1200,3)-mincost,'r--');
legend('FISTA-ADMM','NCG-PR');
ylabel('Relative Square Error'); xlabel('time/s');
saveas(gcf,'costVsTime.eps','psc2');

figure; semilogy(tosave(1:1200,6),tosave(1:1200,4)); hold on; 
semilogy(tosave(1:1200,7),tosave(1:1200,5),'r--');
legend('FISTA-ADMM','NCG-PR');
ylabel('Relative Square Error'); xlabel('time/s');
saveas(gcf,'rmseVsTime.eps','psc2');

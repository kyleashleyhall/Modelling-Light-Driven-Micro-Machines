fileDirectory5=['..' filesep 'RadForce-Y.csv'];
File5=csvread(fileDirectory5,1,0);
figure
hold on;
quiver3(File5(:,1),File5(:,2),File5(:,3),File5(:,5),File5(:,6),File5(:,7),'LineWidth',2);
set(gca,'TickLabelInterpreter','latex')
xlabel('x \(\left(\mu m\right)\)','Interpreter','latex')
ylabel('y \(\left(\mu m\right)\)','Interpreter','latex')
zlabel('z \(\left(\mu m\right)\)','Interpreter','latex')
mylegend=legend({'$Force \left(10^{-13} N\right)$'});
set(mylegend,'Interpreter','latex')
hold off;
view(90,0)
clear
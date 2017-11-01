fileDirectory5=['..' filesep '2_0i' filesep 'RadForce-Y.csv'];
File5=csvread(fileDirectory5,1,0);
fileDirectory6=['..' filesep '2_0.5i' filesep 'RadForce-Y.csv'];
File6=csvread(fileDirectory6,1,0);
fileDirectory7=['..' filesep '2_1i' filesep 'RadForce-Y.csv'];
File7=csvread(fileDirectory7,1,0);
figure
hold on;
k=boundary(File5(:,1),File5(:,2),File5(:,3));
trisurf(k,File5(:,1),File5(:,2),File5(:,3),'Facecolor','blue','FaceAlpha',1)
set(gca,'TickLabelInterpreter','latex')
xlabel('x \(\left(\mu m\right)\)','Interpreter','latex')
ylabel('y \(\left(\mu m\right)\)','Interpreter','latex')
zlabel('z \(\left(\mu m\right)\)','Interpreter','latex')
hold off;
view(45,45)
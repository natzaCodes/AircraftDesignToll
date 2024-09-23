% Authors:
% Albert Canovas Cots
% Natalia Zalewska
function PlotPosicions(geo,aerodata,lh,weights,CG,mtom)

wing = [geo.wing.xle, geo.wing.xle + geo.wing.cr];
zer = [0,0]
tail =[geo.wing.xle + lh, geo.wing.xle + lh + geo.htail.cr];
fact_size = 200;
figure
plot(wing,zer,'linewidth',3)
hold on;
grid on;
plot(tail,zer,'linewidth',3)
scatter(CG,0, sum(cell2mat(weights(:,3)))*fact_size,'filled')
scatter(geo.wing.xle +aerodata.xnp,0,500,'d','filled')
scatter(cell2mat(weights(:,2)),zeros(size( cell2mat(weights(:,2)))), cell2mat(weights(:,3))*2*fact_size,'filled')
xlabel("XPosition (m)")
ylabel("Zposition (m)")
title("Mass distribution")
legend("Wing","Tail","CG","NP","Weights")
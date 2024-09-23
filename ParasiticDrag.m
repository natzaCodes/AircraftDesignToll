% Authors:
% Albert Canovas Cots
% Natalia Zalewska
function parDrag = ParasiticDrag(V,geo)
h = 0; % flight altitude set to 0
%% Wing
[parDrag.Df_w] = SkinFrictionDrag(h,V,geo.wing.mac,0,geo.wing.sweep,geo.wing.airfoil.thickness,geo.wing.airfoil.thickpos,geo.wing.sw,'wing');

%% Fuselage
d_fuselage = geo.fuselage.diam;
l_fuselage = geo.fuselage.length;
exc = sqrt((l_fuselage/2)^2-(d_fuselage/2)^2)/(l_fuselage/2);
S_fuselage = 2 * pi * ((d_fuselage/2)^2 + d_fuselage*l_fuselage/(4*exc)*asin(exc));
f_fuselage = l_fuselage/d_fuselage;
[parDrag.Df_f] = SkinFrictionDrag(h,V,l_fuselage,f_fuselage,0,0,0,S_fuselage,'fuselage');


%% Tail H&V
[parDrag.Df_Ht] = SkinFrictionDrag(h,V,geo.htail.mac,0,geo.htail.sweep,geo.htail.airfoil.thickness,geo.htail.airfoil.thickpos,geo.htail.sw,'tail');

[parDrag.Df_Vt] = SkinFrictionDrag(h,V,geo.vtail.mac,0,geo.htail.sweep,geo.vtail.airfoil.thickness,geo.vtail.airfoil.thickpos,geo.htail.sw,'tail');

%% Miscellaneous
dk = 0.2;

parDrag.D = parDrag.Df_w + parDrag.Df_f + parDrag.Df_Ht + parDrag.Df_Vt;
parDrag.Dmisc = parDrag.D * dk;
parDrag.D = parDrag.D + parDrag.Dmisc;

[~,~,rho,~,~,~,~] = ISAfunction(h);
qs = 0.5 * rho .* V .* V * geo.wing.sw;
parDrag.CD0 = parDrag.D ./ qs;
parDrag.CDw = parDrag.Df_w ./ qs; 
parDrag.CDf = parDrag.Df_f ./ qs; 
parDrag.CDh = parDrag.Df_Ht ./ qs;
parDrag.CDv = parDrag.Df_Vt ./ qs;
parDrag.CDmisc = parDrag.Dmisc ./ qs;
end


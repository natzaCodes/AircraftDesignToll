%%%%%%%%%%%%% Aircraft Design Tool %%%%%%%%%%%%%%%%%%%
% Albert Canovas Cots
% Natalia Zalewska
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all;
close all;
rho = 1.3;
%% Wing Geometry
% ----------------------------------------------------------------------- %

% Model 14 Naca TR 2445
% geo.wing.sweep = deg2rad(0);            % Leading edge sweep angle (rad)
% geo.wing.phi = deg2rad(0);  % Dihedral angle (rad)
% geo.wing.b2   = 1.5;            % SEMI-span (m)
% geo.wing.cr  = 0.45;            % Root chord (m)
% geo.wing.ct  = 0.35;%0.254;  % Tip chord (m)
% geo.wing.tr = geo.wing.ct / geo.wing.cr;
% geo.wing.mac = 2/3 * geo.wing.cr * (geo.wing.tr^2 + geo.wing.tr + 1) / (geo.wing.tr + 1);
% geo.wing.sw   = (geo.wing.ct + geo.wing.cr)*geo.wing.b2;            % Wing surface (m)
% geo.wing.ar  = ((2*geo.wing.b2)^2)/geo.wing.sw;            % Wing Aspect Ratio (-)
% geo.wing.airfoil = ReadAirfoil("S1223",12);
% geo.wing.xle = 0.2;
% %%%
% % Tail coefficients
% tailcoef.vv = 0.03; % Vertical tail volume coefficient
% tailcoef.vh = 0.6; % Vertical tail volume coefficient
% tailcoef.lh = 1.25; % Lenght of the fuselage
% 
% %% Tail sizing
% % ----------------------------------------------------------------------- %
% geo = TailSizing(geo,tailcoef);
% 
% %% Fuselage sizing
% % ----------------------------------------------------------------------- %
% geo.fuselage.diam = 0.1;
% geo.fuselage.length = tailcoef.lh +geo.wing.xle;
% 
% %% Aero derivatives
% % ----------------------------------------------------------------------- %
% aeroCoefs.wing = PreVLM(geo.wing);
% aeroCoefs.htail = PreVLM(geo.htail);
% 
% %% Parasitic drag estimation
% 
% parDrag = ParasiticDrag(25,geo)
% 
% aerodata.cl0 = aeroCoefs.wing.cl0 + aeroCoefs.htail.cl0 * geo.htail.sw/geo.wing.sw;
% aerodata.clalpha = aeroCoefs.wing.clalpha + aeroCoefs.htail.clalpha * geo.htail.sw/geo.wing.sw;
% aerodata.cd0 = aeroCoefs.wing.cdi0 + aeroCoefs.htail.cdi0 * geo.htail.sw/geo.wing.sw;
% aerodata.klcdi = aeroCoefs.wing.klcdi + aeroCoefs.htail.klcdi * geo.htail.sw/geo.wing.sw;
% aerodata.kcdi = aeroCoefs.wing.kcdi + aeroCoefs.htail.kcdi * geo.htail.sw/geo.wing.sw;
% aerodata.cd0 = aerodata.cd0 + parDrag.CD0
% aerodata.cmacw = aeroCoefs.wing.cmac;
% 
% %% Neutral point estimation
% 
% 
% 
% 
% %% TakeOff analysis
% 
% toffdata.tolength = 60;
% toffdata.mu = 0.06;
% toffdata.aoa = 0;
% 
% powerdata.pot = 700;
% powerdata.etapto = 0.6;
% powerdata.etap = 0.8;
% 
% mtom = TOSolver(geo,aerodata,toffdata,powerdata)
% mtom = 13;
% vstall = sqrt(2*mtom*9.81/(rho*geo.wing.sw*(2.1/(1.2)^2)))
% %% Bendingmoment analysis
% aerodata.bm = aeroCoefs.wing.bmadim * mtom * 9.81;
% aerodata.bmy = aeroCoefs.wing.bmadimy * mtom * 9.81;
load("aerodat.mat")
%% Weights approximation
% ----------------------------------------------------------------------- %

% payload 10x10x20 cm


construction.wing.rib.spacing = 0.15;
construction.wing.rib.last = 0.2; %location for last rib brfore root/fuselage, first rib is at the tip
construction.wing.rib.location = -(1)*(-geo.wing.b2:construction.wing.rib.spacing:-construction.wing.rib.last); %from tip to root
construction.wing.rib.material.thickness = 0.02;
construction.wing.rib.NACA2412.known.surface = 0.121; %m2 from chat gpt
construction.wing.rib.NACA2412.known.chord = 1; %m from chat gpt
construction.wing.rib.NACA2412.cg =  0.3994; % *chord
construction.wing.spar.location = 0.25; % spar location in % of chord
construction.wing.spar.location_glob = geo.wing.xle + geo.wing.cr * construction.wing.spar.location; % spar location in X, 0=front end of fuselage
construction.wing.D_box.thickness = 0.002; % thickness of D-wall of D-box
construction.wing.joint.thickness = 0.002; %m
construction.wing.joint.length = 0.02; %m

construction.htail.rib.spacing = 0.15;
construction.htail.rib.last = 0.01; %location for last rib brfore root/fuselage, first rib is at the tip
construction.htail.rib.location = -(1)*(-geo.htail.b2:construction.htail.rib.spacing:-construction.htail.rib.last); %from tip to root
construction.htail.rib.material.thickness = 0.02;
construction.htail.rib.NACA2412.known.surface = 0.121; %m2 from chat gpt
construction.htail.rib.NACA2412.known.chord = 1; %m from chat gpt
construction.htail.spar.location = 0.25; % spar location in % of chord
construction.htail.D_box.thickness = 0.002; % thickness of D-wall of D-box

construction.vtail.rib.spacing = 0.15;
construction.vtail.rib.last = 0.01; %location for last rib brfore root/fuselage, first rib is at the tip
construction.vtail.rib.location = -(1)*(-geo.htail.b2:construction.htail.rib.spacing:-construction.htail.rib.last); %from tip to root
construction.vtail.rib.material.thickness = 0.02;
construction.vtail.rib.NACA2412.known.surface = 0.121; %m2 from chat gpt
construction.vtail.rib.NACA2412.known.chord = 1; %m from chat gpt
construction.vtail.spar.location = 0.25; % spar location in % of chord
construction.vtail.D_box.thickness = 0.002; % thickness of D-wall of D-box


construction.fuselage.frame.xthickness = 0.02;
construction.fuselage.frame.zthickness = 0.01;
construction.fuselage.tail.thickness = 0.001; %m
construction.fuselage.tail.Diameter = 0.02; %m
construction.fuselage.stringer.ythickness = 0.02; %m
construction.fuselage.stringer.zthickness = 0.02; %m
construction.fuselage.payload.mass = [0.0, 2]; % in reference to wing spar
construction.fuselage.bottle.mass = [0.0, 2];  % in reference to wing spar


construction.material.balsa.sigma_max = 10e6; %Pa
construction.material.balsa.density = 150; %kg/m3
construction.material.CFRP.density = 2000; %kg/m3
construction.material.CFRP.sigma_max = 5e6; %Pa

construction.fuselage.frame.xlocation_diameter = [  0,  0.05;
                                                  construction.wing.spar.location_glob-0.05,   0.15;
                                                  0.55,  0.15;
                                                  0.7, 0.04];

[construction, geo] = Masses(geo, aerodata, construction);
[allData, allData_sorted, CG] = xBalance(geo,construction,tailcoef);
disp(mtom)
geo.mass = 2.67 + geo.wing.mass + geo.fuselage.mass + geo.htail.mass + geo.vtail.mass
aerodata.xnp = Neutralpoint(geo,aerodata,tailcoef.lh,aeroCoefs)
PlotPosicions(geo,aerodata,tailcoef.lh,allData_sorted,CG,mtom)

disp(geo.mass)
disp(mtom)

v = linspace(6,44,51);
clv = mtom*9.81*2./(rho*v.^2*geo.wing.sw);
cdv = aerodata.cd0 + aerodata.klcdi*clv + aerodata.kcdi*clv.^2;
dv = 0.5 * rho*v.^2*geo.wing.sw.*cdv;
t = powerdata.pot*powerdata.etap./v;

%%% Bending moment plots
figure
grid on;
hold on;
plot(aeroCoefs.wing.mesh.coor.cpbb(1,2:end,2),aerodata.bmy','--','linewidth',1.5);
plot(aeroCoefs.wing.mesh.coor.cpbb(1,2:end,2),aerodata.bmy'*4,'-','linewidth',1.5);
plot(aeroCoefs.wing.mesh.coor.cpbb(1,2:end,2),aerodata.bmy'*(-2),':','linewidth',1.5);
title("Bending moment distirbution")
xlabel("y (m)")
ylabel("Bending moment (Nm)")
legend("n=1","n=4","n=-2")

aerodata.sheary = diff(aerodata.bmy')./diff(aeroCoefs.wing.mesh.coor.cpbb(1,2:end,2))
figure
grid on;
hold on;
plot(aeroCoefs.wing.mesh.coor.cpbb(1,2:end-1,2),aerodata.sheary','--','linewidth',1.5);
plot(aeroCoefs.wing.mesh.coor.cpbb(1,2:end-1,2),aerodata.sheary'*4,'-','linewidth',1.5);
plot(aeroCoefs.wing.mesh.coor.cpbb(1,2:end-1,2),aerodata.sheary'*(-2),':','linewidth',1.5);
title("Shear stress distirbution")
xlabel("y (m)")
ylabel("Shear stress (N)")
legend("n=1","n=4","n=-2")

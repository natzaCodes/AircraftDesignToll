% Authors:
% Albert Canovas Cots
% Natalia Zalewska
function aeroCoefs = PreVLM(geo)
	
	%% Mesh Configuration
	% ----------------------------------------------------------------------- %
	mesh.npx = 15; % Number of panels in the streamwise direction
	mesh.npy = 30; % Number of panels in the SEMI-spanwise direction
	mesh.npt = mesh.npx*mesh.npy;
	
	pointsAlpha = [0,3,6]
	%% First analysis point 
	flow.alpha = deg2rad(pointsAlpha(1)); % Angle of Attack (rad)
	flow.tas  = 53;       % Wind Speed (m/s)
	flow.rho   = 1.225;    % Air Density (kg/m3)
	[p1,mesh] = VLM(geo,flow,mesh)
	
	%% Second analysis point 
	flow.alpha = deg2rad(pointsAlpha(2)); % Angle of Attack (rad)
	[p2,~] = VLM(geo,flow,mesh)

	%% Third analysis point 
	flow.alpha = deg2rad(pointsAlpha(3)); % Angle of Attack (rad)
	[p3,~] = VLM(geo,flow,mesh)
	
	
	%% Derivatives
	aeroCoefs.CLs = [p1.CL,p2.CL,p3.CL];
	aeroCoefs.CDs = [p1.CD,p2.CD,p3.CD];
	aeroCoefs.CMs = [p1.CMLe,p2.CMLe,p3.CMLe];
	auxcl = polyfit(pointsAlpha,aeroCoefs.CLs,1);
	aeroCoefs.clalpha = auxcl(1);
	aeroCoefs.cl0 = auxcl(2);
	auxcdi = polyfit(aeroCoefs.CLs,aeroCoefs.CDs,2);
	aeroCoefs.kcdi = auxcdi(1);
	aeroCoefs.klcdi = auxcdi(2);
	aeroCoefs.cdi0 = auxcdi(3);
	aeroCoefs.bmadim = p1.BM
	aeroCoefs.bmadimy = p1.bmy;
	aeroCoefs.xac = (aeroCoefs.CMs(3)-aeroCoefs.CMs(1))/(aeroCoefs.CLs(3)-aeroCoefs.CLs(1))*geo.mac;
	aeroCoefs.cmac = (aeroCoefs.CMs(3)-aeroCoefs.CLs(3)*aeroCoefs.xac)
	aeroCoefs.mesh = mesh;
end
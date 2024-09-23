% Authors:
% Albert Canovas Cots
% Natalia Zalewska
function geo = TailSizing(geo,tailcoef)
	%% Horizontal tail
	
	geo.htail.sweep = geo.wing.sweep;
	geo.htail.phi = geo.wing.phi;
	geo.htail.ar = 1/2 * geo.wing.ar;
	geo.htail.tr = geo.wing.tr;
	geo.htail.sw = (tailcoef.vh * geo.wing.sw * geo.wing.mac) / (tailcoef.lh) ;
	geo.htail.b2 = sqrt(geo.htail.ar * geo.htail.sw) / 2;
	geo.htail.cr = geo.htail.sw / (geo.htail.b2 * (1+geo.htail.tr));
	geo.htail.ct = geo.htail.cr * geo.htail.tr;
	geo.htail.mac = 2/3 * geo.htail.cr * (geo.htail.tr^2 + geo.htail.tr + 1) / (geo.htail.tr + 1);
	geo.htail.airfoil =  ReadAirfoil("NACA0012",12);
	%% Vertical tail
	
	geo.vtail.sweep = geo.wing.sweep;
	geo.vtail.phi = deg2rad(90);
	geo.vtail.ar = 2;
	geo.vtail.tr = geo.wing.tr;
	geo.vtail.sw = (tailcoef.vv * geo.wing.sw * 2* geo.wing.b2) / (tailcoef.lh);
	geo.vtail.b2 = sqrt(geo.vtail.ar * geo.vtail.sw);
	geo.vtail.cr = 2 * geo.vtail.sw / (geo.vtail.b2 * (1+geo.vtail.tr));
	geo.vtail.ct = geo.vtail.cr * geo.vtail.tr;
	geo.vtail.mac = 2/3 * geo.vtail.cr * (geo.vtail.tr^2 + geo.vtail.tr + 1) / (geo.vtail.tr + 1);
	geo.vtail.airfoil =  ReadAirfoil("NACA0012",12);

end
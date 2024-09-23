% Authors:
% Albert Canovas Cots
% Natalia Zalewska
function np = Neutralpoint(geo,aeroData,lh,aeroCoefs)
	alphadummy = 5;
	
	clw1 = aeroCoefs.wing.cl0;
	clt1 = aeroCoefs.htail.cl0 * geo.htail.sw/geo.wing.sw;
	
	clw2 = aeroCoefs.wing.cl0 + alphadummy * aeroCoefs.wing.clalpha;
	clt2 = (aeroCoefs.htail.cl0 + alphadummy * aeroCoefs.htail.clalpha) * geo.htail.sw/geo.wing.sw;

	np = (clw2 *aeroCoefs.wing.xac + clt2 * lh - clw1*aeroCoefs.wing.xac - clt1 * lh) / (clw2 + clt2 - clw1 -clt1)


end
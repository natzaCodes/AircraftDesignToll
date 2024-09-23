% Authors:
% Albert Canovas Cots
% Natalia Zalewska
function airfoil = ReadAirfoil(filename,thickn)
	coor = readmatrix(strcat("airfoils/",filename));
	
	airfoil.x = coor(:,1);
	airfoil.z = coor(:,2);
	airfoil.thickness = thickn/100;
	airfoil.thickpos = 0.4;


end

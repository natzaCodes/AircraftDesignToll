% Authors:
% Albert Canovas Cots
% Natalia Zalewska
function [t,p,rho,vmach,mu,liuid1,liuid2] = ISAfunction(hg)
	% @brief calculates pressure, density, dynamic viscosity,
	% sound speed as a function of geometric altitude (up to 86 km)
	%@param h geometric altitude [m]
	%@return t Temperature [k]
	%@return p Pressure [Pa]
	%@return rho Density [kg(m^3]
	%@return mu Dynamic viscosity [kg/(ms)]
	%@return vmach Sound speed [m/s]
	%@return liuid1 LiU ID of author 1
	%@return liuid2 LiU ID of author 2
	
	%% Variable defiinition
	t = zeros(size(hg));%temperature
	p = zeros(size(hg));%pressure
	rho = zeros(size(hg));%density
	mu = zeros(size(hg));%dynamic viscosity
	vmach = zeros(size(hg));%sound speed
	radius = 6.356766e6;
	h = radius./(radius+hg).*hg;
	lambda = [-6.5,0,1,2.8,0,-2.8,-2]* 1e-3;%Temperature gradient at each section
	hlim = [0,11,20,32,47,51,71,84.852]*1e3;%Altitude where layers start/end
	tbase = zeros(size(lambda));%bottom temperature at each section
	pbase = zeros(size(lambda));%bottom pressure at each section
	
	tsea = 288.15;%Temperature at sea level
	psea = 101325;%Pressure at sea level
	g = 9.80665;%% TODO: revise expression
	[R,M] = averageM();
	gamma = 1.4;% adiabatic constant
	beta = 1.458e-6;%Coefficient 1 for dynamic viscosity
	k = 110.4;%Coefficient 2 for dynamic viscosity
	
	%% Computation of base properties of each atmosphere section
	tbase(1) = tsea;
	pbase(1) = psea;
	for i=2:size(lambda,2)
		tbase(i) = tbase(i-1) + lambda(i-1) * (hlim(i)-hlim(i-1));
		if(abs(lambda(i-1))<1e-8)
			pbase(i) = pbase(i-1) * exp(-g * (hlim(i) - hlim(i-1)) / (R * tbase(i)));% TODO: revise expression
		else
			pbase(i) = pbase(i-1) * (tbase(i) / tbase(i-1))^(-g / (R * lambda(i-1)));% TODO: revise expression
		end
	end
	
	%% LiU ID
	liuid1 = "albca693";
	liuid2 = "natza019";
	
	%% Computation of vector
	for i=1:size(h,2)
		if (h(i)<0) 
			throw(MException("myComponent:inputError","Input altitude must be above 0"));
		elseif (h(i)>90e3)
			throw(MException("myComponent:inputError","Input altitude must be below 90 km"));
		end
		% Atmoshpere ayer identification
		if(h(i)<hlim(2))
			idxb = 1;
		elseif(h(i)<hlim(3))
			idxb = 2;
		elseif(h(i)<hlim(4))
			idxb = 3;
		elseif(h(i)<hlim(5))
			idxb = 4;
		elseif(h(i)<hlim(6))
			idxb = 5;
		elseif(h(i)<hlim(7))
			idxb = 6;
		end
		
		%Temperature calculation
		t(i) = tbase(idxb) + lambda(idxb) * (h(i) - hlim(idxb));
		
		%Pressure calculation
		if(abs(lambda(idxb))<1e-8)
			p(i) = pbase(idxb) * exp(-g * (h(i) - hlim(idxb)) / (R * tbase(idxb)));% TODO: revise expression
		else
			p(i) = pbase(idxb) * (t(i) / tbase(idxb))^(-g / (R * lambda(idxb)));% TODO: revise expression
		end
		
		%Sound spedd calculation
		vmach(i) = sqrt(gamma * R * t(i));
		
		%Density calculation
		rho(i) = p(i) / (R * t(i));%% TODO: revise expression

		% Dynamic vsicosity calculation
		mu(i) = beta * t(i)^(3/2) / (t(i) + k);

	end
end
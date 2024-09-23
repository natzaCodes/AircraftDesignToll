% Authors:
% Albert Canovas Cots
% Natalia Zalewska
function mtom = TOSolver(geo,aerodata,toffdata,powerdata,rho)
	m0 = 10;
	g = 9.81;
	clmax = 2.1;
	clrot = clmax/1.2^2;
	
	clto = aerodata.cl0 + aerodata.clalpha*toffdata.aoa;
	cdto = aerodata.cd0 + aerodata.kcdi*clto^2 + aerodata.klcdi*clto;
	cdg = cdto-clto*toffdata.mu;
	iter_max = 100;
	relax = 0.5;
	v_to = sqrt(2*m0*g/(rho*geo.wing.sw*clrot));
	aux = powerdata.etapto*powerdata.pot / (g * v_to);
	aux2 = 0.6*rho*g*cdg*toffdata.tolength*geo.wing.sw/(m0*g);
	for iter=1:iter_max

		m = aux*(1-exp(aux2))/(toffdata.mu-(toffdata.mu+cdg/clrot)*exp(aux2));
		m0 = m*relax + m0*(1-relax);
		v_to = sqrt(2*m0*g/(rho*geo.wing.sw*clmax));
		aux = powerdata.etapto*powerdata.pot / (g * v_to);
		aux2 = 0.6*rho*g*cdg*toffdata.tolength*geo.wing.sw/(m0*g);
		res =(m0-aux*(1-exp(aux2))/(toffdata.mu-(toffdata.mu+cdg/clrot)*exp(aux2)))/m0;
		if(abs(res)<1e-7)
			break
		end
	end
	mtom = m;
end
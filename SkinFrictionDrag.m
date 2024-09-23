% Authors:
% Albert Canovas Cots
% Natalia Zalewska
function [Df] = SkinFrictionDrag(h,V,l,f,Lambda,TC,xc,S,name)
[~,~,rho,vmach,mu,~,~] = ISAfunction(h);                                    
M = V/vmach;                                                                % Mach number

Re = rho*V*l/mu;                                                            % Reynolds number 
Cf = 0.455./(log10(Re).^2.58 .* (1+0.144*M.^2).^0.65);                             %skin friction coefficient (turbulent)


%% Form factors and Interference Factor

if strcmp(name,'wing')  
    FF = (1+ 0.6/xc *TC + 100*TC^4) * (1.34*M.^0.18 * cos(Lambda)^0.28);
    Q = 1.0;
elseif strcmp(name,'tail')
    FF = (1+ 0.6/xc *TC + 100*TC^4) * (1.34*M.^0.18 * cos(Lambda)^0.28);
    Q = 1.06;                                                                  
elseif strcmp(name,'fuselage')
    FF = 1 + 5/f^1.5 + f/400;
    Q = 1.0;             
elseif strcmp(name,'nacelle')
    FF = 1 + 0.35/f;  
    Q = 1.3;
end

Cfe = Cf .* FF * Q;
Df = 0.5*rho.*V.^2*S.*Cfe;
end
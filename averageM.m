% Authors:
% Albert Canovas Cots
% Natalia Zalewska
function [R,M] = averageM();
%% Calculates the average molecular weight
% @return R constant R [Nm/(kg K)]
% @return M molecular weights [Kg/kmol]

%molecular weights
Mi=[28.0134, 31.9988, 39.948, 44.00995, 20.183, 4.0026, 83.80, 131.30, 16.04303, 2.01594];
% fractions
Fi =[0.78084, 0.209476, 0.00934, 0.000314, 0.00001818, 0.00000524, 0.00000114, 0.000000087, 0.000002, 0.0000005];

M = sum(Mi.*Fi);
constant = 8.31432e3;%gas constant
R = constant/M;
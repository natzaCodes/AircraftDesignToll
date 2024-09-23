% Authors:
% Albert Canovas Cots
% Natalia Zalewska
function [construction,geo] = Masses(geo, aerodata, construction)
% MASSES calculates masses of main airplane parts and gives their local location
% .mass = [local_location; mass]
%   Focus is on the parts that are along x-axis (like fuselage) bc we need
%   masses distribution for balance calculations

%% ------------ Wing 
% wing is considered as D-box with spar, ribs + joints

% --------------Spar
rib_chord = geo.wing.cr + ((construction.wing.rib.location)/geo.wing.b2 * (geo.wing.ct-geo.wing.cr)); %chord length for ribs locations
construction.wing.rib.chord = rib_chord;
y = construction.wing.spar.location * rib_chord; % chord locations for each rib to calculate thickness in these places


% airfoil thickness for y chord location, which is height of the spar
construction.wing.rib.thickness = geo.wing.airfoil.thickness*rib_chord; 


construction.wing.spar.thickness = (6*aerodata.bm/construction.material.balsa.sigma_max)^(1/3);
construction.wing.spar.height = construction.wing.rib.thickness;
construction.wing.spar.length = geo.wing.b2; %on one side

n = size(y);
construction.wing.spar.volume = 0.5*(construction.wing.spar.height(1) + construction.wing.spar.height(n(2))) * construction.wing.spar.length * construction.wing.spar.thickness; %one side
construction.wing.spar.volume = construction.wing.spar.volume * 2 ; % full spar
construction.wing.spar.mass(1,1) = 0;  % local location
construction.wing.spar.mass(1,2) = construction.wing.spar.volume * construction.material.balsa.density;

%------------ D-box: (big elipse-smaller elipse)/2 - cross-section; 
% big pyramid-small pyramid  - this is volume
P_big_elipse = pi/2 * y .* construction.wing.spar.height/2; %surface of half big elipse
P_small_elipse = pi/2 * (y-construction.wing.D_box.thickness ) .* (construction.wing.spar.height/2 - construction.wing.D_box.thickness);%surface of half small elipse

pyramid_big = P_big_elipse(1) * construction.wing.rib.spacing/3;
pyramid_small = P_small_elipse(1) * (construction.wing.rib.spacing-construction.wing.D_box.thickness)/3;


for i= 1:size(y')-1
    % volumes between each two neighbour ribs
    D_box_volume0 = pyramid_big - pyramid_small;
    pyramid_big = P_big_elipse(i+1) * (construction.wing.rib.spacing + construction.wing.rib.location(1) - construction.wing.rib.location(i+1))/3;
    pyramid_small = P_small_elipse(i+1) * (construction.wing.rib.spacing + construction.wing.rib.location(1) - construction.wing.rib.location(i+1)-construction.wing.D_box.thickness)/3;
 
    %volumes for half of wing
    construction.wing.D_box.volume(i) = (pyramid_big-pyramid_small) - D_box_volume0; 

end

local_mass_layout = y + 3*y/8 ;
N = size(y);
mass_location = -sum(local_mass_layout)/N(2) ; %location in reference to spar

construction.wing.D_box.volume_total = sum(construction.wing.D_box.volume) * 2 ;
construction.wing.D_box.mass(1,1) = mass_location;
construction.wing.D_box.mass(1,2) = construction.wing.D_box.volume_total * construction.material.balsa.density;

%-------------------- Ribs: 
% surface of the known rib (with known chord size) * scaling factor = our rib
scaling_rib = rib_chord/construction.wing.rib.NACA2412.known.chord;
construction.wing.rib.surface = scaling_rib .*scaling_rib *construction.wing.rib.NACA2412.known.surface;

local_mass_layout = rib_chord*construction.wing.rib.NACA2412.cg - y;
N = size(y);
mass_location = sum(local_mass_layout)/N(2); %location in reference to spar

ribs_volume = sum(construction.wing.rib.surface)*construction.wing.rib.material.thickness;
construction.wing.rib.mass(1,1) = mass_location;
construction.wing.rib.mass(1,2) = ribs_volume*construction.material.balsa.density *2;

%-------------------- Joints:
% calc dimeter based on material strenth and bending moment
syms x positive; % outher Diameter of "empty" joint
eqn = (128/pi*(x^4 - (x-construction.wing.joint.thickness)^4))-aerodata.bm*x/construction.material.CFRP.sigma_max == 0;
Diameter = solve(eqn,x,'Real',true, 'MaxDegree', 3);
construction.wing.joint.D=eval(Diameter);
construction.wing.joint.mass(1,1) = 0;
construction.wing.joint.mass(1,2) = pi* (construction.wing.joint.D/2)^2 * construction.wing.joint.length *construction.material.CFRP.density *2;

geo.wing.mass = construction.wing.spar.mass(1,2) + construction.wing.D_box.mass(1,2) + construction.wing.rib.mass(1,2) + construction.wing.joint.mass(1,2);

%% Vertical stabilizer
% stabilizer is considered as D-box with ribs

% --------------Spar
rib_chord = geo.vtail.cr + ((construction.vtail.rib.location)/geo.vtail.b2 * (geo.vtail.ct-geo.vtail.cr)); %chord length for ribs locations
construction.vtail.rib.chord = rib_chord;
y = construction.vtail.spar.location * rib_chord; % chord locations for each rib to calculate thickness in these places


% airfoil thickness for y chord location, which is height of the spar   
construction.vtail.rib.thickness = geo.vtail.airfoil.thickness*rib_chord; 


construction.vtail.spar.thickness = (6*aerodata.bm/construction.material.balsa.sigma_max)^(1/3);
construction.vtail.spar.height = construction.vtail.rib.thickness;
construction.vtail.spar.length = geo.vtail.b2; %on one side

n = size(y);
construction.vtail.spar.volume = 0.5*(construction.vtail.spar.height(1) + construction.vtail.spar.height(n(2))) * construction.vtail.spar.length * construction.vtail.spar.thickness; %one side
construction.vtail.spar.volume = construction.vtail.spar.volume ; % full spar
construction.vtail.spar.mass(1,1) = 0; % in reference to vtail spar
construction.vtail.spar.mass(1,2) = construction.vtail.spar.volume * construction.material.balsa.density;

%------------ D-box: (big elipse-smaller elipse)/2 - cross-section; 
% big pyramid-small pyramid  - this is volume
P_big_elipse = pi/2 * y .* construction.vtail.spar.height/2; %surface of half big elipse
P_small_elipse = pi/2 * (y-construction.vtail.D_box.thickness ) .* (construction.vtail.spar.height/2 - construction.vtail.D_box.thickness);%surface of half small elipse

pyramid_big = P_big_elipse(1) * construction.vtail.rib.spacing/3;
pyramid_small = P_small_elipse(1) * (construction.vtail.rib.spacing-construction.vtail.D_box.thickness)/3;


for i= 1:size(y')-1
    % volumes between each two neighbour ribs
    D_box_volume0 = pyramid_big - pyramid_small;
    pyramid_big = P_big_elipse(i+1) * (construction.vtail.rib.spacing + construction.vtail.rib.location(1) - construction.vtail.rib.location(i+1))/3;
    pyramid_small = P_small_elipse(i+1) * (construction.vtail.rib.spacing + construction.vtail.rib.location(1) - construction.vtail.rib.location(i+1)-construction.vtail.D_box.thickness)/3;
 
    %volumes for half of wing
    construction.vtail.D_box.volume(i) = (pyramid_big-pyramid_small) - D_box_volume0; 
end
local_mass_layout = y + 3*y/8 ;
N = size(y);
mass_location = -sum(local_mass_layout)/N(2) ; %location in reference to spar

construction.vtail.D_box.volume_total = sum(construction.vtail.D_box.volume) ;
construction.vtail.D_box.mass(1,1) = mass_location;
construction.vtail.D_box.mass(1,2) = construction.vtail.D_box.volume_total * construction.material.balsa.density;

%-------------------- Ribs: 
% surface of the known rib (with known chord size) * scaling factor = our rib
scaling_rib = rib_chord/construction.vtail.rib.NACA2412.known.chord;
construction.vtail.rib.surface = scaling_rib .*scaling_rib *construction.vtail.rib.NACA2412.known.surface;

local_mass_layout = rib_chord*construction.wing.rib.NACA2412.cg - y;
N = size(y);
mass_location = sum(local_mass_layout)/N(2); %location in reference to spar
ribs_volume = sum(construction.vtail.rib.surface)*construction.vtail.rib.material.thickness;
construction.vtail.rib.mass(1,1) = mass_location;
construction.vtail.rib.mass(1,2) = ribs_volume*construction.material.balsa.density;


geo.vtail.mass = construction.vtail.spar.mass(1,2) + construction.vtail.D_box.mass(1,2) + construction.vtail.rib.mass(1,2);


%% Horizontal stabilizer
% stabilizer is considered as D-box with ribs

% --------------Spar
rib_chord = geo.htail.cr + ((construction.htail.rib.location)/geo.htail.b2 * (geo.htail.ct-geo.htail.cr)); %chord length for ribs locations
construction.htail.rib.chord = rib_chord;
y = construction.htail.spar.location * rib_chord; % chord locations for each rib to calculate thickness in these places

construction.htail.rib.thickness = geo.htail.airfoil.thickness*rib_chord; 

construction.htail.spar.thickness = (6*aerodata.bm/construction.material.balsa.sigma_max)^(1/3);
construction.htail.spar.height = construction.htail.rib.thickness;
construction.htail.spar.length = geo.htail.b2; %on one side

n = size(y);
construction.htail.spar.volume = 0.5*(construction.htail.spar.height(1) + construction.htail.spar.height(n(2))) * construction.htail.spar.length * construction.htail.spar.thickness; %one side
construction.htail.spar.volume = construction.htail.spar.volume * 2 ; % full spar
construction.htail.spar.mass(1,1) = 0; % in reference to vtail spar
construction.htail.spar.mass(1,2) = construction.htail.spar.volume * construction.material.balsa.density;

%------------ D-box: (big elipse-smaller elipse)/2 - cross-section; 
% big pyramid-small pyramid  - this is volume
P_big_elipse = pi/2 * y .* construction.htail.spar.height/2; %surface of half big elipse
P_small_elipse = pi/2 * (y-construction.htail.D_box.thickness ) .* (construction.htail.spar.height/2 - construction.htail.D_box.thickness);%surface of half small elipse

pyramid_big = P_big_elipse(1) * construction.htail.rib.spacing/3;
pyramid_small = P_small_elipse(1) * (construction.htail.rib.spacing-construction.htail.D_box.thickness)/3;


for i= 1:size(y')-1
    % volumes between each two neighbour ribs
    D_box_volume0 = pyramid_big - pyramid_small;
    pyramid_big = P_big_elipse(i+1) * (construction.htail.rib.spacing + construction.htail.rib.location(1) - construction.htail.rib.location(i+1))/3;
    pyramid_small = P_small_elipse(i+1) * (construction.htail.rib.spacing + construction.htail.rib.location(1) - construction.htail.rib.location(i+1)-construction.htail.D_box.thickness)/3;
 
    %volumes for half of wing
    construction.htail.D_box.volume(i) = (pyramid_big-pyramid_small) - D_box_volume0; 
end

local_mass_layout = y + 3*y/8 ;
N = size(y);
mass_location = -sum(local_mass_layout)/N(2) ; %location in reference to spar
construction.htail.D_box.volume_total = sum(construction.htail.D_box.volume) * 2 ;
construction.htail.D_box.mass(1,1) = mass_location;
construction.htail.D_box.mass(1,2) = construction.htail.D_box.volume_total * construction.material.balsa.density;

%-------------------- Ribs: 
% surface of the known rib (with known chord size) * scaling factor = our rib
scaling_rib = rib_chord/construction.htail.rib.NACA2412.known.chord;
construction.htail.rib.surface = scaling_rib .*scaling_rib *construction.htail.rib.NACA2412.known.surface;

local_mass_layout = rib_chord*construction.wing.rib.NACA2412.cg - y;
N = size(y);
mass_location = sum(local_mass_layout)/N(2); %location in reference to spar
ribs_volume = sum(construction.htail.rib.surface)*construction.htail.rib.material.thickness;
construction.htail.rib.mass(1,1) = mass_location;
construction.htail.rib.mass(1,2) = ribs_volume*construction.material.balsa.density *2;


geo.htail.mass = construction.htail.spar.mass(1,2) + construction.htail.D_box.mass(1,2) + construction.htail.rib.mass(1,2);


%% Fuselage
% fuselage is considered as stringer + four frames + carbon tube

%-------------------Frames:
N = size(construction.fuselage.frame.xlocation_diameter);
construction.fuselage.frame.mass = zeros(N);
construction.fuselage.frame.mass(:,1) = construction.fuselage.frame.xlocation_diameter(:,1);

for i=1:N(1)
    Pbig = pi* (construction.fuselage.frame.xlocation_diameter(i,2)/2)^2;
    Psmall = pi*((construction.fuselage.frame.xlocation_diameter(i,2)-2*construction.fuselage.frame.zthickness)/2)^2;
    P_circle = Pbig-Psmall;
    Volume =  P_circle * construction.fuselage.frame.xthickness;
    construction.fuselage.frame.mass(i,2) = Volume * construction.material.balsa.density;
end


%-------------------Carbon tube:
split = 30;
construction.fuselage.tail.mass = zeros(split,2);
construction.fuselage.tail.mass(:,1) = linspace(construction.fuselage.frame.xlocation_diameter(N(1)-1,1),geo.fuselage.length,split);
middle_of_section = (construction.fuselage.tail.mass(2,1)-construction.fuselage.tail.mass(1,1))/2;
N=split-1;
Pbig = pi* (construction.fuselage.tail.Diameter/2)^2;
Psmall = pi*((construction.fuselage.tail.Diameter-2*construction.fuselage.tail.thickness)/2)^2;
P_circle = Pbig-Psmall;
Volume =  P_circle * middle_of_section*2;
construction.fuselage.tail.mass(:,1) = construction.fuselage.tail.mass(:,1) + middle_of_section;

for i=1:N
    construction.fuselage.tail.mass(i,2) = Volume * construction.material.CFRP.density;
end

%-------------------Stringer:
 %P = construction.fuselage.stringer.ythickness*construction.fuselage.stringer.zthickness;

 geo.fuselage.mass = sum(construction.fuselage.frame.mass(: ,2)) + sum(construction.fuselage.tail.mass(:,2));

end
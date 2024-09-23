% Authors:
% Albert Canovas Cots
% Natalia Zalewska
function [allData, allData_sorted, CG] = xBalance(geo,construction,tailcoef)
% Finds CG
%   Based on wanted CG, finds wing spar position

%% Masses in reference to Wing spar

vtail_spar_loc = geo.wing.xle + tailcoef.lh + (construction.vtail.spar.location) * geo.vtail.cr
htail_spar_loc = geo.wing.xle + tailcoef.lh + (construction.vtail.spar.location) * geo.htail.cr

%geo.fuselage.length - (1-construction.vtail.spar.location) * geo.vtail.cr;
%htail_spar_loc = geo.fuselage.length - (1-construction.htail.spar.location) * geo.htail.cr;


% Define the location and mass data for each component
wing_spar = {'Wing spar', [construction.wing.spar.location_glob], [construction.wing.spar.mass(1,2)]};
wing_Dbox = {'Wing D-box', [construction.wing.D_box.mass(1,1)+construction.wing.spar.location_glob], [construction.wing.D_box.mass(1,2)]};
wing_rib = {'Wing rib', [construction.wing.rib.mass(1,1)+construction.wing.spar.location_glob], [construction.wing.rib.mass(1,2)]};
wing_joint = {'Wing joint', [construction.wing.joint.mass(1,1)+construction.wing.spar.location_glob], [construction.wing.joint.mass(1,2)]};

vtail_spar = {'Vtail spar', vtail_spar_loc, [construction.vtail.spar.mass(1,2)]};
vtail_Dbox = {'Vtail D-box', [vtail_spar_loc + construction.vtail.D_box.mass(1,1)], [construction.vtail.D_box.mass(1,2)]};
vtail_rib = {'Vtail rib', [vtail_spar_loc + construction.vtail.rib.mass(1,1)], [construction.vtail.rib.mass(1,2)]};

htail_spar = {'Htail spar', [htail_spar_loc], [construction.htail.spar.mass(1,2)]};
htail_Dbox = {'Htail D-box', [htail_spar_loc + construction.htail.D_box.mass(1,1)], [construction.htail.D_box.mass(1,2)]};
htail_rib = {'Htail rib', [htail_spar_loc + construction.htail.rib.mass(1,1)], [construction.htail.rib.mass(1,2)]};

fuselage_frame1 = {'Fuselage frame 1', [construction.fuselage.frame.mass(1,1)], [construction.fuselage.frame.mass(1,2)]};
fuselage_frame2 = {'Fuselage frame 2', [construction.fuselage.frame.mass(2,1)], [construction.fuselage.frame.mass(2,2)]};
fuselage_frame3 = {'Fuselage frame 3', [construction.fuselage.frame.mass(3,1)], [construction.fuselage.frame.mass(3,2)]};
fuselage_frame4 = {'Fuselage frame 4', [construction.fuselage.frame.mass(4,1)], [construction.fuselage.frame.mass(4,2)]};
flight_controller = {'Flight controller', [0.15], [55e-3]};
battery = {'Battery', [0.2], [1140e-3]};
motors = {'Motors', [0.05], [100e-3]};
ESC = {'Electronic Speed Controllers', [0.2], [20e-3]};
GPS = {'GPS module', [0.2], [30e-3]};
receiver = {'Receiver', [0.2], [12e-3]};
processor = {'Processor', [0.2], [40e-3]};
camera = {'Camera', [0.2], [160e-3]};
transmitter = {'Transmitter', [0.15], [30e-3]};
antennas = {'Antennas', [0.15], [12e-3]};
acuators_bottle = {'Water bottle acuators', [construction.fuselage.bottle.mass(1,1)+construction.wing.spar.location_glob], [120e-3]};
acuators_control = {'Flight control acuators', [0.7], [80e-3]};
propeller = {'Propeller', [0.05], [25e-3]};
%lg = {'Landing Gears', [], [700e-3]};
%container = {'Payload structure', [], [150e-3]};
%payload = {'Payload', [construction.fuselage.payload.mass(1,1)+construction.wing.spar.location_glob], [construction.fuselage.payload.mass(1,2)]};
%bottles = {'Bottles', [construction.fuselage.bottle.mass(1,1)+construction.wing.spar.location_glob], [construction.fuselage.bottle.mass(1,2)]};

% Combine the component data into a single cell array
allData = [wing_spar; wing_Dbox; wing_rib; wing_joint; vtail_spar; vtail_Dbox; vtail_rib; htail_spar; ...
    htail_Dbox; htail_rib; fuselage_frame1; fuselage_frame2; fuselage_frame3; fuselage_frame4; ...
    flight_controller; battery; motors; ESC; GPS; receiver; processor; camera; transmitter; antennas; ...
    acuators_bottle; acuators_control; propeller];

% Add the 30x2 matrix for fuselage
%allData = [allData, fuselage_tail_label, fuselage_frames_label];

% Sort the cell array based on location data
allData_sorted = sortrows(allData,2);

%Calc cg:
Data_sorted = cell2mat(allData_sorted(:,2:3));
CG = sum(Data_sorted(:,1) .*Data_sorted(:,2))/sum(Data_sorted(:,2));

% Extract the sorted data
%sorted_locations = allData_sorted(:, 1);
%sorted_names = allData_sorted(:, 2);
%sorted_masses = allData_sorted(:, 3:end);

% Convert the names back to a cell array
%sorted_names = mat2cell(sorted_names, ones(size(sorted_names,1),1), size(sorted_names,2));

end
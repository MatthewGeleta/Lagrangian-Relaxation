% MATLAB script:
%   -Obtains coordinates of Australian towns from .csv file.
%   -Filters outliers and Western half of country.
%   -Defines an approximate fire-hazard region.
%   -Computes distance matrix using great-circle distances and the
%       Haversine formula.
%   -Displays data on geographical map

% --- Parameters
% File name
fname = 'Austrlian_Towns_and_Coordinates.csv';
% Longitude filters
MinLong = 130;
% Latitude filters
MinLat = -39; MaxLat = -15;

% Read in data from all Australian Towns
Aus_Coords = csvread(fname,0,2);

% Get latitude and longitude info
Lats = Aus_Coords(:,1);
Longs = Aus_Coords(:,2);

% Filter incorrect geographical coordinates
[Lats,I] = sort(Lats);
Longs = Longs(I);
I = Longs > MinLong;
Longs = Longs(I);
Lats = Lats(I);
% Filter out regions off mainland Australia
I = MinLat < Lats;
Longs = Longs(I);
Lats = Lats(I);
I = Lats < MaxLat;
Longs = Longs(I);
Lats = Lats(I);
N_towns = length(Lats);

% Function defining approximate location of dessert region %144,145.5
f_desert = @(lat,lng) -37<lat & lat <-30 & 142<lng & lng <147;
% Get indices corresponding to approximate desert region
I_desert = f_desert(Lats,Longs);

% Get Southern-most and Northern-most points
Lat_South = Lats(1);
Long_South = Longs(1);
Lat_North = Lats(end);
Long_North = Longs(end);
% Northern- and Southern-most points
N_North = length(Lats);
N_South = 1;
N_Start = 1113; % Start node (Fowler's Bay, Southern Australia)
N_End = N_North - 3; % End node (Palm Cove, Queensland (North-East))

% Compute distance matrix
Aus_Dist_Mat_East = Haversine(Lats, Longs);
% Matrix Des(i,j) = 1 if (i,j) or (j,i) enters desert, 0 otherwise
Des = zeros(N_towns);
for i = 1:N_towns
    for j = 1:N_towns
        Des(i,j) = I_desert(i)*I_desert(j);
    end
end

% Save all data to .mat file
save('Aus_Coords_East')


myfig = figure();
landareas = shaperead('landareas.shp','UseGeoCoords',true);
geoshow(landareas,'FaceColor',[0.5 1.0 0.5],'EdgeColor',[.6 .6 .6]);
geoshow(Lats,Longs, 'DisplayType', 'point')
% Desert region
geoshow(Lats(I_desert), Longs(I_desert),'Color','k');
title('Australian East-Coast towns')
xlabel('Latitude')
ylabel('Longitude')
%tissot;



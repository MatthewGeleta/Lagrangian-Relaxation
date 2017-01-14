function D = Haversine(lats, longs)
% Returns distance matrix calculated using great-circle distances,
%   computer with the Haversine formula.

% Radius of earth in km
radius=6371;
% Number of points
N = length(lats);
D = zeros(N);
% Scale to radians
lats = lats*pi/180;
longs = longs*pi/180;

for ii = 1:N
    for jj = 1:N
        dLat=lats(jj)-lats(ii);
        dLon=longs(jj)-longs(ii);
        a = sin(dLat/2)^2 + cos(lats(jj))*cos(lats(ii)) * sin(dLon/2)^2;
        c=2*atan2(sqrt(a),sqrt(1-a));
        D(ii,jj) = c;
    end
end
D = D*radius;
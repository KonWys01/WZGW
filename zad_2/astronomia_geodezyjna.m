clear;

% rektascensja i deklinacja 
% gwiazdy Shaula z gwiazdozbioru Skorpiona
rekt = 17 + 35/60 + 3.4/3600; 
dek = -(37 + 7/60 + 4.8/3600);

% Sztokholm
phi = 59.3325800;
lambda = 18.0649000;

% % Mogadiszu - Stolica Somalii
% phi = 2.0371100;
% lambda = 45.3437500; 

% % Hobart - Stolica Tasmanii
% phi = -42.8793600;
% lambda = 147.3294100;
        
h = [0 1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]';

% Kąt dla kolejnych godzin
hour_angle = katgodz(2022, 01, 4, h, lambda, rekt);

% Poprawa kąta
hour_size = size(hour_angle, 1);
for i = 1:hour_size
    if hour_angle(i) > 360
        hour_angle(i) = hour_angle(i) - 360;
    end
end

% odleglosc zenitalna
cos_z = sind(phi).*sind(dek) + cosd(phi).*cosd(dek).*cosd(hour_angle);
z = acosd(cos_z);

% Azymut gwiazdy
nominator = -cosd(dek).*sind(hour_angle);
denominator = cosd(phi).*sind(dek) - sind(phi).*cosd(dek).*cosd(hour_angle);
A = atan2d(nominator, denominator);
for i=1:length(A)
   if A(i) < 0 
       A(i) = A(i) + 360; 
   elseif A(i) > 360
       A(i) = A(i) - 360;
   end
end 

% transformacja współrzednych
x = 1.4.*sind(z).*cosd(A);
y = 1.4.*sind(z).*sind(A);
z = 1.4.*cosd(z);

% Rysowanie półkuli 
[X,Y,Z] = sphere(16);
X = X(9:end,:);
Y = Y(9:end,:);
Z = Z(9:end,:);
surf(X,Y,Z,'FaceColor','black','FaceAlpha',0.3)
axis equal, hold on;
% Rysowanie gwiazdy
scatter3(x,y,z, 160, 'yellow', '*')

% wysokosc od czasu
plot(h, 90-z)

% zenit od czasu
plot(h,z)


function [t] = katgodz(y, m, d, h, lambda, alfa)
    jd = juliandate(datetime(y, m, d)); % dni
    g = GMST(jd); % stopnie
    UT1 = h * 1.002737909350795; % godziny

    % obliczenie czasu gwiazdowego (w stopniach)
    S = UT1*15 + lambda + g; 
    % obliczenie kąta godzinowego (w stopniach)
    t = S - alfa*15;
end

function g = GMST(jd)
    T = (jd - 2451545) / 36525;
    g = 280.46061837 + 360.98564736629  * (jd - 2451545.0) + 0.000387933*T.^2 - T.^3/38710000;
    g = mod(g, 360);
end
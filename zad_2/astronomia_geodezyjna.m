clear;

% rektascensja i deklinacja 
% gwiazdy Shaula z gwiazdozbioru Skorpiona
rekt = 17 + 35/60 + 3.4/3600; 
dek = -(37 + 7/60 + 4.8/3600);
disp(rekt)
disp(dek)

% Sztokholm
phi = 59.3325800;
lambda = 18.0649000;

% % Mogadiszu - Stolica Somalii
% phi = 2.0371100;
% lambda = 45.3437500; 

% % Hobart - Stolica Tasmanii
% phi = -42.8793600;
% lambda = 147.3294100;

h = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24]';
disp(h)

% Kąt dla kolejnych godzin
hour_angle = katgodz(2022, 01, 04, h, lambda, rekt);
disp(hour_angle)

% Poprawa kąta
hour_size = size(hour_angle, 1);
for i = 1:hour_size
    if hour_angle(i) > 360
        hour_angle(i) = hour_angle(i) - 360;
    end
end

% odleglosc zenitalna
cos_z = sind(phi).*sind(dek) + cosd(phi).*cosd(hour_angle);
z = acosd(cos_z);

% Azymut gwiazdy
tg_A = (-cosd(dek).*sind(hour_angle))./(cosd(phi).*sind(dek) - sind(phi).*cosd(dek).*cosd(hour_angle));
A = atand(tg_A);

% transformacja współrzednych
x = 1.*sind(z).*cosd(A);
y = 1.*sind(z).*sind(A);
z = 1.*cosd(z);












% d = datetime('now');
% d = datetime(2016, 07, 30);
% jd = juliandate(d);
% disp(jd)
% disp(datetime(jd,'ConvertFrom','juliandate'))
% katgodz(2021, 10, 10, 10, 10, 10)
function [t] = katgodz(y, d, m, h, lambda, alfa)
    jd = juliandate(datetime(y, m, d)); % dni
%     disp('jd=')
    g = GMST(jd); % stopnie
    UT1 = h * 1.002737909350795; % godziny

    % obliczenie czasu gwiazdowego (w stopniach)
    S = UT1*15 + lambda + g; 
    % obliczenie kąta godzinowego (w stopniach)
    t = S - alfa*15;
end

function g = GMST(jd)
    T = (jd - 2451545) / 36525;
    g = 280.4606837 + 360.98564736629 * (jd - 2451545.0) + 0.00387933*T.^2 - T.^3/3870000;
    g = mod(g, 360);
end




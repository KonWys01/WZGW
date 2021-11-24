clear;

% GRS 80
a=6378137;
e2=0.00669437999013;

% Dane lotniska początkowego
phiB = 52.45066;
lambdaB = 20.65099;
hB = 341;

% Dane dot. lotu
dane = load('dane.txt');
phi = dane(:,1);
lambda = dane(:,2);
h = dane(:,3);
rows = size(dane, 1);

% XYZ samolotu
N = (a./sqrt(1-e2 .* sind(phi) .* sind(phi)));
x = ((N+h) .* cosd(phi) .* cosd(lambda));
y = ((N + h) .* cosd(phi) .* sind(lambda));
z = ((N .* (1-e2) + h) .* sind(phi));

% XYZ lotniska
NB = (a./sqrt(1-e2 .* sind(phiB) .* sind(phiB)));
xB = ((NB+hB) .* cosd(phiB) .* cosd(lambdaB));
yB = ((NB + hB) .* cosd(phiB) .* sind(lambdaB));
zB = ((NB .* (1-e2) + hB) .* sind(phiB));

% NEU
obrotu = [-sind(phiB).*cosd(lambdaB) -sind(lambdaB) cosd(phiB).*cosd(lambdaB)
          -sind(phiB).*sind(lambdaB) cosd(lambdaB) cosd(phiB).*sind(lambdaB)
          cosd(phiB) 0 sind(phiB)];

transponowana = transpose(obrotu);
delty = [(x-xB)'
         (y-yB)'
         (z-zB)'];
neu = transponowana * delty;

n = neu(1,:)';
e = neu(2,:)';
u = neu(3,:)';

% Odległość skośna
s = sqrt(n.^2 + e.^2 + u.^2);
km = s ./ 1000;

% Wyznaczanie kąta zenitalnego
cos_z = u ./ s;
zenit = acosd(cos_z);

% Azymut A
tan_A = e ./ n;
azymut = atand(tan_A);

% Dostosowanie azymutu względem ćwiartki
for i = 1:rows
    if n(i) > 0 && e(i) < 0
        azymut(i) = azymut(i) + 360;
    elseif ((n(i) < 0 && e(i) > 0) || (n(i) < 0 && e(i) < 0))
        azymut(i) = azymut(i) + 180;
    end
end

% lot z phi i lambda
% geoscatter(phi,lambda,5,'ro');

% lot XYZ
% plot3(x,y,z);

% lot NEU
% plot3(n,e,u);
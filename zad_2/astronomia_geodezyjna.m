clear;

% zamiana na kalendarz juliański
% function jul_date = julday(y,m,d,h)
%     if m <= 2
%         y = y - 1;
%         m = m + 12;
%     end
%     julian_date = floor(365.25*(y+4716))+floor(30.6001*(m+1))+d+h/24-1537.5; 
% end


d = datetime('now');
% d = datetime(2016, 07, 30);
% jd = juliandate(d);
% disp(jd)
disp(datetime(jd,'ConvertFrom','juliandate'))
function [t] = katgodz(y, d, m, h, lambda, alfa)
    jd = juliandate(datetime(y, m, d, 0)); % dni
    g = GMST(jd); % stopnie
    UT1 = h * 1.002737909350795; % godziny

    % obliczenie czasu gwiazdowego (w stopniach)
    S = UT1*15 + lambda + g; 
    % obliczenie kąta godzinowego (w stopniach)
    t = S - alfa*15;
end

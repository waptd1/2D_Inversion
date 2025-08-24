function  e = Misfit(pm,c_obs)

% Number of data points in theoretical and experimental dispersion curves.
Q = length(pm); 

b=~isnan(pm);

% misfit
e = (1/Q)*sum(abs(c_obs(b)-pm(b))./c_obs(b))*100;

end

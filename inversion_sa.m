%% Initial parameters

lb = input('Lower bounds on profile = ');
ub = input('upper bounds on profile = ');

x0 = input('Initial profile = ');

max_iter = input('Maximum iterations (provide between 1000-5000) = ');
ini_temp = input('Initial temperature (provide between 50-250) = ');
reanl_int = input('Reanneal interval (provide between 50-150) = ');

%% Options and inversion

options=optimoptions(@simulannealbnd,'Display','iter','PlotFcn','saplotbestf','MaxIterations',max_iter);
options=optimoptions(options,'InitialTemperature',ini_temp,'TemperatureFcn','temperatureboltz','ReannealInterval',reanl_int); 

[x,misfit,stop_reason,output] = simulannealbnd(@Inversion_CFEM,x0,lb,ub,options);
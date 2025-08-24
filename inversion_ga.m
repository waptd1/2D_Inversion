%% Initial parameters

lb = input('Lower bounds on profile = ');
ub = input('upper bounds on profile = ');
nvar = length(lb);

A=[]; b=A; Aeq=A; beq=A; nonlcon=A;
gens = input('Number of generations (provide between 50-100) = ');
xcf = input('Crossover fraction (provide between 0.5-0.9) = ');

%% Set GA options

options = optimoptions('ga','UseParallel', true,'Display','iter','PlotFcn','gaplotbestf','OutputFcn',@gaoutfun);
options = optimoptions(options,'Generations',gens,'CrossoverFraction',xcf);

rng default;

%% Inversion and output 

[x,misfit,stop_reason,output] = ga(@Inversion_CFEM,nvar,A,b,Aeq,beq,lb,ub,nonlcon,options);
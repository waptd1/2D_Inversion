%% Input 

sg=input('experimental or synthetic signal matrix = ');

fs=input('acquistion frequency = ');

sz=size(sg); M=sz(2);
so=input('source offset = '); 
dx=input('receiver interval ='); 
x=so+((1:M)-1)*dx;

dt=1/fs;
td=sz(1)*dt;
t=(0:dt:td-dt)';

%% time intercept of actual signal 

% By plotting the signal at first receiver, find manually the time at which the actual impact signal is received

figure; plot(sg(:,1),t);
hold on;
plot([0,0],[t(1),t(end)],'k')
ax=gca;
ax.YDir='reverse';

% for example
tc=find(round(t,4)==.2245); % where 0.2245 is the time of signal received at 1st geophone

%% Vertical displacement

sg1(:,1:sz(2))=sg(:,1:end)-mean(sg(1:tc,1:end));

sg2=zeros(sz(1),sz(2));

for j=1:sz(2)
    for i=2:sz(1)
        sg2(i,j)=sg2(i-1,j)+(sg1(i-1,j)+sg1(i,j))*0.5*dt;
    end
end


mxsg2=min((sg2)); % maximum vertical downward displacement

%% for synthetic data
mxsg2_fea=mxsg2;

%% for experimental data

% sens is only for experimental data (For RTC 4.5HZ 395OHM VERTICAL GEOPHONE, sens = 23.4)
sens=23.4; % V/m/s = Vs/m. 

mxsg2_exp=mxsg2/sens;

%% Plot

figure;
hold on;
% plt1=plot(f,phv.','-k'); % to plot all rows at once

idx=1;

plt1=plot(x(1:1:end),mxsg2_exp(1:1:end));
lw=0.7;
plt1(idx).LineWidth =lw;plt1(idx).LineStyle='-'; plt1(idx).Color='k';
plt1(idx).Marker='o'; plt1(idx).MarkerSize=2.5; 
plt1(idx).MarkerEdgeColor='k';plt1(idx).MarkerFaceColor='k';
plt1(1).MarkerIndices = 1:2:length(x);

plt3=plot(x(1:1:end),mxsg2_fea(1:1:end));
% idx=1;
lw=1;
plt3(idx).LineWidth =lw; plt3(idx).LineStyle='--'; plt3(idx).Color='k';
plt3(idx).Marker='o'; plt3(idx).MarkerSize=3.5; 
plt3(idx).MarkerEdgeColor='k';plt3(idx).MarkerFaceColor='w';
plt3(idx).MarkerIndices = 1:2:length(x);



xlabel('Distance from source (m)','FontName','Times','FontUnit','points','FontSize',12,...
    'FontWeight','bold','FontAngle', 'normal','Interpreter','tex');
%
ylabel('Downward maximum displacement (m)','FontName','Times','FontUnit','points','FontSize',13,...
    'FontWeight','bold','FontAngle', 'normal','Interpreter','tex');


%%
lgd=legend([plt1(1),plt3(1)],{'\rmExperimental','\rmDynamic FE analysis'});

lgd.Units='normalized';
lgd.Position=[.462 .615 .07 .07]; % for portrait for 3 lines in annotation
lgd.FontSize=11;
lgd.FontName='Times';
lgd.FontUnits='normal';
lgd.FontWeight='bold';

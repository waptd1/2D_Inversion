%% velocity with depth matrix calculation

x = input('Obtained inverted Vs profile/vector = ');
lr=(length(x)+1)/2;
vs_temp(1:lr,1)=[x(lr+1:end)'; 5]; % h on column 1
vs_temp(1:lr,2)=x(1:lr)'; % vs on column 2

vs_temp(:,3)=vs_temp(:,1); % for consequent summation
for i=2:lr
    vs_temp(i,1)=vs_temp(i-1,1)+vs_temp(i,3);
end

b=size(vs_temp(:,1:2)); nb=2*b(1)-1; % nb is the new length of f and Cm matrix
vs(1:2:nb,:)=vs_temp(1:b(1),1:2); vs=[0 0;vs]; vs(1,2)=vs(2,2);

for i=1:2:(lr*2-3)
    vs(i+2,1)=vs(i+1,1);
    vs(i+2,2)=vs(i+3,2);
end


%% vs15 calculation

index = find(vs_temp(:,1)>=15,1);

den=sum(vs_temp(1:index-1,3)./vs_temp(1:index-1,2)); % denominator
den=den+(15-vs_temp(index-1,1))/vs_temp(index,2); % denominator last term

vs15=15/den;

%% Plotting Vs profile 

figure;fg=plot(vs(:,2),vs(:,1)); fg.Color=[0 0 0]; hold on; 

% 2. fig dimensions
set(gcf,'Units', 'centimeters'); %setting units
afFigurePosition=[20 15 6.5 5.1]; %setting fig pos [pos_x pos_y width_x width_y] % for landscape
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto'); %link the dim of the fig ON THE PAPER in such a way that it is equal to the dim on the screen


% 3. Axes
ax=gca; 
ax.Units='normalized';
ax.Position=[.162 .06 .78 .73]; % for landscape ... manipulated
axpos = ax.Position;
ax.TickDir='out';
ax.TickLength=[.005 .005];
ax.XColor='k';
ax.YColor='k';
ax.XAxisLocation='top';
ax.YAxisLocation='left';
ax.XDir='normal'; % ... axis increasement direction [{normal} | reverse]
ax.YDir='reverse'; %... axis increasement direction [{normal} | reverse]
ax.FontName='Times'; %... kind of fonts of labels
ax.FontSize=10;  %... size of fonts of labels
ax.FontUnits='points'; %... units of the size of fonts
ax.FontWeight='bold'; %... weight of fonts of labels
ax.FontAngle='normal';  %... inclination of fonts of labels
ax.LineWidth=0.5;

xlabel('Shear wave velocity (m/s)','FontName','Times','FontUnit','points','FontSize',11,...
    'FontWeight','bold','FontAngle', 'normal','Interpreter','tex');
ylabel('Depth (m)','FontName','Times','FontUnit','points','FontSize',11,...
    'FontWeight','bold','FontAngle', 'normal','Interpreter','tex');
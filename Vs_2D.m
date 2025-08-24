%% 1. input

x = input('Obtained inverted Vs profile/vector = ');
lr=(length(x)+1)/2;


%% 2. Vs profile

vs1(1:lr,1)=[x(lr+1:end)'; 5]; % h on column 1
vs1(1:lr,2)=x(1:lr)'; % vs on column 2
vs1(:,3)=vs1(:,1); % for consequent summation
for i=2:lr
    vs1(i,1)=vs1(i-1,1)+vs1(i,3);
end
vs1(:,3)=[];
vs1(:,1)=round(vs1(:,1),2);

%% 3. x-y coordinates

max_depth = input('maximum depth of 2D profile, check from vs1 matrix (in m) = ');

x1=input('vector contaning distance of centre of profiles from source = ');
hint=linspace(0,max_depth,1000)';
vint=zeros(length(hint),length(x1));

%% 4. Vs interpolation in y-coordinate

vs1(lr,1)=max_depth;

hv=[0; vs1(:,1)];
vv=[vs1(:,2);vs1(lr,2)];

file=input('file_no.?='); % file_no is the Vs profile number from the start to end
vint(:,file)=interp1(hv, vv, hint, 'pchip', 'extrap');

%% Important note

% Repeat step 1-4 for all the inverted profiles for a given site until-
% vint contains all profiles on its columns.

%% 5. Vs interpolation in x-coordinate

xint= linspace(x1(1),x1(end),4*length(x1));
vint_interp=griddata(x1, hint, vint, xint, hint, 'cubic');

%% 6. 2D Vs profile plot

figure;fg=contourf(xint,hint,vint_interp, 'LineColor', 'none'); colormap jet; shading interp; cb=colorbar; axis xy; %view(2); 

% fig dimensions
set(gcf,'Units', 'centimeters'); %setting units
afFigurePosition=[20 5 13 8]; %setting fig pos [pos_x pos_y width_x width_y] % for portrait mode
set(gcf, 'Position', afFigurePosition); % [left bottom width height]
set(gcf, 'PaperPositionMode', 'auto'); %link the dim of the fig ON THE PAPER in such a way that it is equal to the dim on the screen


% 3. Axes
ax=gca; 
ax.Units='normalized';
ax.Position=[.08 .048 .79 .83]; % for portrait
axpos = ax.Position;
cb.Position(3) = 0.5*cb.Position(3);
cb.Position(1) = .99*cb.Position(1);
cmax=round(max(max(vint)),2);
cmin=round(min(min(vint)),2);
cb.Limits=[cmin cmax];
ax.TickDir='out';
ax.TickLength=[.004 .004];
ax.XColor='k';
ax.YColor='k';
ax.XAxisLocation='top';
ax.YAxisLocation='left';
ax.XDir='normal'; % ... axis increasement direction [{normal} | reverse]
ax.YDir='reverse'; %... axis increasement direction [{normal} | reverse]
ax.FontName='Times'; %... kind of fonts of labels
ax.FontSize=11;  %... size of fonts of labels
ax.FontUnits='points'; %... units of the size of fonts
ax.FontWeight='normal'; %... weight of fonts of labels
ax.FontAngle='normal';  %... inclination of fonts of labels
ax.LineWidth=0.5;
% 
xlabel('Offset to centre of inverted profiles (m)','FontName','Times','FontUnit','points','FontSize',12,...
    'FontWeight','bold','FontAngle', 'normal','Interpreter','tex');

ylabel('Depth (m)','FontName','Times','FontUnit','points','FontSize',12,...
    'FontWeight','bold','FontAngle', 'normal','Interpreter','tex');

ylabel(cb, 'Shear wave velocity (m/s)','FontName','Times','FontUnit','points','FontSize',12,...
    'FontWeight','bold','FontAngle', 'normal','Interpreter','tex');


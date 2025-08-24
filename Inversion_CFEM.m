%% This function is taken from the following references:
%
% http://www.vaziri.info/wavedisp/
% 
% 1. A Vaziri Astaneh, MN Guddati, Efficient computation of
% dispersion curves for multilayered waveguides and half-spaces,
% Comput. Meth. Appl. Mech. Eng. 300 (2016) 27-46.

% 2. A Vaziriastaneh, On the Forward and Inverse Computational
% Wave Propagation Problems, Ph.D. Thesis, North Carolina
% State University, USA, 2016.
% (chapter 2)

% 3. A Vaziri Astaneh, MN Guddati, Improved Inversion
% Algorithms for Near Surface Characterization,
% Geophys. J. Int.  206.2 (2016) 1410-1423.


function [misfit] = Inversion_CFEM(x)


%% Observed predominant mode

% Provide frequency vector
f_obs=[10.0097656250000	10.9863281250000	11.9628906250000	12.9394531250000	13.9160156250000	14.8925781250000	15.8691406250000	16.8457031250000	17.8222656250000	18.7988281250000	19.7753906250000	20.7519531250000	21.7285156250000	22.7050781250000	23.6816406250000	24.6582031250000	25.6347656250000	26.6113281250000	27.5878906250000	28.5644531250000	29.5410156250000	30.5175781250000	31.4941406250000	32.4707031250000	33.4472656250000	34.4238281250000	35.4003906250000	36.3769531250000	37.3535156250000	38.3300781250000	39.3066406250000	40.2832031250000	41.2597656250000	42.2363281250000	43.2128906250000	44.1894531250000	45.1660156250000	46.1425781250000	47.1191406250000	48.0957031250000	49.0722656250000	50.0488281250000];

% Provide predominant mode phase velocity vector
c_obs=[386	331	514	333	431	346	467	382	389	366	404	400	390	402	401	408	432	420	429	443	446	450	448	450	454	460	456	459	465	466	465	468	471	472	476	467	463	470	1053	1034	1025	969]; % M=48, so=5


%% Input parameters

n=floor(length(x)/2);

% Provide layers' Poisson's ratio vector
nu=ones(1,n+1)*0.3; 


% Provide layer density vector (kg/m^3)
rho=[1700 1700 1750 1750 1800 1800 1850 1850 1900 1950]; 

Pr.r0=  5;     % Set Source offset to nearest geophone
Pr.dr=  1;    % Set Receiver interval
Pr.rN=  48;    % Set Number of Receivers



Pr.cs= x(1:n);      % Layer Shear Wave Velocity (m/s)
Pr.cp(1:n)= Pr.cs(1:n)*sqrt(2*(1-nu(1:n))/(1-2*nu(1:n)));      % Layer Pressure Wave Velocity (m/s)

Pr.roS = rho(1:end-1);
Pr.h=     x(n+2:end);    % Layer Thickness (m)
Pr.nDivS= ones(1,n)*2;            % Number of Quartic (5-noded) Elements per Layer

%% BOTTOM HALF-SPACE [using PMDLs]

Pr.hsB=    'yes';    % Bottom Half-Space 'yes','no'
Pr.csB=    x(n+1) ;     % HS Shear Wave Velocity (m/s)
Pr.cpB=    Pr.csB*sqrt(2*(1-nu(n+1))/(1-2*nu(n+1)));     % HS Pressure Wave Velocity (m/s)
Pr.roB=    rho(end);     % HS Density (kg/m3)
Pr.nDivB=  10;       % Number of (linear) PMDL elements
Pr.pL1=    1;        % 1st Layer Lenght: Lj=  L1  * alpha ^ (j-1)
Pr.pAlp=   2;        % Increase Ratio:   Lj=  L1  * alpha ^ (j-1)

%% EFFECTIVE DISPERSION CURVE

Pr.eff= 'yes';  % Effective Dispersion Curve Calculation: 'yes','no'
Pr.pad= 2048;   % Number of Padding Layers in FFT Calculation (0: No Padding)
Pr.q=   1.4e6;    % Magnitude of Circular Distributed load (N/m2)   1.4e6 for 20lbs hammer with 15 cm dia plate;
Pr.R=   .075;%1e-2;   % Radius of Circular Distributed load (m)     1.2490e+07 for 3.3 lbs (1.5 kg) with 3 cm dia hammer

% DISPLAY SETTINGS

% Pr.w=      [11.0473632812500	12.0849609375000	13.1225585937500	14.1601562500000	15.1977539062500	16.2353515625000	17.2729492187500	18.3105468750000	19.3481445312500	20.3857421875000	21.4233398437500	22.4609375000000	23.4985351562500	24.5361328125000	25.5737304687500	26.6113281250000	27.6489257812500	28.6865234375000	29.7241210937500	30.7617187500000];    % Frequency (Hz)
Pr.w= f_obs;
Pr.cg=     'no';         % Group Velocity Calculation: 'yes','no'
Pr.kzTol=  1e-3;          % Filtering Tolerance for imag(eigenvalue)
Pr.cpMax=  max(x(1:n+1))*1.25;           % Maximum Phase Velocity (m/s)
Pr.cpMin=  min(x(1:n+1))*0.5;           % Minimum Phase Velocity (m/s)
Pr.trace=  'no';          % Trace Curves: 'yes','no'
Pr.nMode=  10;            % Number of Modes to Trace

%% SOLID

% MATERIAL

nenS=5;  nL=length(Pr.cs);   nelS=sum(Pr.nDivS); numnodS = nelS*(nenS-1)+1;
D=cell(3,nL);
Lz=zeros(3,2); Lz(1,1)=1; Lz(3,2)=1; Ly=zeros(3,2); Ly(2,2)=1; Ly(3,1)=1;
for i=1:nL
    Dmat=zeros(3,3);
    Dmat(1,1)=Pr.roS(i)*Pr.cp(i)^2; Dmat(2,2)=Dmat(1,1);  Dmat(3,3)=Pr.roS(i)*Pr.cs(i)^2;
    Dmat(1,2)=Pr.roS(i)*(Pr.cp(i)^2-2*Pr.cs(i)^2);        Dmat(2,1)=Dmat(1,2);
    D{1,i}=Lz'*Dmat*Lz;             D{2,i}=Lz'*Dmat*Ly;   D{3,i}=Ly'*Dmat*Ly;
end

% MESH

matelem=zeros(1,nelS);  x1=zeros(1,nelS*(nenS-1)+1); connS=zeros(nenS,nelS);
j=0;
for i=1:nL
    for k=1:Pr.nDivS(i)
        j=j+1;   matelem(j)=i;
        b=(j-1)*(nenS-1)+1; e=j*(nenS-1)+1;
        x1(b:e)=x1(b)+((1:nenS)-1)*((Pr.h(i)/Pr.nDivS(i))/(nenS-1));
    end
end
for i=1:nenS
    connS(i,1:nelS)=i:nenS-1:numnodS-(nenS-i);
end

% STIFFNESS

GP=[-sqrt(5+2*sqrt(10/7))/3 -sqrt(5-2*sqrt(10/7))/3 0 sqrt(5-2*sqrt(10/7))/3 sqrt(5+2*sqrt(10/7))/3;
    (322-13*sqrt(70))/900 (322+13*sqrt(70))/900 128/225 (322+13*sqrt(70))/900 (322-13*sqrt(70))/900];
Kzz=zeros(2*numnodS);Kzy=Kzz;Kyy=Kzz;M=Kzz;
j=0;
for i=1:nL
    for k=1:Pr.nDivS(i)
        j=j+1;
        Dzz=D{1,matelem(j)};  Dzy=D{2,matelem(j)}; Dyy=D{3,matelem(j)};
        nde=connS(:,j); ndeE=[2*nde-1 2*nde]'; ndeE=ndeE(:); % for elasticity 2 dof per node
        for gpos=GP
            s=gpos(1); wt=gpos(2);
            N=[-(s*(- 4*s^3 + 4*s^2 + s - 1))/6;(4*s*(- 2*s^3 + s^2 + 2*s - 1))/3;
                4*s^4 - 5*s^2 + 1; (4*s*(- 2*s^3 - s^2 + 2*s + 1))/3;
                -(s*(- 4*s^3 - 4*s^2 + s + 1))/6].';
            N=kron(N,eye(2));
            Bs=[ (8*s^3)/3 - 2*s^2 - s/3 + 1/6; - (32*s^3)/3 + 4*s^2 + (16*s)/3 - 4/3;
                16*s^3 - 10*s; - (32*s^3)/3 - 4*s^2 + (16*s)/3 + 4/3;
                (8*s^3)/3 + 2*s^2 - s/3 - 1/6].';
            L=Pr.h(matelem(j))/Pr.nDivS(matelem(j));
            detJ=L/2;      B=Bs/detJ; B=kron(B,eye(2));
            Kzz(ndeE,ndeE)=Kzz(ndeE,ndeE)    + (N.'*Dzz*N)*(detJ*wt);
            Kzy(ndeE,ndeE)=Kzy(ndeE,ndeE)    - (N.'*Dzy*B)*(detJ*wt);
            Kzy(ndeE,ndeE)=Kzy(ndeE,ndeE)    + (B.'*Dzy.'*N)*detJ*wt;
            Kyy(ndeE,ndeE)=Kyy(ndeE,ndeE)    + (B.'*Dyy*B)*(detJ*wt);
            M(ndeE,ndeE)=M(ndeE,ndeE)        + (N.'*eye(2)*N)*(detJ*wt*Pr.roS(matelem(j)));
        end
    end
end

% BOTTOM HALF-SPACE -------------------------------------------------------

if strcmp(Pr.hsB,'yes')
    Dmat=zeros(3,3);
    Dmat(1,1)=Pr.roB*Pr.cpB^2; Dmat(2,2)=Dmat(1,1);  Dmat(3,3)=Pr.roB*Pr.csB^2;
    Dmat(1,2)=Pr.roB*(Pr.cpB^2-2*Pr.csB^2);          Dmat(2,1)=Dmat(1,2);
    Dzz=Lz'*Dmat*Lz;             Dzy=Lz'*Dmat*Ly;    Dyy=Ly'*Dmat*Ly;
    nenB=2; nelB=Pr.nDivB;  numnodB = nelB*(nenB-1)+1;  connB=zeros(nenB,nelB);
    Lpmdl=(1)*Pr.pL1*Pr.pAlp.^(0:nelB-1);
    xFo=cumsum([0 Lpmdl])+sum(Pr.h);
    for i=1:nenB
        connB(i,1:nelB)=i:nenB-1:numnodB-(nenB-i);
    end
    connB=connB+numnodS-1; GP=[0;2];
    id=(numnodS+numnodB-1)*2;Kzz(id,id)=0;Kzy(id,id)=0;Kyy(id,id)=0;M(id,id)=0;
    for k=1:Pr.nDivB
        nde=connB(:,k); ndeE=[2*nde-1 2*nde]'; ndeE=ndeE(:); % for elasticity 2 dof per node
        for gpos=GP
            s=gpos(1); wt=gpos(2); L= Lpmdl(k);
            N=1/2*[(1-s) (1+s)];  N=kron(N,eye(2)); Bs=1/2*[-1 1];
            detJ=L/2;             B=Bs/detJ;        B=kron(B,eye(2));
            Kzz(ndeE,ndeE)=Kzz(ndeE,ndeE)    + (N.'*Dzz*N)*(detJ*wt);
            Kzy(ndeE,ndeE)=Kzy(ndeE,ndeE)    - (N.'*Dzy*B)*(detJ*wt);
            Kzy(ndeE,ndeE)=Kzy(ndeE,ndeE)    + (B.'*Dzy.'*N)*detJ*wt;
            Kyy(ndeE,ndeE)=Kyy(ndeE,ndeE)    + (B.'*Dyy*B)*(detJ*wt);
            M(ndeE,ndeE)=M(ndeE,ndeE)        + (N.'*eye(2)*N)*(detJ*wt*Pr.roB);
        end
    end
    Kzz=Kzz(1:end-2,1:end-2); Kzy=Kzy(1:end-2,1:end-2); Kyy=Kyy(1:end-2,1:end-2); M=M(1:end-2,1:end-2);
end

Z=zeros(size(Kzz,1)/2); z=1:2:size(Kzz,1)-1; y=2:2:size(Kzz,1);
K2=([Kzz(z,z) Z; -Kzy(y,z) Kzz(y,y)]);
K0=([Kyy(z,z) Kzy(z,y);Z Kyy(y,y)]);
M=([M(z,z) Z;Z M(y,y)]);
clear Kzz Kzy Kyy
sizeS=size(K2,1);

% EIGENVALUE PROBLEM ------------------------------------------------------

kz=zeros(sizeS,length(Pr.w));
if (strcmp(Pr.cg,'yes') || strcmp(Pr.eff,'yes'))
    evcR=zeros(sizeS,sizeS,length(Pr.w)); evcL=evcR;
end
nrm=norm(K2(:)); K2=K2/nrm; K0=K0/nrm;  M=M/nrm;
for i=1:length(Pr.w)
    w=Pr.w(i)*2*pi;
%     clc; disp ([num2str(ceil(i/length(Pr.w)*100)) '% Completed']);
    if (strcmp(Pr.cg,'no') && strcmp(Pr.eff,'no'))
        kzi = sqrt(eig(K0-w^2*M, -K2));
    elseif (strcmp(Pr.cg,'yes') || strcmp(Pr.eff,'yes'))
        [evcRi,kzi,evcLi]  = eig(K0-w^2*M, -K2);  kzi=sqrt(diag(kzi));
        evcR(:,:,i)=(evcRi);evcL(:,:,i)=(evcLi);
    end
    kz(1:sizeS,i)=kzi;
end
if strcmp(Pr.eff,'yes'); kzE=kz; end

% FILTERING ---------------------------------------------------------------

for i=1:length(Pr.w)
    w=Pr.w(i)*2*pi;     kzi=kz(:,i);
    kzi(kzi==Inf)=0;   kzi(real(kzi)<0)=0;
    kzi(abs(imag(kzi))>Pr.kzTol)=0; kzi=real(kzi);
    kz(:,i)=kzi;
end

% GROUP VELOCITY-----------------------------------------------------------

if strcmp(Pr.cg,'yes')
%     disp('Group Velocity Calculation')
    cg=zeros(size(kz));
    for i=1:length(Pr.w)
        w=Pr.w(i)*2*pi;
        kzi=kz(:,i);
        if (strcmp(Pr.hsB,'yes'))
            kzi(real(kzi)<w/Pr.csB)=0;
        end
        [iloc]=find(kzi~=0);
        for j=iloc'
            cg(j,i)= real ( - ( evcL(:,j,i)' *(2*kzi(j)*K2)* evcR(:,j,i) )/ ...
                ( evcL(:,j,i)' *(-2*w*M)* evcR(:,j,i)) ) ;
        end
    end
end

% EFFECTIVE CURVE ---------------------------------------------------------

if strcmp(Pr.eff,'yes')
%     disp('Effective Curve Calculation')
    cpE=zeros(1,size(w,2));
    r=(Pr.r0 : Pr.dr : Pr.r0 + (Pr.rN-1) * Pr.dr);
    xmax=max(r);xmin=min(r);dx=Pr.dr;  x1=(xmin:dx:xmax-dx);
    dk = 2*pi/(xmax-xmin); kmin = -pi/dx; kmax = pi/dx;
    k = (kmin:dk:kmax-dk);
    for i=1:length(Pr.w)
        w=Pr.w(i)*2*pi; kzi=kzE(:,i);   evcRi=evcR(:,:,i);
        evcLi=evcRi*0;evcLi(1:end/2,:)=evcRi(1:end/2,:)*diag(kzi);
        evcLi(end/2+1:end,:)=evcRi(end/2+1:end,:)*diag(1./kzi);
        nr= diag(evcLi.'*(K2)*evcRi)./kzi(:);
        evcRi=evcRi*diag(1./sqrt(nr));evcLi=evcLi*diag(1./sqrt(nr));
        fiz=evcRi(1:end/2,:);fiy=evcLi(end/2+1:end,:);
        evcRi(1:end/2,:,i)=fiz;  evcRi(end/2+1:end,:,i)=fiy;
        kzi=kz(:,i);
        idk=find(kzi>0);
        uxw=zeros(1,length(r));
        bsl=zeros(1,length(idk));
        for j=1:length(idk)
            kzin=kzi(idk(j));
            bsl(j)=besselj(1,kzin*Pr.R);
        end
        j=0;
        for ir=r
            u=0;  j=j+1;
            for l=1:length(idk)
                kzin=kzi(idk(l));  bsl(l)=besselj(1,kzin*Pr.R);
                phi=real(evcRi(end/2+1,idk(l),i)^2);
                u=u+(phi*bsl(l)*besselh(0,2,kzin*ir)*(-1i*Pr.q*Pr.R*pi /(2*kzin)));
            end
            uxw(j)=u;
        end
        uxw(end)=[];
        if Pr.pad~=0
            Pr.rN=Pr.pad;
            xmax= Pr.r0 + (Pr.rN-1) * Pr.dr;
            dk = 2*pi/(xmax-xmin); kmin = -pi/dx; kmax = pi/dx;
            k = (kmin:dk:kmax-dk);             uxw(Pr.rN-1)=0;
        end
        ukw=fftshift(fft(uxw)); [~,kE]=max(abs(ukw)); cpE(1,i)=w/abs(k(kE));
    end
end

%% PLOT  -------------------------------------------------------------------

%{
if strcmp(Pr.trace,'no')
    figure(2);
    for i=1:size(kz,1)
        [iloc]=find(kz(i,:)~=0);
        pT=plot(Pr.w(iloc),2*pi*Pr.w(iloc)./real(kz(i,iloc)),'b.-'); hold on;
    end
    if strcmp(Pr.eff,'yes')
        [iloc]=find(cpE~=0);
        plot(Pr.w(iloc),cpE(iloc),'ro');
        title('B: Theoretical    R: Effective');
    end
    xlabel('Frequency (Hz)');ylabel('Phase Velocity (m/s)');
    ylim([Pr.cpMin Pr.cpMax]); box on;
    figure(3);
    for i=1:size(kz,1)
        [iloc]=find(kz(i,:)~=0);
        plot(Pr.w(iloc),real(kz(i,iloc)),'k.'); hold on;
    end
    xlabel('Frequency (Hz)');ylabel('Wavenumber (rad/m)'); box on;
    if strcmp(Pr.cg,'yes')
        figure(4);
        for i=1:size(kz,1)
            [ilocK]=find(kz(i,:)~=0); [ilocG]=find(cg(i,:)~=0);
            iloc=intersect(ilocK,ilocG);
            plot(Pr.w(iloc),cg(i,iloc),'r.'); hold on;
        end
        xlabel('Frequency (Hz)');ylabel('Group Velocity (m/s)'); box on;
    end
    figure(2);
elseif strcmp(Pr.trace,'yes')
    
    % TRACING
    
    if strcmp(Pr.cg,'yes');cg1=cell(1,Pr.nMode); end
    cp1=cell(1,Pr.nMode);  kz1=cp1;
    cp=repmat(Pr.w*2*pi,size(kz,1),1)./real(kz);
    for i=1:Pr.nMode
        k=zeros(1,length(Pr.w));
        [cpm,cpl]=min(cp);
        for j=1:length(cpl);  k(j)=kz(cpl(j),j);   end
        [~,iloc]=find(cpm==Inf);
        cpm(iloc)=-1;         k(iloc)=0;    cp1{1,i}=cpm;
        if strcmp(Pr.cg,'yes')
            for j=1:size(cp,2); cg1{1,i}(j)=cg(cpl(j),j); end
        end
        kz1{1,i}=k; for j=1:size(cp,2); cp(cpl(j),j)=Inf; end
    end
    cp1=cp1(1,1:min(Pr.nMode,size(cp1,2)));
    
    % PLOT
    
    c_l={'k','r','b','g','m'};for cli=1:20; c_l=[c_l c_l]; end
    figure(2)
    for i=1:size(cp1,2)
        w=Pr.w;
        c=cp1{1,i};
        [~,iloc]=find(c>0);
        c=c(iloc);w=w(iloc);
        if strcmp(Pr.eff,'no')
            plot(w,c,'-','Color',c_l{1,i},'LineWidth',1.5);
        elseif strcmp(Pr.eff,'yes')
            pT=plot(w,c,'-','Color','b','LineWidth',1.5);
        end
        hold on;
    end
    if strcmp(Pr.eff,'yes')
        [iloc]=find(cpE~=0);
        pE=plot(Pr.w(iloc),cpE(iloc),'ro');
        legend([pT,pE],'Theoretical','Effective', 'Location','southeast')
    end
    xlabel('Frequency (Hz)');ylabel('Phase Velocity (m/s)');
    ylim([Pr.cpMin Pr.cpMax]); box on;
    figure(3);
    for i=1:size(cp1,2)
        w=Pr.w;
        k=kz1{1,i};
        [~,iloc]=find(k~=0);k=k(iloc);w=w(iloc);
        plot(w,real(k),'-','Color',c_l{1,i},'LineWidth',1.5);
        hold on
    end
    xlabel('Frequency (Hz)');ylabel('Wavenumber (rad/m)'); box on;
    if strcmp(Pr.cg,'yes')
        figure(4)
        for i=1:size(cp1,2)
            w=Pr.w;
            c=cp1{1,i};   g=cg1{1,i};
            [~,ilocC]=find(c>0); [~,ilocG]=find(g>0); iloc=intersect(ilocC,ilocG);
            c=c(iloc);w=w(iloc); g=g(iloc);
            plot(w,g,'-','Color',c_l{1,i},'LineWidth',1.5);
            hold on
        end
        xlabel('Frequency (Hz)');ylabel('Group Velocity (m/s)'); box on;
    end
    figure(2);
end

%}

misfit = Misfit(cpE,c_obs);

end
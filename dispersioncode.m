%% Input 

sg=input('signal matrix = ');
fs = input('acquistion frequency = ');
disp_freq=input('dispersion maximum frequency = ');
max_vel = input('dispersion maximum phase velocity = ');
so=input('source offset = '); 
dx=input('receiver interval ='); 

sg = [sg;sg;sg;sg;sg]; % padding

sz=size(sg); M=sz(2);
x=so+((1:M)-1)*dx;

L =2^nextpow2(sz(1));
f = fs*linspace(0,1,L+1);
w=2*pi*f;
N=2000;

%% FFT
X=fft(sg,L);
Xp=X./abs(X);
kk=ceil((L+1)/fs*disp_freq);
V=zeros(kk,N);
c=zeros(1,N);

%% wavefield transform
for j=1:N
    c(j)=max_vel*j/N;
    for k=1:kk
        for g=1:M
            V(k,j)=V(k,j)+exp(1i*(w(k)/c(j))*x(g))*Xp(k,g);%.*dx;
        end
    end
end
ff=ceil((L+1)/(fs)*1);%round(1/df)+1;
v1=V(ff:kk,:);
Z=abs(v1')./max(max(abs(v1)));

%% dispersion plot

figure;imagesc(f(ff:kk),c,Z); colormap jet; shading interp; cb=colorbar; axis xy; xlim([1 disp_freq]); ylim([0 max_vel]);


%% Predominant mode determination

min_vel=input('minimum cutoff velocity to avoid lower velocity spurious modes (e.g. 50) = ');
min_vel_idx=find(abs(min_vel-c)==min(abs(min_vel-c)));

[zmax, idx]=max(Z(min_vel_idx:end,:));
cid=c(idx+min_vel_idx-1);
fd=f(ff:kk);
ints=10:1.5:100;  % provide the vector of interested frequency values for predominant mode determination
[~, fid] = arrayfun(@(x) find(abs(fd - x) == min(abs(fd - x))), ints);
f_obs=fd(fid); % predominant mode frequency values
c_obs=cid(fid); % % predominant mode phase velocity values


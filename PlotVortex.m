%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S1 = imread('June7/V1.bmp');
S1=im2double(S1);
S1=rgb2gray(S1(:,:,1:3));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dz = 0.001872;
dx = 0.1;
W = 40; % simulation window
theta = 40; % for the interference pattern
Input = 'What input? ';
Inp = input(Input);
display(Inp*dz/1000) %mm
Inp = sprintf('%d',Inp);
InpAmp = strcat(Inp,'.dat'); 
InpPhase = strcat(Inp,'phase','.dat');
Eamp = load(InpAmp);
sEamp = size(Eamp);
Ephase = load(InpPhase);
E1 = Eamp.*exp(j*Ephase);
figure(1); surf(Eamp./max(max(Eamp))); shading interp; view([0,90]); colormap gray; axis square; colorbar
figure(2);surf(Ephase); shading interp;view([0,90]);colormap gray;axis square; colorbar 
N = load('N5E.dat');
%figure(3);surf(N);shading interp;view([0,90]);axis square;colorbar
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%Interference with a plane wave
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nx = W/dx;
[X Y] = meshgrid(-W/2:dx:W/2,-W/2:dx:W/2);
X = X(1:Nx,1:Nx); 
Y = Y(1:Nx,1:Nx);
Wg = 2;
Plane = 1*exp(j*25*sin(theta).*X);
G0 = exp(-(X.^2+Y.^2)/Wg^2);
clear Y; 
Gauss = G0.*exp(j*30*sin(theta).*X);
Inter = abs(Plane+E1).^2; 
figure(4); surf(Inter./max(max(Inter))); 
shading interp; view([0,90]); axis tight; axis square;colormap gray;colorbar;
%annotation('ellipse',[0.5 0.5 0.2 0.2],'LineWidth',2,'Color','r')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%finding phase singularitis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MinPhase = min(min(abs(Ephase)))
[r c] = find(abs(Ephase)==MinPhase)
phaseE = (Ephase);
xp = linspace(-20,20,sEamp(2));
figure(5); plot(xp,Ephase(:,250),'Linewidth',2); set(gca,'fontsize',20); xlabel('X(\mum)');ylabel('Phase(rad)')
name=strcat(Inp,'.png')
saveas(gca,name,'png')




function [Hz,magsweepup,magsweepdown,j] = Hsweep(percK2,thetaoff)
% A color plot of the resultant magnetization if we sweep field
% The two inputs are:
% thetaoff: the offset angle of the easy axis
% percK2: K2's percentage of K1 

tic;
StartOfProg = toc;

% Initial conditions in time and magnetization
mz_0 = .9; mx_0 = 0; my_0 = .1;
% Magnitude of the magnetization vector
m = sqrt((mx_0)^2+(my_0)^2+(mz_0)^2);
% Make the magnetization vector unit sized
mx_0 = mx_0 / m; my_0 = my_0 / m; mz_0 = mz_0 / m; 

% Integration times
t_0 = 0; t_f = 2.5e-8;

% Mutual Demagnetization Tensor
X = 0; Y = 0; Z = 7;
% circle
dx=50; dy=50; dz=3;
% ellipse
% dx = 100; dy = 50; dz = 3;
NR = demagtensor(X,Y,Z,dx,dy,dz);
% For infinite lateral size use NR_xx = 0; NR_yy = 0; NR_zz = 0;

% Self Demagnetization Matrix
X = 0; Y = 0; Z = 0;
% circle
dx=50; dy=50; dz=3;
% ellipse
% dx = 100; dy = 50; dz = 3;
ND = demagtensor(X,Y,Z,dx,dy,dz);
% For infinite lateral size use ND_xx = 0; ND_yy = 0; ND_zz = 1;

% Define the normal LLG parameters for the problem
Gamma = 2.21e5; 		% Gyromagnetic ratio, 2.2127616 x 1e5 s^(-1).A^(-1).m
Alpha = 0.01; 			% Gilbert's damping factor, dimensionless
d = 3e-9;               % Free layer thickness, m
Ms = 650000; 			% Free layer saturation magnetization, A/m
K1 = 3e5;               % Free layer anisotropy constant, J/m^3
K2 = (percK2/100)*K1;   % Free layer anisotropy constant
mu0 = 4*pi*1e-7; 		% Magnetic permeability of vaccuum, H/m
e = 1.602e-19; 			% Electron charge, -1.602176 x 1e-19 C
hbar = 1.054e-34; 		% Reduced Planck constant, 1.054571 x 1e-34 J.s

% Other constants
MRs = 600000;                % Reference layer saturation magnetization, A/m
px = 0; py = 0; pz = 1;		 % Pinned layer initial state
p = [px; py; pz];
    
% Slonczewski's spin transfert torque coefficient Beta=P*j*Beta_div_P_and_j
Beta_div_P_and_j = (Gamma*hbar)/(2*d*mu0*Ms*e);

% Some useful parameters
H_K1 = 2*K1/mu0/Ms;                   % Free layer anisotropy field, A/m
H_K2 = 4*K2/mu0/Ms;                   % Free layer anisotropy field, A/m
norm = H_K1/sqrt(H_K1^2+H_K2^2);
H_K1 = norm*H_K1;
H_K2 = norm*H_K2;
H_K2 = 0;
Gamma_bar = Gamma/(1+Alpha^2);
Alpha_bar = Alpha*Gamma/(1+Alpha^2);

% Initialize our applied field vector
Hzmin = -1e6;%-0.6e6;
Hzmax = 1e6;%0.3e6;
nfield = 50;
Hz = linspace(Hzmin,Hzmax,nfield);

% Current Density Vector
jmin = -1e11;%-0.5e11;
jmax = 1.5e11;%1e11;
n = 50;
j = linspace(jmin,jmax,n); 

% Final magnetization
magsweepdown = zeros(n,nfield);
magsweepup = zeros(n,nfield);

% Offset for easy axis of free FM layer
thetaoffset = thetaoff*pi/180;

StartOfFor = toc;

for i=1:n

    % Sweep down
    IC = [mx_0 my_0 mz_0];
    for k=nfield:-1:1
        [T,Y] = ode45(@llg1,[t_0 t_f],IC);
        % Record final magnetization
        magsweepdown(i,k) = Y(end,3);
    end
    
    % Sweep up
    IC = [mx_0 my_0 -mz_0];
    for k=1:nfield
        options = odeset('RelTol',1e-4,'AbsTol',5e-6);
        [T,Y] = ode45(@llg1,[t_0 t_f],IC,options);
        % Record final magnetization
        magsweepup(i,k) = Y(end,3);
    end

end

EndOfFor = toc;
fprintf(1, 'Elapsed time in For loop: %f seconds\n', (EndOfFor-StartOfFor));

    function dy = llg1(t,y)

    % Time depdendent parameters
    Beta = Beta_div_P_and_j*Polarisation(y(1),y(2),y(3),px,py,pz)*j(i);
    Beta_bar1 = Alpha*Beta/(1+Alpha^2);
    Beta_bar2 = Beta/(1+Alpha^2);
    
     % Effective field is the sum of the anisotropy field (Hani), the self (HD)
     % and mutual (HR) demag. fields, and the external field (H)
     % Applied Field
     H = [0; 0; Hz(k)];
     % Anistropic field
     Hani = [H_K1*y(1)*sin(thetaoffset); 0; H_K1*y(3)*cos(thetaoffset)+H_K2*y(3)*(y(1)^2+y(2)^2)];
     % Self-demag field
     HD = -ND*y*Ms;
     % Mutual demag field
     HR = -NR*p*MRs;
     % Effective field
     Heff = H + Hani + HD + HR;

     % In fact, one has to integrate three differential equations,
     % one for each component of m
     % dm_x/dt=fx(m_x,m_y,m_z,t)
     dy = zeros(3,1);

     dy(1) = -Gamma_bar*(y(2)*Heff(3)-y(3)*Heff(2))-Alpha_bar*((y(1)*y(1)-1)*Heff(1)+...
         y(1)*y(2)*Heff(2)+y(1)*y(3)*Heff(3))+Beta_bar1*(y(2)*p(3) - y(3)*p(2))...
           -Beta_bar2*((y(1)*y(1)-1)*p(1)+y(1)*y(2)*p(2)+y(1)*y(3)*p(3));

     dy(2) = -Gamma_bar*(y(3)*Heff(1)-y(1)*Heff(3))-Alpha_bar*((y(2)*y(2)-1)*Heff(2)+...
         y(1)*y(2)*Heff(1)+y(2)*y(3)*Heff(3))+Beta_bar1*(y(3)*p(1)-y(1)*p(3))...
           -Beta_bar2*((y(2)*y(2)-1)*p(2)+y(1)*y(2)*p(1)+y(2)*y(3)*p(3));

     dy(3) = -Gamma_bar*(y(1)*Heff(2)-y(2)*Heff(1))-Alpha_bar*((y(3)*y(3)-1)*Heff(3)...
         +y(1)*y(3)*Heff(1)+y(2)*y(3)*Heff(2))+Beta_bar1*(y(1)*p(2)-y(2)*p(1))....
           -Beta_bar2*((y(3)*y(3)-1)*p(3)+y(1)*y(3)*p(1)+y(2)*y(3)*p(2));
    
    end

% Make up some names for writing out the files
a = strcat('magsweepup',num2str(percK2));
b = strcat('magsweepdown',num2str(percK2));
c = strcat('current',num2str(percK2));
d = strcat('field',num2str(percK2));

% Get rid of weird integration errors
magsweepdown(find(magsweepdown>1))=1;
magsweepup(find(magsweepup>1))=1;
magsweepdown(find(magsweepdown<-1))=-1;
magsweepup(find(magsweepup<-1))=-1;

% Write the files out
dlmwrite(a,magsweepup);
dlmwrite(b,magsweepdown);
dlmwrite(c,j);
dlmwrite(d,Hz);

% More concatenation
a = strcat('PAP',num2str(percK2),'.jpg');
b = strcat('APP',num2str(percK2),'.jpg');
c = strcat('sum',num2str(percK2),'.jpg');
d = strcat('diff',num2str(percK2),'.jpg');

% Plot the resultant color maps
% Sweep from Parallel to Anti-Parallel
hd = pcolor(Hz, j, magsweepdown);
colorbar;
title(['Parallel to Antiparallel Switching - K2 is ',num2str(percK2),' percent of K1']);
xlabel('Applied Field (A/m)');
ylabel('Current (A/m^2)');
saveas(hd,a) 

% Sweep from Anti-Parallel to Parallel
hu = pcolor(Hz, j, magsweepup);
colorbar;
title(['Antiparallel to Parallel Switching - K2 is ',num2str(percK2),' percent of K1']);
xlabel('Applied Field (A/m)');
ylabel('Current (A/m^2)');
saveas(hu,b) 

% The sum of the two - spinwave region
sum = magsweepup+magsweepdown;
hs = pcolor(Hz, j, sum);
colorbar;
title(['Sum of Switchings - K2 is ',num2str(percK2),' percent of K1']);
xlabel('Applied Field (A/m)');
ylabel('Current (A/m^2)');
saveas(hs,c) 

% The difference of the two - hysteretic part
diff = magsweepdown-magsweepup;
hs = pcolor(Hz, j, diff);
colorbar;
title(['Difference of Switchings - K2 is ',num2str(percK2),' percent of K1']);
xlabel('Applied Field (A/m)');
ylabel('Current (A/m^2)');
saveas(hs,d) 

EndOfProg = toc;
fprintf(1, 'Elapsed time in program loop: %f seconds\n', (EndOfProg-StartOfProg));

end

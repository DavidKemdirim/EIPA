% ELEC 4700 PA 5 Harmonic Wave Equation in 2D FD and Modes

close all
clear 
clc

global C

C.q_0 = 1.60217653e-19;             % electron charge
C.hb = 1.054571596e-34;             % Dirac constant
C.h = C.hb * 2 * pi;                % Planck constant
C.m_0 = 9.10938215e-31;             % electron mass
C.kb = 1.3806504e-23;               % Boltzmann constant
C.eps_0 = 8.854187817e-12;          % vacuum permittivity
C.mu_0 = 1.2566370614e-6;           % vacuum permeability
C.c = 299792458;                    % speed of light

% dE^2/dx^2 + dE^2/dy^2 = aE

nx = 50;
ny = 50;
G = sparse(nx*ny,nx*ny);
d = 1; % dx = dy
d2 = d^2; % dx^2 = dy^2

for i = 1:nx
    for j = 1:ny
        
        if i == j
            G(i*j,i*j) = 1;
        elseif i == 1
            G(i*j,i*j) = 0;
        elseif j == 1
            G(i*j,i*j) = 0;
        elseif i == nx*ny
            G(i*j,i*j) = 0;
        elseif j == nx*ny
            G(i*j,i*j) = 0;
        else 
%             G(i,:) = 0;
        
            G(i*j,i*j) = -4/d2;
            G(i*j,i*j-1) = 1/d2;
            G(i*j,i*j+1) = 1/d2;
            G(i*j-1,i*j) = 1/d2;
            G(i*j+1,i*j) = 1/d2;
        end                
        % [G(j-1,k) -2*G(j,k) + G(j+1,k)]/dx^2;
        % [G(j,k-1) -2*G(j,k) + G(j,k+1)]/dy^2;        
    end 
end

figure('name','Matrix')
figure(1)
spy(G)

modes = 20;
[E,D] = eigs(G, modes);

figure('name','Eigenvalues')
figure(2)
plot(diag(D),'*')

np = ceil(sqrt(modes));

figure('name','Modes')
V = zeros(nx,ny);
for k = 1:modes
   M = E(:,k);
   for i = 1:nx
       for j = 1:ny           
           n = j + (i-1)*ny;
           V(i,j) = M(n);           
       end
       subplot(np,np,k), surf(V)
       title(['EV = ', num2str(D(k,k))])
   end      
end

clear
clc
format shorte
close all

% Parameters for the Schnakenberg model
D_u = 0.05;
D_v = 1;
gamma = 100;
A = 0.1305;
B = 0.7695;

% Time parameters
T = 2;
dt = 1e-5;
Mt = T / dt;

% Set here the parameters of the square box domain and mesh :
LX=1;
LY=1;
NELX = 6; NELY = 6; P = 6;
dxe = LX/NELX;
dye = LY/NELY;
NEL = NELX*NELY;
NGLL = P+1; % number of GLL nodes per element
[iglob,x,y] = MeshBox(LX,LY,NELX,NELY,NGLL);
nglob = length(x);
[xgll,wgll,H] = GetGLL(NGLL);
Ht = H';
wgll2 = wgll * wgll' ;

% (x,y)->(xi,eta)
dx_dxi  = 0.5*dxe;
dy_deta = 0.5*dye;
jac = dx_dxi*dy_deta;

% global mass matrix, diagonal
M = zeros(nglob,1);		
for e=1:NEL
    ig = iglob(:,:,e);
    % Diagonal mass matrix
    M(ig) = M(ig) + wgll2 *jac;
end 

% local contributions for stiffness: set material properties here.
W = WMatrix_New(NELX, NELY, NGLL, wgll2, dxe, dye);
% assemble stiffness as a sparse matrix
D = assemble_K_matrix_local_2d(NELX, NELY, NGLL, dxe, dye, nglob, iglob,W);
M = diag(M);

% for solve system by C-N
A_left_u  = M + 0.5 * dt * D_u * D;
A_right_u = M - 0.5 * dt * D_u * D;
A_left_v  = M + 0.5 * dt * D_v * D;
A_right_v = M - 0.5 * dt * D_v * D;

% Initial conditions
u_initial = 1 + B + 1e-3 .* exp(-100 .* ((x - 1/3).^2 + (y - 1/2).^2) );
v_initial = (B / (A + B)^2) .* ones(size(u_initial));

U = u_initial; 
V = v_initial; 

for k = 1:Mt
    t = k * dt
    
    % reaction terms
    R_u = gamma * (A - U + (U.^2) .* V);
    R_v = gamma * (B - (U.^2) .* V);
    
    % right-hand side for u and v
    F_u = A_right_u * U + dt .* M * R_u;
    F_v = A_right_v * V + dt .* M * R_v;

    % Solve the systems
    U = A_left_u\ F_u;
    V = A_left_v\ F_v;
end

[X,Y] = meshgrid(0:0.01:LX);
Z1 = griddata(x,y,U,X,Y);
pcolor(X,Y,Z1), 
shading interp 

Z2 = griddata(x,y,V,X,Y);
pcolor(X,Y,Z2), 
shading interp 
%% Initialization
clear; close all; clc;

%% Defining the symbolic parameters for the Simulation
syms t
syms m l positive
syms x(t) y(t) theta(t)

dx = diff(     x, t );
dy = diff(     y, t );
 w = diff( theta, t );

% The center of mass position 
xc = x + l * cos( theta );
yc = y + l * sin( theta );
I = 1/12 * (2*m) * (2*l)^2;

dxc = diff(    xc, t );
dyc = diff(    yc, t );

% The kinetic energy (and thus, the Lagrangian)
L = 1/2 * (2*m) * (dxc^2 + dyc^2) + 1/2 * I * w^2;

% The Equations of motion
eq1 = simplify( diff( diff( L, dx ), t ) - diff( L,     x ) );
eq2 = simplify( diff( diff( L, dy ), t ) - diff( L,     y ) );
eq3 = simplify( diff( diff( L,  w ), t ) - diff( L, theta ) );

% The mass matrix 
M = sym( zeros( 3, 3) );
M( 1, 1 ) = 2 * m;
M( 2, 2 ) = 2 * m;
M( 1, 3 ) = -2*m*l*sin( theta );
M( 3, 1 ) = -2*m*l*sin( theta );
M( 2, 3 ) =  2*m*l*cos( theta );
M( 3, 2 ) =  2*m*l*cos( theta );
M( 3, 3 ) =  I + 2*m*l^2;

% Coriolis
C = sym( zeros( 3, 1 ) );
C( 1 ) = -2*m*l*cos( theta ) * w^2;
C( 2 ) = -2*m*l*sin( theta ) * w^2;

% Sanity Check
tmpv = [ dx; dy; w ];
tmp = simplify( tmpv.' * diff( M , t )  * tmpv - 2 * tmpv.' * C );

% The Force input is simply 
syms k b
syms r0 w0 F0

% x0 = r0 * cos( w0 * t );
% y0 = r0 * sin( w0 * t );

Fx = F0 * cos( w0 * t );
Fy = F0 * sin( w0 * t );

% dx0 = diff( x0, t );
% dy0 = diff( y0, t );

% ddx0 = diff( dx0, t );
% ddy0 = diff( dy0, t );

F = [ -k*x - b*dx + Fx; ...
      -k*y - b*dy + Fy; 0 ];

%% Substitute the value
M_mat = subs( M, { m,l, theta }, {1,1, 'theta' } );
C_mat = subs( C, { m,l, diff( theta, t ), theta }, {1,1, 'w', 'theta' } );
F_mat = subs( F, { k,b, F0, w0, diff( x, t ), diff( y, t ), x, y }, { 4, 4, 3, 2, 'dx', 'dy', 'x', 'y' } );

M_func = matlabFunction( M_mat );
C_func = matlabFunction( C_mat );
F_func = matlabFunction( F_mat );

dt = 1e-3;
T  = 20;
t_arr = 0:dt:T;
N = length( t_arr );

Ntrial = 120;
x0_arr =  6*rand( 6,    Ntrial );
x_arr  = zeros( 6, N, Ntrial );


for j = 1 : Ntrial
    
    x_curr = x0_arr( :, j );
    j
    for i = 1 : N
    
        x_arr( :, i, j ) = x_curr;
    
        t = t_arr( i );
    
        M_mat = M_func( x_curr( 3 ) );
        C_mat = C_func( x_curr( 3 ), x_curr( 6 ) );
        F_mat = F_func( t, x_curr( 4 ), x_curr( 5 ), x_curr( 1 ), x_curr( 2 ) );
    
        dx = zeros( 6, 1 );
    
        dx( 1:3 ) = x_curr( 4:6 );
        dx( 4:6 ) = inv( M_mat ) * ( -C_mat + F_mat );
    
        x_curr = x_curr + dx * dt;
        
    end

end

%%

f = figure( ); a = axes( 'parent', f );
hold on; axis equal
for i = 1 : Ntrial
    plot( a, 2*i + x_arr( 1, : , i ), x_arr( 2, : , i ), 'linewidth', 3 )
end
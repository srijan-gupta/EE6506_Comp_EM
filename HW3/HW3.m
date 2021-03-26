tic
% 1. Defining the parameters
lambda = 1;
k0 = 2*pi/lambda;
a = 3*lambda;               %edge length of pentagon
da = lambda/15;           %discretization length, please ensure that rem(a,da) = 0
theta_i = 45*(pi/180); %angle of incidence
a_ff = 4*lambda;          %far field radius
n_ff = 180;                    %discretizations for far field
tolabs = 1e-9;               %absolute tolerance in integral
tolrel = 1e-6;                %relative tolerance in integral

% 2. Defining the variables 
% 2.1. Related to the pentagon
[n_e, test_pt, strt_pt, normals] = get_shape_coords(a,da, 4);  
figure;
scatter(test_pt(1,:), test_pt(2,:))
hold on
scatter(strt_pt(1,:), strt_pt(2,:))

params = [k0, da, theta_i, a_ff, n_ff, tolabs, tolrel, n_e];

% 2.2. Far field points
th_ff = 0:360/n_ff:(360-360/n_ff);
ff_pt = a_ff*[cosd(th_ff); sind(th_ff)];

% 3. Calculating RCS
% For PEC, please set eps_r as Inf
% 3.1. For PEC
eps_r = Inf;
fields_bndry = solve_on_boundary(eps_r, test_pt, strt_pt, normals, params);
RCS_PEC = get_RCS(eps_r, fields_bndry, ff_pt, strt_pt, normals, params);

%3.2. For carbon fibre
eps_r = 2;
fields_bndry = solve_on_boundary(eps_r, test_pt, strt_pt, normals, params);
RCS_CF = get_RCS(eps_r, fields_bndry, ff_pt, strt_pt, normals, params);

toc

figure;
polarplot(th_ff*pi/180, RCS_PEC);
title('RCS for PEC')
figure;
polarplot(th_ff*pi/180, RCS_CF);
title('RCS for carbon fibre')



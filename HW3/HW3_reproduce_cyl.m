tic
% 1. Defining the parameters
lambda = 1;
k0 = 2*pi/lambda;
n_sides = 90;
a = 0.0698*lambda;               %edge length of pentagon
da = 0.0698*lambda;           %discretization length, please ensure that rem(a,da) = 0
theta_i = 0*(pi/180); %angle of incidence
a_ff = 4*lambda;          %far field radius
n_ff = 180;                    %discretizations for far field
tolabs = 1e-9;               %absolute tolerance in integral
tolrel = 1e-6;                %relative tolerance in integral

% 2. Defining the variables 
% 2.1. Related to the pentagon
[n_e, test_pt, strt_pt, normals] = get_shape_coords(a, da, n_sides);  
figure;
scatter(test_pt(1,:), test_pt(2,:))
hold on
scatter(strt_pt(1,:), strt_pt(2,:))

params = [k0, da, theta_i, a_ff, n_ff, tolabs, tolrel, n_e];

% 2.2. Far field points
th_ff = 0:360/n_ff:(360-360/n_ff);
ff_pt = a_ff*[cosd(th_ff); sind(th_ff)];

% 3. Calculating RCS
eps_r = 4;
%fields_bndry = solve_on_boundary(eps_r, test_pt, strt_pt, normals, params);
n = length(test_pt);
strt_pt(:,end+1) = strt_pt(:,1); %for ease while calculating Dr for the last starting point
A = zeros(2*n, 2*n);
b = zeros(2*n,1);
for i = 1:n
    for j = 1:n
        Dr = strt_pt(:,j+1) - strt_pt(:,j);
        n_hat_id = fix((j-1)/n_e) + 1;
        n_hat = normals(:, n_hat_id);
        A(i, j)     = integral(@(d)green2d(k0, test_pt(:,i), strt_pt(:,j), Dr, d),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);            %The g1 term
        A(i+n, j)   = integral(@(d)green2d(sqrt(eps_r)*k0, test_pt(:,i), strt_pt(:,j), Dr, d),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);%The g2 term
        A(i, j+n)   = -1*integral(@(d)gradgreen2d_dot_n(k0, test_pt(:,i), strt_pt(:,j), Dr, d, n_hat),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);            %The grad(g1).n term
        A(i+n, j+n) = -1*integral(@(d)gradgreen2d_dot_n(sqrt(eps_r)*k0, test_pt(:,i), strt_pt(:,j), Dr, d, n_hat),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);%The grad(g2).n term
    end
    b(i) = inc_field(test_pt(:,i), theta_i, k0);
end
A = da*A; %scaling by the length of the segment since the integral iterated over d from 0 to 1
x = A\b;
RCS = get_RCS(eps_r, fields_bndry, ff_pt, strt_pt, normals, params);

toc

figure;
polarplot(th_ff*pi/180, sqrt(RCS));
title(['sqrt(RCS) for n_sides = ',num2str(n_sides),', eps_r = ',num2str(eps_r)])



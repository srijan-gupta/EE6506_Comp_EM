tic
% 1. Defining the parameters
global k0 da theta_i a_ff n_ff tolabs tolrel n_e
lambda = 1;
k0 = 2*pi/lambda;
a = 3*lambda; %edge length of pentagon
da = lambda/15; %discretization length, please ensure that rem(a,da) = 0
theta_i = 45*(pi/180); %angle of incidence
a_ff = 4*lambda; %far field radius
n_ff = 180; %discretizations for far field
tolabs = 1e-9; %absolute tolerance in integral
tolrel = 1e-6; %relative tolerance in integral

% 2. Defining the variables 
% 2.1. Related to the pentagon
c = a/(2*sind(36)); %length of line segment connecting centre to vertex
vrtx_A = [c*cosd(18); c*sind(18)];
dir_AB = [-cosd(36); sind(36)];
rotat72 = [cosd(72) -sind(72); sind(72) cosd(72)]; 
n_e = round(a/da);
test_pt = zeros(2,5*n_e);
strt_pt = zeros(2,5*n_e);
test_pt(:,1:n_e) = vrtx_A + [da/2:da:a].*dir_AB;
strt_pt(:,1:n_e) = vrtx_A + [0:da:a-da].*dir_AB;
for i = 1+n_e:n_e:5*n_e
    test_pt(:,i:(i+n_e-1)) = rotat72*test_pt(:,i-n_e:i-1);
    strt_pt(:,i:(i+n_e-1)) = rotat72*strt_pt(:,i-n_e:i-1);
end
% figure;
% scatter(test_pt(1,:), test_pt(2,:))
% hold on
% scatter(strt_pt(1,:), strt_pt(2,:))
th_nm_deg = 54:72:(54+4*72);
normals = [cosd(th_nm_deg);sind(th_nm_deg)];

% 2.2. Far field points
th_ff = 0:360/n_ff:(360-360/n_ff);
ff_pt = a_ff*[cosd(th_ff); sind(th_ff)];

% 3. Calculating RCS
% For PEC, please set eps_r as Inf
% 3.1. For PEC
eps_r = Inf;
fields_bndry = solve_on_boundary(eps_r, test_pt, strt_pt, normals);
RCS_PEC = get_RCS(eps_r, fields_bndry, ff_pt, strt_pt, normals);

%3.2. For carbon fibre
eps_r = 2;
fields_bndry = solve_on_boundary(eps_r, test_pt, strt_pt, normals);
RCS_CF = get_RCS(eps_r, fields_bndry, ff_pt, strt_pt, normals);

toc

figure;
polarplot(th_ff*pi/180, RCS_PEC);
title('RCS for PEC')
figure;
polarplot(th_ff*pi/180, RCS_CF);
title('RCS for carbon fibre')

function g = green2d(k, r_p, r0, Dr, d)
% r_p = [r_p_x r_p_y] is the primed coordinate
% r0 = [r0_x r0_y] is the starting point
% Dr = [Dr_x Dr_y] is the total change in r0
% d \in [0,1) gives r = r0 + d*Dr

% Green's function, g(r,r_p) = -j/4*H_0^(2)(k|r' - r|)

R_x = (r0(1) + d.*Dr(1)) - r_p(1);
R_y = (r0(2) + d.*Dr(2)) - r_p(2);
norm_R = sqrt( R_x.^2 + R_y.^2 );
g = -1j/4*besselh(0, 2, k.*norm_R);
end

function grad_g_dot_n = gradgreen2d_dot_n(k, r_p, r0, Dr, d, n_hat)
% n = [nx ny] is the direction of the normal
% rest of the variables are same as in green2d

% grad(g).n = jk/4 H_1^(2)(k|r'-r|)hat(r-r').n

R_x = (r0(1) + d.*Dr(1)) - r_p(1);
R_y = (r0(2) + d.*Dr(2)) - r_p(2);
norm_R = sqrt( R_x.^2 + R_y.^2 );
hat_R_dot_n = (n_hat(1).*R_x + n_hat(2).*R_y)./norm_R;

grad_g_dot_n = 1j*k/4*besselh(1, 2, k.*norm_R).*hat_R_dot_n;
end

function phi_i = inc_field(r_p, theta_i, k0)
phi_i = exp( -1j*k0*( cos(theta_i)*r_p(1) + sin(theta_i)*r_p(2) ) );
end

function x = solve_on_boundary(eps_r, test_pt, strt_pt, normals)
global k0 da theta_i tolabs tolrel n_e
n = length(test_pt);
strt_pt(:,end+1) = strt_pt(:,1); %for ease while calculating Dr for the last starting point
if eps_r == Inf %implies PEC
    A = zeros(n, n);
    b = zeros(n,1);
    for i = 1:n
        for j = 1:n
            Dr = strt_pt(:,j+1) - strt_pt(:,j);
            A(i, j) = integral(@(d)green2d(k0, test_pt(:,i), strt_pt(:,j), Dr, d),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);
        end
        b(i) = inc_field(test_pt(:,i), theta_i, k0);
    end
    A = da*A; %scaling by the length of the segment since the integral iterated over d from 0 to 1
    x = A\b;
else
    A = zeros(2*n, 2*n);
    b = zeros(2*n,1);
    for i = 1:n
        for j = 1:n
            Dr = strt_pt(:,j+1) - strt_pt(:,j);
            n_hat_id = fix((j-1)/n_e) + 1;
            n_hat = normals(:, n_hat_id);
            A(i, j)     = integral(@(d)green2d(k0, test_pt(:,i), strt_pt(:,j), Dr, d),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);            %The g1 term
            A(i+n, j)   = integral(@(d)green2d(sqrt(eps_r)*k0, test_pt(:,i), strt_pt(:,j), Dr, d),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);%The g2 term
            A(i, j+n)   = integral(@(d)gradgreen2d_dot_n(k0, test_pt(:,i), strt_pt(:,j), Dr, d, n_hat),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);            %The grad(g1).n term
            A(i+n, j+n) = integral(@(d)gradgreen2d_dot_n(sqrt(eps_r)*k0, test_pt(:,i), strt_pt(:,j), Dr, d, n_hat),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);%The grad(g2).n term
        end
        b(i) = inc_field(test_pt(:,i), theta_i, k0);
    end
    A = da*A; %scaling by the length of the segment since the integral iterated over d from 0 to 1
    x = A\b;
end
end

function RCS = get_RCS(eps_r, fields_bndry, ff_pt, strt_pt, normals)
global k0 da a_ff n_ff tolabs tolrel n_e
n = length(strt_pt);
strt_pt(:,end+1) = strt_pt(:,1); %for ease while calculating Dr for the last starting point

% phi(p) = phi_inc(p) - oint[g1(r,p) grad(phi).n - grad(g1(r,p).n phi]dr
scat_field_ff = zeros(1,n_ff);
if eps_r == Inf %implies PEC
    for i = 1:n_ff
        for j = 1:n
            Dr = strt_pt(:,j+1) - strt_pt(:,j);
            integral_g1 = da*integral(@(d)green2d(k0, ff_pt(:,i), strt_pt(:,j), Dr, d),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);
            scat_field_ff(i) = scat_field_ff(i) - fields_bndry(j)*integral_g1;
        end
    end
else
    for i = 1:n_ff
        for j = 1:n
            Dr = strt_pt(:,j+1) - strt_pt(:,j);
            n_hat_id = fix((j-1)/n_e) + 1;
            n_hat = normals(:, n_hat_id);
            integral_g1 = da*integral(@(d)green2d(k0, ff_pt(:,i), strt_pt(:,j), Dr, d),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);
            integral_grad_g1_dot_n = da*integral(@(d)gradgreen2d_dot_n(k0, ff_pt(:,i), strt_pt(:,j), Dr, d, n_hat),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);
            scat_field_ff(i) = scat_field_ff(i) - (fields_bndry(j)*integral_g1 - fields_bndry(j+n)*integral_grad_g1_dot_n);
        end
    end
end

RCS = 2*pi*a_ff*abs(scat_field_ff).^2; % (the magnitude of the incident field is just 1)
end


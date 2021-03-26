% 1. Defining the parameters
global k0 eps_r tolabs tolrel n_e
lambda = 1;
k0 = 2*pi/lambda;
a = 3*lambda; %edge length of pentagon
da = lambda/15; %discretization length, please ensure that rem(a,da) = 0
% For PEC, please give eps_r as -1
eps_r = 2; %relative permittivity
theta_i = 45*(pi/180); %angle of incidence
a_obs = 5*lambda; %observation radius
tolabs = 1e-9; %absolute tolerance in integral
tolrel = 1e-6; %relative tolerance in integral

% 2. Defining the variables related to the pentagon
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


% 3. Function for calculating \phi_i(r')
phi_i = @(r_p, theta_i, k0) exp( -1j*k0( cos(theta_i)*r_p(1) + sin(theta_i)*r_p(2) ) );

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
hat_R_dot_n = (n_hat(1).*r_x + n_hat(2).*r_y)./norm_R;

grad_g_dot_n = 1j*k/4*besselh(1, 2, k.*norm_R).*hat_R_dot_n;
end

function x = solve_on_boundary(test_pt, strt_pt, normals)
global k0 eps_r tolabs tolrel n_e
n = length(test_pt);
strt_pt(:,end+1) = strt_pt(:,1); %for ease while calculating Dr for the last starting point
if eps_r == -1 %implies PEC
    A = zeros(n, n);
    b = zeros(1,n);
    for i = 1:n
        for j = 1:n
            Dr = strt_pt(:,j+1) - strt_pt(:,j);
            A(i, j) = integral(@(d)green2d(k0, test_pt(:,i), strt_pt(:,j), Dr, d),0.0,1,'AbsTol',tolabs,'RelTol',tolrel);
        end
        b(i) = phi_i(test_pt(:,i), theta_i, k0);
    end
    x = A\b;
else
    A = zeros(2*n, 2*n);
    b = zeros(1,2*n);
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
        b(i) = phi_i(test_pt(:,i), theta_i, k0);
    end
    x = A\b;
end
end
    
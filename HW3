% 1. Defining the parameters
lambda = 1;
k0 = 2*pi/lambda;
a = 3*lambda; %edge length of pentagon
da = lambda/15; %discretization length
theta_i = 45*(pi/180); %angle of incidence

% 2. Defining the testing sites on the pentagon
c = a/(2*sind(36)); %length of line segment connecting centre to vertex
vrtx_A = [c*cosd(18); c*sind(18)];
dir_AB = [-cosd(36); sind(36)];
rotat72 = [cosd(72) -sind(72); sind(72) cosd(72)];
tst_pt_AB = vrtx_A + [da/2:da:a].*dir_AB;
tst_pt_BC = rotat72*tst_pt_AB;
tst_pt_CD = rotat72*tst_pt_BC;
tst_pt_DE = rotat72*tst_pt_CD;
tst_pt_EA = rotat72*tst_pt_DE;
test_pt = [tst_pt_AB, tst_pt_BC, tst_pt_CD, tst_pt_DE, tst_pt_EA];
% scatter(test_pt(1,:), test_pt(2,:))

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

function grad_g_dot_n = gradgreen2d_dot_n(k, r_p, r0, Dr, d, n)
% n = [nx ny] is the direction of the normal
% rest of the variables are same as in green2d

% grad(g).n = jk/4 H_1^(2)(k|r'-r|)hat(r-r').n

R_x = (r0(1) + d.*Dr(1)) - r_p(1);
R_y = (r0(2) + d.*Dr(2)) - r_p(2);
norm_R = sqrt( R_x.^2 + R_y.^2 );
hat_R_dot_n = (n(1).*r_x + n(2).*r_y)./norm_R;

grad_g_dot_n = 1j*k/4*besselh(1, 2, k.*norm_R).*hat_R_dot_n;
end

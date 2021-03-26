function phi_i = inc_field(r_p, theta_i, k0)
        phi_i = exp( -1j*k0*( cos(theta_i)*r_p(1) + sin(theta_i)*r_p(2) ) );
end
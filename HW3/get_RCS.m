function RCS = get_RCS(eps_r, fields_bndry, ff_pt, strt_pt, normals, params)

        [k0, da, ~, a_ff, n_ff, tolabs, tolrel, n_e] = feval (@(x) x{:} , num2cell(params));        
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

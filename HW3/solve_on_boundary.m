function x = solve_on_boundary(eps_r, test_pt, strt_pt, params)

        [k0, da, theta_i, ~, ~, tolabs, tolrel] = feval (@(x) x{:} , num2cell(params));         
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
                    t_hat = Dr/norm(Dr);
                    n_hat = [-t_hat(2), t_hat(1)];
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
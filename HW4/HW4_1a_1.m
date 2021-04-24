tic

n = 7;
seq =1;
eps_n = 2.4^2; %relative permitivitty
air_thickness = 1;
ratio = ((sqrt(5) + 1)/2);


k_vec = 0:pi/1000:pi;
len_vec = length(k_vec);
tou_arr = zeros(1,len_vec);
ref_arr = zeros(1,len_vec);

eps_arr = get_multilayer_eps(seq, n, eps_n);
wid_arr = get_width(eps_arr, air_thickness, ratio);
len = size(eps_arr');

tou = @(eps1, eps2) 2*eps1 ./ (eps1 + eps2); 
ref = @(eps1, eps2)  (eps1 - eps2) ./ (eps1 + eps2); 

I_mat  =@(t,r) (1/t) .* [1 r; r 1;];
delta = @(eps,wid, k) k*sqrt(eps)*wid;
P_mat = @(delta) [exp(1i*delta) 0; 0 exp(-1i*delta)];

eps_arr(end+1) = 1;

for k_id = 1:len_vec

        T_mat = I_mat( tou(1, eps_arr(1)), ref(1, eps_arr(1)) );
        for i = 1:len
              T_mat = T_mat * P_mat( delta(eps_arr(i), wid_arr(i), k_vec(k_id)) ) * I_mat( tou(eps_arr(i), eps_arr(i+1)), ref(eps_arr(i), eps_arr(i+1)) );     
        end

        net_tou = 1 / abs(T_mat(1,1));
        net_ref = abs( T_mat(2,1) / T_mat(1,1) );

        tou_arr(k_id) = net_tou;
        ref_arr(k_id) = net_ref;

end

toc

figure;
hold on;
plot(k_vec, tou_arr);
plot(k_vec, ref_arr);
xticks(k_vec(1:100:len_vec))
xticklabels(strcat(string(k_vec(1:100:len_vec)./pi), '\pi'))
xlabel('k')
hold off;
legend("transmission","reflection");


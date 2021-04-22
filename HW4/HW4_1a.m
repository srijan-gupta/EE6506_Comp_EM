lambda = 0.1;
global k0;
k0 = 2*pi / lambda;
n = 7;
seq =1;
eps_n = 3.5;
air_thickness = 1.5;

eps_arr = get_multilayer_eps(seq, n, eps_n);
wid_arr = get_width(eps_arr, air_thickness);
len = size(eps_arr');

tou = @(eps1, eps2) 2*eps1 ./ (eps1 + eps2); 
ref = @(eps1, eps2)  (eps1 - eps2) ./ (eps1 + eps2); 

I_mat  =@(t,r) (1/t) .* [1 r; r 1;];
delta = @(eps,wid) 2*pi*k0*eps*wid;
P_mat = @(delta) [exp(1i*delta) 0; 0 exp(-1i*delta)];

T_mat = I_mat( tou(1, eps_arr(1)), ref(1, eps_arr(1)) );
eps_arr(end+1) = 1;
for i = 1:len
      T_mat = T_mat * P_mat( delta(eps_arr(i), wid_arr(i)) ) * I_mat( tou(eps_arr(i), eps_arr(i+1)), ref(eps_arr(i), eps_arr(i+1)) );     
end

net_tou = 1 / abs(T_mat(1,1));
net_ref = abs( T_mat(2,1) / T_mat(1,1) );
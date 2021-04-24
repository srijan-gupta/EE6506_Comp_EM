tic

n_layers = 3;
seq =1;
eps_r = 1.4^2; %relative permitivitty
n = eps_r^0.5; %refractive index
air_thickness = 1;
ratio = ((sqrt(5) + 1)/2);

num_pts = 1000;

k_vec = 0:pi/100:pi;
len_vec = length(k_vec);
tau_arr = zeros(1,len_vec);
ref_arr = zeros(1,len_vec);

global wid_arr DL

n_obj_arr = [1 (get_multilayer_eps(seq, n_layers, eps_r)).^0.5 1];
wid_arr = get_width(n_obj_arr, air_thickness, ratio);
wid_arr([1 end]) = 10;
len = size(n_obj_arr');

DL = sum(wid_arr)/(num_pts-1);
n_node_arr = zeros(1,num_pts);

id_obj = 1;
for i = 1:num_pts
    if (i-1)*DL > sum(wid_arr(1:id_obj))
        id_obj = id_obj+1;
    end
    n_node_arr(i) = n_obj_arr(id_obj);
end

DLsq = DL^2;
intN1N1 = @(x1, x2) (x2^3-x1^3)/(3*DLsq) - (x2^2-x1^2)/DL + (x2-x1);
intN1N2 = @(x1, x2) -(x2^3-x1^3)/(3*DLsq) + (x2^2-x1^2)/(2*DL);
intN2N2 = @(x1, x2) (x2^3-x1^3)/(3*DLsq);

for k_id = 1:len_vec
    
    k = k_vec(k_id);
    Uin = @(x) exp(-1j*k.*x);
    
    alpha_in = -1j*k;
    alpha = @(n) 1j*k*n;
    
    A = zeros(num_pts);
    b = zeros(num_pts, 1);
    
    %assuming DL<wid_arr(1)
    id_obj = 1;
    
    n = n_node_arr(1);
    A(1,1) = 1/DL + (k*n)^2*intN1N1(0,DL) + alpha(n);
    A(1,2) = -1/DL + (k*n)^2*intN1N2(0,DL);
    b(1) = -(alpha_in - alpha(n));
    
    f = 1;
    for i = 2:(num_pts-1)
        n_im1 = n_node_arr(i-1);
        n_i   = n_node_arr(i);
        n_ip1 = n_node_arr(i+1);
        
        %A(i,i-1)
        A(i,i-1) = -1/DL + (k*n_im1)^2*intN1N2(0,f*DL) + (k*n_i)^2*intN1N2(f*DL,DL);
        
        %A(i,i):
        int1 = (k*n_im1)^2*intN2N2(0,f*DL) + (k*n_i)^2*intN2N2(f*DL,DL);
        
        if (i-1)*DL > sum(wid_arr(1:id_obj))
            id_obj = id_obj+1;
        end
        f = k1k2split(i, n_i, n_ip1, id_obj);
        
        int2 = (k*n_i)^2*intN1N1(0,f*DL) + (k*n_ip1)^2*intN2N2(f*DL,DL);
        A(i,i) = 2/DL + int1 + int2;
        
        %A(i,i+1)
        A(i,i+1) = -1/DL + (k*n_i)^2*intN1N2(0,f*DL) + (k*n_ip1)^2*intN1N2(f*DL,DL);
    
    end
    
    n = n_node_arr(end);
    A(num_pts,num_pts-1) = -1/DL + (k*n)^2*intN1N2(0,DL);
    A(num_pts,num_pts) = 1/DL + (k*n)^2*intN2N2(0,DL) - alpha(n);
    b(end) = -(alpha_in - alpha(n));
    
    U = A\b;
    
    ref_arr(k_id) = abs(U(1)-Uin(0));
    tau_arr(k_id) = abs(U(end)-Uin(sum(wid_arr)));
end

toc

figure;
hold on;
plot(k_vec, tau_arr);
plot(k_vec, ref_arr);
xticks(k_vec(1:fix(len_vec/10):len_vec))
xticklabels(strcat(string(k_vec(1:100:len_vec)./pi), '\pi'))
xlabel('k')
hold off;
legend("transmission","reflection");

x_arr = 0:DL:sum(wid_arr);
Uin_arr = Uin(x_arr);
figure
hold on
plot(x_arr, abs(U))
plot(x_arr, abs(Uin_arr))
plot(x_arr, abs(U'-Uin_arr))
plot(x_arr, real(U))
plot(x_arr, real(Uin_arr))
plot(x_arr, real(U'-Uin_arr))
plot(x_arr, n_node_arr)
legend('abs(U)','abs(Uin)','abs(Us)','Re(U)','Re(Uin)','Re(Us)', 'refr. index')

function frac_k1 = k1k2split(i, n1, n2, id_obj)
    global wid_arr DL
    if n1 == n2
        frac_k1 = 1;
    else
        wid_1 = sum(wid_arr(1:id_obj));
        wid_2 = DL*i;
        excess = wid_2 - wid_1;
        if excess>DL
            disp("excess width is more than permissible")
        else
            frac_k1 = excess/DL;
        end
    end
end
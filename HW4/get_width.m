function [wid_arr] = get_width(eps_arr,wid)
        %generates array of widths of each layer
        
         wid_a = wid;
         wid_n = wid * ((sqrt(5) + 1)/2);
         
         wid_arr = zeros(size(eps_arr));
         for i = 1 : size(eps_arr')
                 if eps_arr(i) == 1 
                         wid_arr(i) = wid_a;
                 else
                         wid_arr(i) = wid_n;
                 end
         end        
         
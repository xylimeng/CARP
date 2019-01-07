% function to generate the true 3D images 
% f1: f1 in Peihua Qiu's paper
% f2: f2 in Peihua Qiu's paper
% f3: 3D phantom 

% Input: n - length of each dimensiona; 
%        sigma_noise - standard deviation of Gaussian noise; 
%        type - type of simulated data, could be f1, f2 or f3

function [true, obs] = sim_data(n, sigma_noise, type)

if type ~= 3
    
n1 = n;
n2=  n;
n3 = n; 
z = zeros([1,n1*n2*n3]);

if (type == 1)
for (i = 1:n1) 
    for (j = 1:n2) 
        for (k = 1:n3)
            z((i-1)*n2*n3+(j-1)*n3+k) = -(i/n1-0.5)^2 -(j/n2-0.5)^2 -(k/n3-0.5)^2;
        if (  ((abs(i/n1-0.5)<=0.25) & (abs(j/n2-0.5)<=0.25) & (abs(k/n3-0.5)<=0.25)) )
            z((i-1)*n2*n3+(j-1)*n3+k) = -(i/n1-0.5)^2 -(j/n2-0.5)^2 -(k/n3-0.5)^2 + 1;
        end
        if ( (((i/n1-0.5)^2+(j/n2-0.5)^2 <=0.15^2) & (abs(k/n3-0.5)<=0.35))  )
            z((i-1)*n2*n3+(j-1)*n3+k) = -(i/n1-0.5)^2 -(j/n2-0.5)^2 -(k/n3-0.5)^2 + 1;
        end
        end
    end
end      
end

if (type == 2)
for (i = 1:n1) 
    for (j = 1:n2) 
        for (k = 1:n3)
            z((i-1)*n2*n3+(j-1)*n3+k) = 1/4 * sin(2*pi*(i/n1 + j/n2 + k/n3) + 1) + 0.25;
        if (  (((i/n1-0.5)^2 + (j/n2  - 0.5)^2 <= 1/4 * (k/n3 - 0.5)^2)) & (k/n3 >= 0.2) & (k/n3 <= 0.5) )
            z((i-1)*n2*n3+(j-1)*n3+k) = 1/4 * sin(2*pi*(i/n1 + j/n2 + k/n3) + 1) + 1.25; 
        end
        if ( (((i/n1-0.5)^2+(j/n2-0.5)^2 + (k/n3 - 0.5)^2) <= 0.4^2) & ((i/n1-0.5)^2+(j/n2-0.5)^2 + (k/n3 - 0.5)^2 > 0.2^2) & (k/n3 < 0.45)  )
            z((i-1)*n2*n3+(j-1)*n3+k) = 1/4 * sin(2*pi*(i/n1 + j/n2 + k/n3) + 1) + 1.25; 
        end
        end
    end
end      
end    


    g = z + randn(size(z)) * sigma_noise;   
     
      true = zeros([n1,n2,n3]); 
      obs = true;       
      for (i = 1:n1) 
          for (j = 1:n2)
              for (k = 1:n3)
                true(i,j,k) = z((i-1)*n2*n3+(j-1)*n3+k);
                obs(i,j,k)  = g((i-1)*n2*n3+(j-1)*n3+k);
              end
          end
      end
end

if type == 3
    true = phantom3d(n); 
    obs = true + randn([n,n,n]) .* sigma_noise;
end

end

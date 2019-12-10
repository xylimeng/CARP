function output = my_idwt(c, LoR, HiR) 

% input c 
atom_R = [LoR; HiR]; 
c3 = c(1:2);

cA2 = atom_R * [c3; c3]; % recover the approximation coefficients
c2 = reshape([cA2'; c(3:4)'], [4, 1]); 

W2_R = [circshift(atom_R, 2, 2); atom_R]; 
cA1 = W2_R * c2; 
c1 = reshape([cA1'; c(5:8)'], [8, 1]); 

W_R = zeros(8); 
W_R(1:2, [7, 8, 1, 2]) = atom_R; 
W_R(3:4, 1:4) = atom_R; 
W_R(5:6, 3:6) = atom_R; 
W_R(7:8, 5:8) = atom_R; 

output = W_R * c1; 

end 
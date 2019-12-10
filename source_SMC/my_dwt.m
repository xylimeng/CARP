function c = my_dwt(obs_vec, Lo, Hi) 

atom = [Lo; Hi]; 

W = zeros(8); 
W(1:2, 1:4) = atom; 
W(3:4, 3:6) = atom; 
W(5:6, 5:8) = atom; 
W(7:8, [7, 8, 1, 2]) = atom; 

c1 = W * obs_vec'; 
% apprxomation, detail, approximation, detail, ...
cA1 = c1(1:2:8); 
cD1 = c1(2:2:8); 

W2 = [atom; circshift(atom, 2, 2)];  
c2 = W2 * cA1;
cA2 = c2(1:2:4); 
cD2 = c2(2:2:4); 

% last level: per 
W3 = atom; 
c3 = W3 * [cA2;cA2]; 
c = [c3; cD2; cD1]; 
end 
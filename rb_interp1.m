function z = rb_interp1(y,x, energy_grid)

[C, ia, ic] = unique(y); 
z = interp1(y(ia), x(ia), energy_grid); 
end 
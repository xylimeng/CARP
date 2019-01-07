% 'path' is the absolute path of the directory containing the header 'armadillo'  
%  For example, it could be '/usr/local/include' in macOS

function mex_this(path)

if nargin == 0
    path = '/usr/local/include'; % default 
end 

ipath = ['-I' path];

mex('-v', ipath, '-larmadillo', '-lgfortran', 'treeLikelihood.cpp', 'tree_class.cpp','helper.cpp'); 
mex('-v', ipath, '-larmadillo', '-lgfortran', 'treeFit.cpp', 'tree_class.cpp','helper.cpp'); 
% mex -I/usr/local/include -larmadillo -lgfortran treeLikelihood.cpp tree_class.cpp helper.cpp -v
% mex -I/usr/local/include -larmadillo -lgfortran treeFit.cpp tree_class.cpp helper.cpp -v

end 

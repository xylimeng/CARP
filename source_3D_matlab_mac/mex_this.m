% mex -larmadillo -lgfortran treeLikelihood.cpp tree_class.cpp helper.cpp -v
% mex -larmadillo -lgfortran treeFit.cpp tree_class.cpp helper.cpp -v
% mex -larmadillo -lgfortran treeLikelihood_2.cpp tree_class.cpp helper.cpp -v
% mex -larmadillo -lgfortran treeFit_2.cpp tree_class.cpp helper.cpp -v
% mex -larmadillo -lgfortran rank2tube.cpp tree_class.cpp helper.cpp -v
% mex -larmadillo -lgfortran Intermediate.cpp tree_class.cpp helper.cpp -v
% mex -larmadillo -lgfortran log_expey.cpp tree_class.cpp helper.cpp -v

mex -larmadillo treeLikelihood.cpp tree_class.cpp helper.cpp -v
mex -larmadillo treeFit.cpp tree_class.cpp helper.cpp -v
mex -larmadillo treeLikelihood_2.cpp tree_class.cpp helper.cpp -v
mex -larmadillo treeFit_2.cpp tree_class.cpp helper.cpp -v
mex -larmadillo DrawPosition.cpp tree_class.cpp helper.cpp -v
mex -larmadillo MAP.cpp tree_class.cpp helper.cpp -v
mex -larmadillo MAPpar.cpp tree_class.cpp helper.cpp -v
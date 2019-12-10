function output = re_assign(input, flag)
output = input;
nflag = length(flag);

if min(size(input)) == 1
    for i = 1:nflag
        output(i) = input(flag(i));
    end
else
    for i = 1:nflag
        output(:, i) = input(:, flag(i));
    end
end





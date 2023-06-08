function [dec] = bi2dec(bi)
[N_num,N_bi] = size(bi);
dec = NaN(N_num,1);
for inum = 1 : N_num
    this_dec = 0;
    for ibi = 1 : N_bi
        if bi(inum,ibi)
            this_dec = this_dec + 2^(ibi-1);
        end
    end
    dec(inum) = this_dec;
end
end
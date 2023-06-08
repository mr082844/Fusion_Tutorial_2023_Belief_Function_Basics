function [bi] = dec2bi(dec)
[N_num] = length(dec);
max_dec = max(dec);
N_bi = floor(log(max_dec)/log(2)) + 1;
bi = false(N_num,N_bi);
for inum = 1 : N_num
    this_bi = false(1,N_bi);
    this_dec = dec(inum);
    for ibi = N_bi : -1 : 1
        if 2^(ibi-1) <= this_dec
            this_bi(ibi) = true;
            this_dec = this_dec - 2^(ibi-1);
        end
    end
    bi(inum,:) = this_bi;
end
end
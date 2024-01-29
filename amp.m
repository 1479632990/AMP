function  xmap = amp(M, N, H, y, pntset, nspw)
H1 = H';
H2 = abs(H) .^ 2;
H3 = abs(H1) .^ 2;
L = length(pntset);
u1 = zeros(M * N, 1);
v1 = ones(M * N, 1);
u3 = zeros(M * N, 1);
const = 1e-10;


for num = 1 : 20

    v2 = H2 * v1;
    u2 = H * u1 - v2 .* u3;  
    
    v3 = 1 ./ (v2 + nspw);
    u3 = (y - u2) .* v3;
    
    v4 = 1 ./ (H3 * v3);
    u4 = u1 + v4 .* (H1 * u3);
    
    logpp = - abs(u4 - pntset.') .^ 2 ./ v4;
    logpp = logpp - max(logpp,[],2);
    pp = exp(logpp);
    for j = 1 : M * N
        sumpp(j) = sum(pp(j,:));
        for a = 1 : L
            p(j, a) = pp(j, a) / sumpp(j);
        end
    end

    u1 = p * pntset;
    for j = 1 : M * N
        v1(j) = p(j, :) * abs(pntset - u1(j)) .^ 2;
    end


    if sum(v1) / (M * N) < const
       break;
    end   
end

xmap =  u1;
end

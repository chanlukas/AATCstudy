function [v] = sub(l,c,l1,c1)
    vec = 1:(l*c);
    vec=reshape(vec,c,l)';
    v=vec(l1,c1);
    v=v(:);
    subplot(l,c,v)
    
    %%%
    v = subplot(l,c,v);
    
    
end
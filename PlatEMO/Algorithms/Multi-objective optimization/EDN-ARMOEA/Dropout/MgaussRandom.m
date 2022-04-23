function [R1,R2] = MgaussRandom()

    r=0;
    while (r==0||r>1)
        u=rand*2-1; 
        v=rand*2-1;
        r=u^2+v^2;
    end

    c=sqrt(-2*log(r)/r);
    R2=v*c;
    R1=u*c;
end
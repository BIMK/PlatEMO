function array = iniA(Nx, Ny, flag, c)

    mu=0;
    if flag==0    
        std=sqrt(1/Nx);
        n=Nx*Ny;
        round=n/2;% n has to be an integer multiple of 2
        if round~=fix(round)
            printf('Warning: not a integer');
        end
        R=[];
        for i=1:round
            [r1,r2]=MgaussRandom;
            R=[R, r1, r2];
        end
        R=reshape(R, Nx, Ny);
        array=mu+R*std;

    else
        std=sqrt(1/Ny);
        n=Ny;   
        if c==0
            round=ceil(n/2);
            R=[];
            for i=1:round
                [r1,r2]=MgaussRandom;
                R=[R, r1, r2];
            end
        else
            R=ones(1,Ny)*c;
        end
        array=mu+R*std;
    end
end
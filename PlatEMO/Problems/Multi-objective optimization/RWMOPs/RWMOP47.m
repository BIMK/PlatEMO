classdef RWMOP47 < PROBLEM
% <multi> <real> <constrained>
% Optimal droop setting for minimizing active and reactive power loss

%------------------------------- Reference --------------------------------
% A. Kumar, G. Wu, M. Ali, Q. Luo, R. Mallipeddi, P. Suganthan, and S. Das,
% A benchmark-suite of real-world constrained multi-objective optimization
% problems and some baseline results, Swarm and Evolutionary Computation,
% 2021, 67: 100961.
%------------------------------- Copyright --------------------------------
% Copyright (c) 2023 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    methods
        %% Initialization
        function Setting(obj)
            obj.M        = 2;
            obj.D        = 18;
            obj.lower    = [-ones(1,10), zeros(1,2),   zeros(1,6)];
            obj.upper    = [ ones(1,10),2*ones(1,2),500*ones(1,6)];
            obj.encoding = ones(1,obj.D);
        end
        %% Evaluate multiple solutions
        function Population = Evaluation(obj,varargin)
            x = varargin{1};
            P = [  1.6138989835002382e-01   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000000e+00   0.0000000000000000e+00;
                   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000000e+00   0.0000000000000000e+00;
                   2.1451873678268588e-01   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000000e+00   0.0000000000000000e+00;
                   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000000e+00   0.0000000000000000e+00;
                   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000000e+00   0.0000000000000000e+00;
                   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000000e+00   0.0000000000000000e+00];
            Q = [  1.0680528035555388e-01   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000000e+00   0.0000000000000000e+00;
                   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000000e+00   0.0000000000000000e+00;
                   1.5161777012574437e-01   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000000e+00   0.0000000000000000e+00;
                   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000000e+00   0.0000000000000000e+00;
                   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000000e+00   0.0000000000000000e+00;
                   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   0.0000000000000000e+00   1.0000000000000000e+00   0.0000000000000000e+00];
            L = [  1.0000000000000000e+00   2.0000000000000000e+00   2.6660053320106641e-01   7.4329468658937317e-02   0.0000000000000000e+00   1.0000000000000000e+00;
                   1.0000000000000000e+00   4.0000000000000000e+00   1.8600037200074401e-01   8.1809163618327255e-02   0.0000000000000000e+00   1.0000000000000000e+00;
                   2.0000000000000000e+00   5.0000000000000000e+00   1.2400024800049600e-01   5.8435116870233741e-02   0.0000000000000000e+00   1.0000000000000000e+00;
                   2.0000000000000000e+00   3.0000000000000000e+00   9.3000186000372007e-02   4.3078368156736313e-01   0.0000000000000000e+00   1.0000000000000000e+00;
                   3.0000000000000000e+00   6.0000000000000000e+00   3.1000062000124000e-02   1.1687023374046750e-02   0.0000000000000000e+00   1.0000000000000000e+00];
            % Voltage initilization
            V    = zeros(6,1);
            V(1) = 1;
            Pc   = zeros(6,1);
            Qc   = zeros(6,1);
            for i = 1 : size(x,1)
                V(2:6)      = x(i,1:5)+1j*x(i,6:10);
                w           = x(i,11);
                V(1)        = x(i,12)+1e-5;
                Pc([4,5,6]) = x(i,13:15);
                Qc([4,5,6]) = x(i,16:18);
                % Current calculation
                Y     = ybus(L,w);
                I     = Y*V;
                Ir    = real(I);
                Im    = imag(I);
                Vr    = real(V);
                Vm    = imag(V);
                Psp   = Pc.*(1-w)-P(:,1).*(abs(V)./P(:,5)).^P(:,6);
                Qsp   = Qc.*(1-sqrt(Vr.^2+Vm.^2))-Q(:,1).*(abs(V)./Q(:,5)).^Q(:,6);
                spI   = conj((Psp+1j*Qsp)./V);
                spIr  = real(spI);
                spIm  = imag(spI);
                delIr = Ir-spIr;
                delIm = Im-spIm;
                % Objective calculation
                f(i,1) = sum(Psp) ;
                f(i,2) = sum(Qsp) ;
                h(i,:) = [delIr(1:6)',delIm(1:6)'];
            end
            Population = SOLUTION(varargin{1},f,abs(h)-1e-4,varargin{2:end});
            obj.FE     = obj.FE + length(Population);
        end
        %% Generate a point for hypervolume calculation
        function R = GetOptimum(obj,~)
            R = [9.3764479e-03   2.4737407e-03];
        end
    end
end

function Y = ybus(linedata,f)
    linedata(:,4) = linedata(:,4).*f;
    fb = linedata(:,1);             % From bus number...
    tb = linedata(:,2);             % To bus number...
    r  = linedata(:,3);             % Resistance, R...
    x  = linedata(:,4);             % Reactance, X...
    b  = linedata(:,5);             % Ground Admittance, B/2...
    a  = linedata(:,6);             % Tap setting value..
    z  = r + 1j*x;                  % z matrix...
    y  = 1./z;                      % To get inverse of each element...
    b  = 1j*b;                      % Make B imaginary...

    nb = max(max(fb),max(tb));      % No. of buses...
    nl = length(fb);                % No. of branches...
    Y = zeros(nb,nb);               % Initialise YBus...

    % Formation of the Off Diagonal Elements...
    for k = 1 : nl
        Y(fb(k),tb(k)) = Y(fb(k),tb(k)) - y(k)/a(k);
        Y(tb(k),fb(k)) = Y(fb(k),tb(k));
    end

    % Formation of Diagonal Elements....
    for m = 1 : nb
        for n = 1 : nl
            if fb(n) == m
                Y(m,m) = Y(m,m) + y(n)/(a(n)^2) + b(n);
            elseif tb(n) == m
                Y(m,m) = Y(m,m) + y(n) + b(n);
            end
        end
    end
end
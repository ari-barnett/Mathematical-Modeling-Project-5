%Project 5 Unsteady Aerodynamics 
%Barnett, Ari 


%Define Parameters for Amplitude, Omega, Time Span, and N
N = 100;
tspan = 0:0.01:10;
A = 0.1;
omega = 2;


%==========Define LHS and RHS matrices=====================

%LHS
[xn_grid, x_prime_k_grid] = meshgrid(xn(1:(N-1),N), x_prime_k(1:(N-1),N));
delta_xn_grid = repmat(delta_xn(1:(N-1),N), N-1, 1);
LHS = delta_xn_grid ./ (x_prime_k_grid - xn_grid);

%Exact Delta Xn values
dx = delta_xn_grid(1,:)';

%==========RHS========================

%Set inital values for vectors/matrices
GammaVector(1) = 0;
vortex_accumulation(1) = 0;

SumGamma(1) = 0;
SumGammaDeltaXn(1) = 0;

x_prime_k_values = x_prime_k_grid(:,1);

for m = 1:length(tspan)
    for M = 1:m
        xm_params = [0.01, m , M, 1, 0, 2.5];
        xm_collect(M) = xm(xm_params);
    end 


    vortex_accumulation = vortex_accumulation(:) + (GammaVector(m)./(x_prime_k_values(:) - xm_collect(m)));
    RHS_value= RHS(A,omega,vortex_accumulation,tspan(m));
    
    %Solve for little gamma_n
    gamma_n = LHS\RHS_value; 

    %Sovle and Update Captial Gamma Values
    SumGammaDeltaXn = sum(gamma_n.*dx);
    SumGamma = SumGamma + GammaVector(m); 

    GammaVector(m+1) = - SumGammaDeltaXn - SumGamma;

    %Collect Lift and Moment Coefficients
    LC_Vector(m) = LiftCoeffiecient(SumGammaDeltaXn);
    LM_Vector(m) = MomentCoeffiecient(xn_grid(1,:),gamma_n, delta_xn_grid(1,:));

end

figure()
plot(tspan,-LC_Vector,'r')
title("Time Accurate Lift")
xlabel("Time")
ylabel("CL")

hold on
plot(tspan,-2*pi*(A*omega*cos(omega*tspan)),'b')
grid on
legend("Unsteady", "Quasi-Steady")

figure()
plot(tspan,-LM_Vector,'r')
title("Time Accurate Moment")
xlabel("Time")
ylabel("CM")

hold on
plot(tspan,-((-0.5+0.5)*pi*A*omega*cos(omega*tspan)),'b')
grid on
legend("Unsteady", "Quasi-Steady")


%=========Collection of Functions to Use============ 
function xn_value = xn(n,N)
xn_value = -cos( ...
    (n-1)*pi / (N - 1));
end

function delta_xn_value = delta_xn(n,N)
delta_xn_value = -cos( ...
    ((n-1)*pi) / (N - 1)) + cos(n*pi/ (N-1));
end

function x_prime_k_value = x_prime_k(K,N)
x_prime_k_value = -(1/2)*(cos( ...
    (K-1)*pi / (N - 1)) + cos((K*pi) / (N-1)));
end

function rhs_value = RHS(A,omega,Gamma, t)
    rhs_value = (2*pi*A*omega*cos(omega*t)) - Gamma;
end

function xm_values = xm(parameters)
    time_delta = parameters(1);
    m = parameters(2);
    M = parameters(3);
    U_infinity = parameters(4);
    alpha = parameters(5);
    xs = parameters(6);

    xm_values = (time_delta*(M-m+1)*U_infinity*cos(alpha)) + xs;
end

function LC = LiftCoeffiecient(SumGammaDeltaXN)
    LC = SumGammaDeltaXN;
end

function LM = MomentCoeffiecient(xn,gamma,deltaXn)
    LM = (1/2)*sum((-0.5 - xn)*(gamma.*(deltaXn')));
end



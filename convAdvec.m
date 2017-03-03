function [p] = convAdvec(h,sigma,mu,q)

%   Input:
%           h: scalar. cellsize
%       sigma: vector. spatially varying coefficient, evaluated at cell centres.
%          mu: scalar. constant coefficient
%           q: vector. Source function, evaluated at nodes.
%  
%   Output: 
%           p: Approximate value of p on the nodes.
%           
%   You must write your discretization as a matrix equation A*p = q.
%   Once you have formed the matrix A, you can solve for p using the 
%   Matlab backslash operator. Make sure that you form A as a sparse
%   matrix. If A is dense the backslash operator will be very inefficient.

%% Creating the mesh with nodes and cell centers
n = 1/h;

xN = [0:h:1];

xC = [(h/2):h:(1-(h/2))];

xi = xN(2:n);

%% For the diffusion segment
Dnc=1/h*spdiags(ones((n),1)*[-1 1],[0,1],(n)-1,(n));

Dcn = -(Dnc)';

dif = Dnc*sigma*Dcn;

%% For the convection segment
Avg = (1/2)*spdiags(ones((n),1)*[1 1],[0,1],(n)-1,(n));

conv = mu.*Avg*Dcn;
%% Finding A
A = dif + conv;

% Finding P
p = A\q'


end

% Writing and running the test function from the other uploaded file, 
% it can be seen that the values of p outputted by the function code are very 
% close to the expected values given by the function p = convAdvec(h,sigma,mu,q)
% However, if h, the distance between nodes is increased to be greater than
% 10^-4, the function no longer runs as h has reached a value where the
% pressure at a point between nodes varies unexpectedly.
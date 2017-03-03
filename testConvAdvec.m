% Note: h needs to remain less than or equal to 10^-4 or else the function does
% not work

for i = 2:4
h=10.^(-i);
n=1/h;
% For x on the end nodes
xN = [0:h:1];

% For x on the cell centers
xC = [(h/2):h:(1-(h/2))];

% For x on the nodes except for the end nodes
xi = xN(2:n);

mu = 0.1;
%mu = 10;

% Creating sigma points
sig = (1+xC.^2)';
Sig = spdiags(sig, [0], n, n);
% 
dif = 4*pi*xi.*cos(2*pi*xi)-4*pi^2*(xi.^2 + 1).*sin(2*pi*xi);
adv = mu*2*pi*cos(2*pi*xi);

% q, the source function is the sum of the diffusion and the advection
% values
q =  dif + adv;
pExpected = convAdvec(h,Sig,mu,q);
pActual = sin(2*pi*xC);
i+1;

figure
hold on
plot(pExpected)
plot(pActual)
hold off
pause(10);
end
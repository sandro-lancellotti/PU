% Example 1
d = 2; n = [65].^d; % Define the space dimension and the number of data
p = haltonset(d); x = net(p,n); % Define the interpolation data (Halton points)
phi = @(epsilon,r) (1+epsilon*r).*exp(-epsilon*r); epsilon =1;% Define the kernel
franke = @(x1,x2) 0.75 * exp(-(9*x1-2).^2/4 - (9*x2-2).^2/4)+ 0.75 * exp(-(9*x1+1).^2/49 - ...
(9*x2+1)./10)+0.5 * exp(-(9*x1-7).^2/4 - (9*x2-3).^2/4)-0.2 * exp(-(9*x1-4).^2 - (9*x2-7).^2);
f = franke(x(:,1),x(:,2)); % Define the test function and the function values
m_d = floor(n^(1/d)/2); s_d = 60; % Define the number of patches and evaluation data in one direction
w = @(supp,r) (max(1-(supp*r),0).^4).*(4*(supp*r)+1); % Define the PU weights (Wendland C^2)
bar_x = MakeSDGrid(d,s_d); % Create s_d^d equally spaced test data
Pf = PU(d,x,bar_x,m_d,phi,w,f,epsilon); % Compute the PU interpolant
% Compute the maximum error
MAE = max(abs(franke(bar_x(:,1),bar_x(:,2))-Pf));
% Plot the result
figure 
view(3)
grid on; hold on; axis square 
xx = linspace(0,1,s_d);
[X,Y] = meshgrid(xx,xx);
surf(X,Y,reshape(Pf,s_d,s_d),'EdgeColor','none')
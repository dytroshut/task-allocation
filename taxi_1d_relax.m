%% taxi allocation one-dimensional case
clear all
clc

%%Random distribution generator
n1 = 60;  % (not exceed 1000)
n2 = 60;  % (not exceed 1000)

x1 = (0.01:n1-0.01)'/n1; % x \in [-1,0]
x2 = (0.01:n2-0.01)'/n2;     % y \in [ 0,1]

nx1 = length(x1);
nx2 = length(x2);

% Useful functions 
normalize = @(a)a/sum(a,'all');
Gaussian = @(x,t0,sigma)exp( -(x-t0).^2/(2*sigma^2) );
sigma0 = 0.01; sigma1 = .02; sigma2 = .03; sigma3 = .04; sigma4 = .05;

w = 3.*Gaussian(x1, .5, .07);
w1 = 2.*Gaussian(x1, .3, .03);
w2 = 5.*Gaussian(x1, .7, .08);
w = w*w' + w1*w1' + w2*w2'+w1*w2';
% w = w*w' + w2*w2'+w1*w2';
w = normalize(w);

ODconstraint = w;

p = sum(ODconstraint',2)';
q = sum(ODconstraint',1);
% Plot the two distributions p and q

figure(1)
heatmap(ODconstraint)

%%
% Discretization of the time mariginal by nt (not exceed 500)
nt = 60;

% Define the time interval [0,T], notice that T is not neccesary to be 1.
y = (0.01:nt-0.01)'/nt;
ny = length(y);

% wt = 5.*Gaussian(y, .3, .4) + 3.*Gaussian(y, .8, .1);
wt = 3.*ones(nt,1);
wt = normalize(wt);

% Compute the 3-dimensional Cost matrix
% [X1,X2,Y] = meshgrid(x1,x2,y);
[X1,X2,Y] = ndgrid(x1,x2,y);
% C = (Y-X1).^2 + (X2-Y).^2 + (X1-X2).^2;
C = abs(Y-X1) + abs(X2-Y) + abs(X1-X2);

% The order of the x-y axis need to be changed by the permute function,
Cost = permute(C,[2,1,3]);

%%
cvx_begin
cvx_solver mosek
cvx_precision best
    variable Pi(nx1,nx2,ny)
    variable rk(ny)
    minimize (sum(sum(sum(Cost.*Pi))))
    subject to  
    % 1.Marinigal on t-coordinate
    squeeze(sum(sum(Pi))) <= 2.*wt;
    % 2.Element-wise lower bound for \Pi
    Pi >= 0;
    % 3. DA constraint
    sum(Pi,3) == ODconstraint;
cvx_end


% Check 3 marginal constraints:
fprintf('Constraints deviation (should be 0):'); 

p = squeeze(sum(sum(permute(Pi,[2,3,1]))));
sum(p)
q = squeeze(sum(sum(permute(Pi,[3,1,2]))));
sum(q)

rk = squeeze(sum(sum(Pi)));


% Plots of rk, evolutionary probability distributions at moment t
% 3D plot with 3 marginals
% Generate new coordinates for 3D plot
xx1 = x1;
xx2 = zeros(size(xx1));
xx = [xx1,xx2,p];
yy1 = x2;
yy2 = ones(size(yy1));
yy = [yy1,yy2,q];
tt1 = y;
tollheight = w;

%%
list = find (Pi>0.000001);
[I,J,KK] = ind2sub(size(Pi),list);
tt2 = ones(size(tt1));

[X,Y] = meshgrid(x1,x2);
%
figure(2)
scalar = 4;
p1 = plot3(xx1,xx2,scalar.*p,'LineWidth',2,'Color','g');
hold on
p2 = plot3(yy2,yy1,scalar.*q,'LineWidth',2,'Color','r');
S = surf(X,Y,20.*ODconstraint,'FaceAlpha',0.95);
S.EdgeColor = 'none';

tplot = plot3(5*rk,tt2,tt1,'LineWidth',2,'Color','b');

 for i=1:length(I)
    pt_assignement = plot3(x1(I(i)),x2(J(i)),y(KK(i)),'o','MarkerSize',...
    2000*Pi(I(i),J(i),KK(i)),'MarkerFaceColor','none', 'LineWidth',4);
 end

hold off
% xlabel('x^1','FontSize',20)
% ylabel('x^2','FontSize',20)
% zlabel('y','FontSize',20)
axis tight
view([43 26])

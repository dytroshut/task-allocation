%% taxi allocation one-dimensional case
clear all
clc

%%Random distribution generator
n = 10;  % (not exceed 1000)

x = (0.01:n-0.01)'/n; 

nx = length(x);

% Useful functions 
normalize = @(a)a/sum(a,'all');
Gaussian = @(x,t0,sigma)exp( -(x-t0).^2/(2*sigma^2) );
sigma0 = 0.01; sigma1 = .02; sigma2 = .03; sigma3 = .04; sigma4 = .05;

w0 = 3.*Gaussian(x, .2, .08);
w0 = normalize(w0);
w1 = 2.*Gaussian(x, .3, .1);
w1 = normalize(w1);
w2 = 5.*Gaussian(x, .5, .2);
w2 = normalize(w2);
w3 = 2.*Gaussian(x, .7, .15);
w3 = normalize(w3);
w4 = 2.*Gaussian(x, .8, .09);
w4 = normalize(w4);

od0 = tensorprod(w2,w2,2);
od0 = tensorprod(od0,w2,3);
od0 = tensorprod(od0,w2,4);

od1 = tensorprod(w2,w2,2);
od1 = tensorprod(od1,w3,3);
od1 = tensorprod(od1,w3,4);

od2 = tensorprod(w3,w3,2);
od2 = tensorprod(od2,w3,3);
od2 = tensorprod(od2,w4,4);

od3 = tensorprod(w4,w4,2);
od3 = tensorprod(od3,w4,3);
od3 = tensorprod(od3,w4,4);

wod = od0 + od1 + od2 + od3;

% awod = od0;
wod = normalize(wod);
wo = sum(wod,[3,4]);
wd = squeeze(sum(wod,[1,2]));

wt = tensorprod(w1,w2,2) + tensorprod(w2,w3,2);
% w2 = Gaussian(x, .5, .5);
% wt = tensorprod(w2,w2,2);
wt = normalize(wt);

figure(1)
x = 1:n;
y = 1:n;
hold on





%%
% Compute the 3-dimensional Cost matrix
[X1,X2,X3,X4,Y1,Y2] = ndgrid(x,x,x,x,x,x);

% 
% cost = -2.*Y1.*X1 - 2.*X2.*Y2 - 2.*Y1.*X3 - 2.*X4.*Y2;
cost = (Y1-X1).^2 + (X2-Y2).^2 + (X1-X3).^2 + (X2-X4).^2 + (Y1-X3).^2 + (X4-Y2).^2;

% The order of the x-y axis need to be changed by the permute function,
cost = permute(cost,[3,4,1,2,5,6]);
% cost = permute(cost,[6,5,4,3,2,1]);

cvx_begin
% cvx_solver SDPT3
cvx_solver mosek
%cvx_precision best
    variable Pi(n,n,n,n,n,n)
    minimize (sum(sum(sum(sum(sum(sum(cost.*Pi)))))))
    % minimize sum(sum(sum(sum(sum(sum(-entr(Pi)))))))
    subject to  
    % 1.Marinigals
    squeeze(sum(sum(sum(sum(permute(Pi,[1,2,3,4,5,6])))))) == wt;
    squeeze(sum(sum(permute(Pi,[5,6,1,2,3,4])))) == wod;
    % 2.Non-negative element
    Pi >= 0;
cvx_end


% Check 3 marginal constraints:
%fprintf('Constraints deviation (should be 0):'); 

p = squeeze(sum(sum(sum(sum(Pi)))))-wt;
sum(p,'all')
q = squeeze(sum(sum(permute(Pi,[5,6,1,2,3,4]))))-wod;
sum(q,'all')


figure(2)

%% taxi allocation one-dimensional case
clear all
clc
clf 
%%Random distribution generator
n = 10; 

x = (0.01:n-0.01)'/n; 

nx = length(x);

% Useful functions 
normalize = @(a)a/sum(a,'all');

wod = zeros(n,n,n,n);
wod(7,1,7,8) = 2;
wod(9,4,7,10) = 3;
wod(8,2,7,8) = 2;
% wod(10,3,9,8) = 3;
% wod(9,1,7,8) = 2;
% wod(10,2,8,9) = 2;
wod = normalize(wod);


wt = zeros(n,n);
wt(1,3) = 1;
wt(2,5) = 2;
wt(1,7) = 4;
wt = normalize(wt);


figure(1)
beforeplan(wod, wt)


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



%%
figure(2)
afterplan1(wod, wt, Pi)
%% functions

function beforeplan(wod, wt)
wo = sum(wod,[3,4]);
wd = squeeze(sum(wod,[1,2]));

hold on;

linearIndices = find(wod);
[I1, J1, K1, L1] = ind2sub(size(wod), linearIndices);
numLines = length(I1);
for i = 1:numLines
    plot([J1(i), L1(i)],[I1(i),K1(i)],'LineWidth',10.*wod(I1(i),J1(i),K1(i),L1(i)),'Color','k');
end

   
[row, col, values] = find(wo);
markerSize = abs(values) * 1500;
scatter(col, row, markerSize, values, 'o', ...
        'MarkerEdgeColor', [0.6350 0.0780 0.1840], ...
        'MarkerFaceColor', [0.6350 0.0780 0.1840]);

[row, col, values] = find(wd);
markerSize = abs(values) * 1500;
scatter(col, row, markerSize, values, 'o', ...
        'MarkerEdgeColor', [0 0.4470 0.7410], ...
        'MarkerFaceColor', [0 0.4470 0.7410]);

[row, col, values] = find(wt);
markerSize = abs(values) * 1500;
scatter(col, row, markerSize, 0.8.*values, "^", ...
        'MarkerEdgeColor', "#77AC30", ...
        'MarkerFaceColor', "#77AC30");

hold off
axis off
end



function afterplan1(wod, wt, Pi)
wo = sum(wod,[3,4]);
wd = squeeze(sum(wod,[1,2]));

hold on;

linearIndices1 = find(Pi);
[I, J, K, L, H, W] = ind2sub(size(Pi), linearIndices1);


numLines = length(I);
hold on;

for i = 1:numLines
    x = [J(i),L(i),W(i),J(i)];
    y = [I(i),K(i),H(i),I(i)];
    plot( x,y,'LineWidth',10.*Pi(I(i),J(i),K(i),L(i),H(i),W(i))  ) ;
end

linearIndices = find(wod);
[I1, J1, K1, L1] = ind2sub(size(wod), linearIndices);
numLines = length(I1);
for i = 1:numLines
    plot([J1(i), L1(i)],[I1(i),K1(i)],'LineWidth',10.*wod(I1(i),J1(i),K1(i),L1(i)),'Color','k');
end



[row, col, values] = find(wo);
markerSize = abs(values) * 1500;
scatter(col, row, markerSize, values, 'o', ...
        'MarkerEdgeColor', [0.6350 0.0780 0.1840], ...
        'MarkerFaceColor', [0.6350 0.0780 0.1840]);

   
[row, col, values] = find(wd);    
markerSize = abs(values) * 1500;
scatter(col, row, markerSize, values, 'o', ...
        'MarkerEdgeColor', [0 0.4470 0.7410], ...
        'MarkerFaceColor', [0 0.4470 0.7410]);

   
[row, col, values] = find(wt);
markerSize = abs(values) * 1500;
scatter(col, row, markerSize, 0.8.*values, "^", ...
        'MarkerEdgeColor', "#77AC30", ...
        'MarkerFaceColor', "#77AC30");


hold off
axis off
end

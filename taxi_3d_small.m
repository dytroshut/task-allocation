%% taxi allocation one-dimensional case
clear all
clc
clf 
%%Random distribution generator
n = 7; 

x = (0.01:n-0.01)'/n; 

nx = length(x);

% Useful functions 
normalize = @(a)a/sum(a,'all');

wod = zeros(n,n,n,n,n,n);
wod(4,5,4,1,2,5) = 2;
wod(5,6,5,2,3,6) = 1;
wod(6,7,6,3,4,7) = 2;
wod = normalize(wod);

%
wt = zeros(n,n,n);
wt(5,2,1) = 1;
wt(2,1,1) = 2;
wt(1,7,1) = 2;
wt = normalize(wt);


%%
% Compute the 3-dimensional Cost matrix
[X1,X2,X3,X4,X5,X6,Y1,Y2,Y3] = ndgrid(x,x,x,x,x,x,x,x,x);

cost = (Y1-X1).^2 + (X2-Y2).^2 + (X3-Y3).^2 + (X4-Y1).^2 + (X5-Y2).^2 + (X6-Y3).^2;

% The order of the x-y axis need to be changed by the permute function,
cost = permute(cost,[4,5,6,1,2,3,7,8,9]);
% cost = permute(cost,[6,5,4,3,2,1]);

cvx_begin
cvx_solver mosek
    variable Pi(n,n,n,n,n,n,n,n,n)
    minimize (sum(sum(sum(sum(sum(sum(sum(sum(sum(cost.*Pi))))))))))
    subject to  
    % 1.Marinigals
    squeeze(sum(sum(sum(sum(sum(sum(permute(Pi,[1,2,3,4,5,6,7,8,9])))))))) == wt;
    squeeze(sum(sum(sum(permute(Pi,[7,8,9,1,2,3,4,5,6]))))) == wod;
    % 2.Non-negative element
    Pi >= 0;
cvx_end


%%

figure(1)
beforeplan(wod, wt)

figure(2)
afterplan1(wod, wt, Pi)

%% functions

function beforeplan(wod, wt)
wo = sum(wod,[4,5,6]);
wd = squeeze(sum(wod,[1,2,3]));

hold on;
linearIndices = find(wod);
[I1, J1, K1, L1, M1, N1] = ind2sub(size(wod), linearIndices);

numLines = length(I1);

for i = 1:numLines
    plot3([I1(i),L1(i)],[J1(i),M1(i)],[K1(i),N1(i)],'LineWidth',10.*wod(I1(i),J1(i),K1(i),L1(i),M1(i),N1(i)),'Color','k');
end

linearIndices = find(wo);
values = wo(linearIndices);
[row, col, mod] = ind2sub(size(wo), linearIndices);
markerSize = abs(values) * 1500;
scatter3(row, col, mod, markerSize, values, 'o', ...
        'MarkerEdgeColor', [0.6350 0.0780 0.1840], ...
        'MarkerFaceColor', [0.6350 0.0780 0.1840]);


linearIndices = find(wd);
values = wd(linearIndices);
[row, col, mod] = ind2sub(size(wd), linearIndices);
markerSize = abs(values) * 1500;
scatter3(row, col, mod, markerSize, values, 'o', ...
        'MarkerEdgeColor', [0 0.4470 0.7410], ...
        'MarkerFaceColor', [0 0.4470 0.7410]);

linearIndices = find(wt);
values = wt(linearIndices);
[row, col, mod] = ind2sub(size(wt), linearIndices);
markerSize = abs(values) * 1500;
scatter3(row, col, mod, markerSize, 0.8.*values, "^", ...
        'MarkerEdgeColor', "#77AC30", ...
        'MarkerFaceColor', "#77AC30");

hold off
%axis off
view([69 20])
ax = gca;
ax.FontSize = 15; 
end



function afterplan1(wod, wt, Pi)
wo = sum(wod,[4,5,6]);
wd = squeeze(sum(wod,[1,2,3]));

hold on;

linearIndices1 = find(Pi);
[I, J, M, K, L, N, H, W, O] = ind2sub(size(Pi), linearIndices1);


numLines = length(I);

for i = 1:numLines
y = [J(i),L(i),W(i),J(i)];
x = [I(i),K(i),H(i),I(i)];
z = [M(i),N(i),O(i),M(i)];
plot3(x,y,z,'LineWidth',10.*Pi(I(i), J(i), M(i), K(i), L(i), N(i), H(i), W(i), O(i))) ;
end



linearIndices = find(wod);
[I1, J1, K1, L1, M1, N1] = ind2sub(size(wod), linearIndices);

numLines = length(I1);

for i = 1:numLines
    plot3([I1(i),L1(i)],[J1(i),M1(i)],[K1(i),N1(i)],'LineWidth',10.*wod(I1(i),J1(i),K1(i),L1(i),M1(i),N1(i)),'Color','k');
end

%
linearIndices = find(wo);
values = wo(linearIndices);
[row, col, mod] = ind2sub(size(wo), linearIndices);
markerSize = abs(values) * 1500;
scatter3(row, col, mod, markerSize, values, 'o', ...
        'MarkerEdgeColor', [0.6350 0.0780 0.1840], ...
        'MarkerFaceColor', [0.6350 0.0780 0.1840]);


linearIndices = find(wd);
values = wd(linearIndices);
[row, col, mod] = ind2sub(size(wd), linearIndices);
markerSize = abs(values) * 1500;
scatter3(row, col, mod, markerSize, values, 'o', ...
        'MarkerEdgeColor', [0 0.4470 0.7410], ...
        'MarkerFaceColor', [0 0.4470 0.7410]);

linearIndices = find(wt);
values = wt(linearIndices);
[row, col, mod] = ind2sub(size(wt), linearIndices);
markerSize = abs(values) * 1500;
scatter3(row, col, mod, markerSize, 0.8.*values, "^", ...
        'MarkerEdgeColor', "#77AC30", ...
        'MarkerFaceColor', "#77AC30");

hold off
% axis off
view([69 20])
ax = gca;
ax.FontSize = 15; 
end

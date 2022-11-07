close all
Sections = {'Rectangular'; 'Circular'; 'WideFlange'};

np = 10000;
dku = 0.15;
dk = dku/np;

E  = 2.10E11;
ey = 0.002;
Sy = E*ey;

kappa = 0.0:dk:dku;
M = zeros(np, length(Sections));
for j = 1:length(Sections)

    vals = load(Sections{j});
    yi = vals(:,1);
    zi = vals(:,2);
    Ai = vals(:,3);

    %Section's geometric center
    xcm = sum(yi.*Ai)/sum(Ai);
    zcm = sum(zi.*Ai)/sum(Ai);

    %Computes the moment for each curvature
    i = 0;
    for k = kappa
        i = i + 1;
        ei = k*(zi - zcm);
        Si = E*ei.*(abs(ei) <= ey) + Sy*(ei > ey) - Sy*(ei < -ey);
        Fi = Si.*Ai;
        M(i,j) = sum(Fi.*(zi - zcm));
    end

    %Finds the yield moment and normalize
    ind = find(abs(kappa - 1.425E-2) < 10*eps);
    M(:,j) = M(:,j)/M(ind,j);
end

%Plots the Moment-Curvature diagram
figure(1)
plot(0.0:dk:dku ,M, 'LineWidth',2)
xlabel('$\kappa$','Interpreter','latex', 'FontSize', 22)
ylabel('$\frac{M_p}{M_y}$','Interpreter','latex', 'FontSize', 22)
legend(Sections,'Location','southeast','Interpreter','latex', 'FontSize', 16)
%ylim([0,2])
grid on;
box on;
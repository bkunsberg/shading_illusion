% This script creates two illusions based on the concept of 'critical
% contours'.  The critical contours are separatrices of the shaded image
% that perceptually segment the object into 'bumps' or 'dimples'.

% Ben Kunsberg
% 10/17/18
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Single Layer 'Lego' Illusion
% Here, there will be a critical contour loop (circle) around the main
% feature.  This critical contour loop will segment the inner region into
% either a 'bump' or a 'dimple'.  Both of these percepts can be seen in the
% one image.  The key is that the critical contour represents a boundary
% and does not define whether the interior is 'in' or 'out'.  Note that is
% not just the 'concave/convex' illusion as we are only flipping a part of
% the surface.

%Define 'components'
n = 1000;
[X, Y] = meshgrid(linspace(-1, 1, n), linspace(-1, 1, n));
factor = 25;
t = sqrt(X.^2 + Y.^2);
z2 = 1./(5*(1 + exp(-factor*(t-0.4))));
radius = 0.7;


%Define the surface (we make it not symmetrical so that the effect is
%'generic')
sinx = 0.5*sin(Y + .5);
base =  sqrt(max(0, 1 - radius*X.^2 - radius*Y.^2));
Z = -(- z2 - base + sinx);

%Render Image
L1 = [0.5, -0.5, 2];
L2 = [0.3, -0.3, 1];
render_image(L1, L2, Z, n, X, Y, 1)


%% Multi 'Lego' Illusion
% We repeat the above 'lego' illusion with two bump/dimples.  Both
% segmented regions have the bistability.

%Define 'components'
n = 1000;
[X, Y] = meshgrid(linspace(-1, 1, n), linspace(-1, 1, n));
factor = 25;
center1 = [0.5, 0];
t1 = sqrt((X - center1(1)).^2 + (Y - center1(2)).^2);
center2 = [-0.5, 0];
t2 = sqrt((X - center2(1)).^2 + (Y - center2(2)).^2);
z1 = 1./(5*(1 + exp(-factor*(t1-0.3))));
z2 = 1./(5*(1 + exp(-factor*(t2-0.3))));
radius = 0.7;


%Define the surface (we make it not symmetrical so that the effect is
%'generic')
sinx = 0.5*sin(Y + .5);
base =  sqrt(max(0, 1 - radius*X.^2 - radius*Y.^2));
Z = -(- z2 - z1 - base + sinx);

%Render Image
L1 = [0.5, -0.5, 2];
L2 = [0.2, -0.2, 1];
render_image(L1, L2, Z, n, X, Y, 3)

%% Wedding Cake Illusion

% Similar to the above 'lego' illusion but we show the bumps can be nested.
% There are 4 possible 3D perceptions.


%Define 'components'
n = 1000;
[X, Y] = meshgrid(linspace(-1, 1, n), linspace(-1, 1, n));
t = sqrt(X.^2 + Y.^2);
z2 = 1./(5*(1 + exp(-factor*(t-0.6))));
z3 = 1./(7*(1 + exp(-factor*(t-0.2))));
radius = 0.7;

%Define the surface (we make it not symmetrical so that the effect is
%'generic')
sinx = 0.5*sin(Y + .5);
base =  sqrt(max(0, 1 - radius*X.^2 - radius*Y.^2));
Z = -(- z2 - z3 - base + sinx);

%Render the image
L1 = [0, -0.3, 2];
L2 = [0.2, 0, 1];
render_image(L1, L2, Z, n, X, Y, 2)



%%
function render_image(L1, L2, Z, n, X, Y, type)
[Zx, Zy] = gradient(Z, 1/n, 1/n);
denom = sqrt(1 + Zx.^2 + Zy.^2);
N1 = -Zx./denom;
N2 = -Zy./denom;
N3 = 1./denom;
I1 = L1(1) * N1 + L1(2) * N2 + L1(3) * N3;
I2 = L2(1) * N1 + L2(2) * N2 + L2(3) * N3;

fullI = 0*max(I1, 0) + 0.4*max(I2, 0) + 0*(1-(X.^2 + Y.^2));
fullI = (fullI - min(fullI(:)))./(max(fullI(:)) - min(fullI(:)));

figure('units','normalized','outerposition',[0 0 0.9 0.5]); subtightplot(1, 3, 1); surf(Z, 'EdgeColor', 'none'); title('The Surface'); subtightplot(1, 3, 2); imshow(fullI); title('The associated image (two light sources)');



if type == 1
    %Define the critical contours, lego
    t= linspace (0, 2*pi, 1000);
    X = 0.4*n/2 * cos(t) + 500;
    Y = 0.4*n/2 * sin(t) + 500;
    subtightplot(1, 3, 3); imshow(fullI); hold on;  plot(X,Y, 'LineWidth', 6); title('Critical Contour')
elseif type ==2
    %Define the critical contours, wedding cake
    t= linspace (0, 2*pi, 1000);
    X = 0.6*n/2 * cos(t) + 500;
    Y = 0.6*n/2 * sin(t) + 500;
    X2 = 0.2*n/2 * cos(t) + 500;
    Y2 = 0.2*n/2 * sin(t) + 500;
    subtightplot(1, 3, 3); imshow(fullI); hold on;  plot(X,Y, 'LineWidth', 6); plot(X2,Y2, 'LineWidth', 6); title('Critical Contours')
elseif type == 3
    %Define the critical contours, multi-lego
    t= linspace (0, 2*pi, 1000);
    X = 0.3*n/2 * cos(t) + 750;
    Y = 0.3*n/2 * sin(t) + 500;
    X2 = 0.3*n/2 * cos(t) + 250;
    Y2 = 0.3*n/2 * sin(t) + 500;
    subtightplot(1, 3, 3); imshow(fullI); hold on;  plot(X,Y, 'LineWidth', 6); plot(X2,Y2, 'LineWidth', 6); title('Critical Contours')
end

end
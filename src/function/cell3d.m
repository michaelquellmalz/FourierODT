% [x,y,z] = meshgrid(linspace(-1,1,160));
function f = cell3d(x,y,z,radius)

x = x/radius;
y = y/radius;
z = z/radius;

% Cytoplasm (ellipsoid)
r1 = .9;
r2 = .75;
r3 = .75;
c  = 1;
f = (x.^2/r1.^2 + y.^2/r2.^2 + z.^2/r3.^2 <=1) * 0.85;

% Nucleus (ellipsoid)
r1 = .4; % radius
r2 = .35;
r3 = .3;
c1 = .3; % center
c2 = .2;
f((x-c1).^2/r1.^2 + (y-c2).^2/r2.^2 + z.^2/r3.^2 <=1) = .5;

% Nucleolus left (little ball)
r1 = .1; % radius
r2 = .1;
r3 = .1;
f_nucleoleus = 1;
c1 = .15; % center
c2 = .15;
c3 = 0;
f((x-c1).^2/r1.^2 + (y-c2).^2/r2.^2 + (z-c3).^2/r3.^2 <=1) = f_nucleoleus;

% Nucleolus right
c1 = .45; % center
c2 = .16;
c3 = 0;
f((x-c1).^2/r1.^2 + (y-c2).^2/r2.^2 + (z-c3).^2/r3.^2 <=1) = f_nucleoleus;

% Nucleolus top
c1 = .3; % center
c2 = .38;
c3 = 0;
f((x-c1).^2/r1.^2 + (y-c2).^2/r2.^2 + (z-c3).^2/r3.^2 <=1) = f_nucleoleus;

% Mitochondria (peanut, Cassini oval)
c1 = -.5;
a = .25;
b = .08;
f(((y-a).^2+(x-c1).^2+z.^2).*((y+a).^2+(x-c1).^2+z.^2) <= b.^2) = .35;

% Centriole (star cylinder)
% https://nyjp07.com/index_asteroid_E.html
c1 = 0;
c2 = -.4;
r1 = .12;
r3 = .2;
corners = 7;
a = .8;
b = .6;
r = sqrt((x-c1).^2 + (y-c2).^2) / r1;
phi = atan2(y-c2,x-c1);
f((r <= sqrt(-log(2*exp(-a^2)-exp(-b^2.*sin((phi-pi/2)*corners/2).^2)))) & (z.^2<=r3^2)) = .65;

% Ellipsoid
r1 = .4; % radius
r2 = .3;
r3 = .2;
c1 = -.2; % center
c2 = 0;
c3 = .4;
f((x-c1).^2/r1.^2 + (y-c2).^2/r2.^2 + (z-c3).^2/r3.^2 <=1) = .4;

% Rectangle
r1 = .2; % radius
r2 = .2;
r3 = .2;
c1 = -.3; % center
c2 = 0.05;
c3 = -.35;
f(((x-c1).^2<r1.^2) & ((y-c2).^2<r2.^2) & ((z-c3).^2<r3.^2)) = .6;

end % function


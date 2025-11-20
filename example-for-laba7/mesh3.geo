Mesh.MshFileVersion = 2.0;

L=2;
r=0.1;
p1 = L/30;
p2 = r/5;

Point(1) = {0, 0, 0, p1};
Point(2) = {L, 0, 0, p1};
Point(3) = {L, L, 0, p1};
Point(4) = {r, L, 0, p2};
Point(5) = {r, L/2, 0, p2};
Point(6) = {0, L/2, 0, p2};
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};
Line Loop(7) = {4, 5, 6, 1, 2, 3};
Plane Surface(8) = {7};
Physical Line(1) = {4 , 5};
Physical Line(2) = {3, 6, 1, 2};
Physical Surface(1) = {8};
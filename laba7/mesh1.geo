Mesh.MshFileVersion = 2.0;

L=2;
r=0.1;
p1 = L/10;
p2 = r/0.2;

Point(1) = {0, 0, 0, p1};
Point(2) = {L, 0, 0, p1};
Point(3) = {L-r, L, 0, p2};
Point(4) = {r, L, 0, p2};
Point(5) = {r, L/2, 0, p2};
Point(6) = {0, L/2, 0, p2};
Point(7) = {L-r, L/2, 0, p2};
Point(8) = {L, L/2, 0, p2};//+
Line(1) = {1, 6};
//+
Line(2) = {6, 5};
//+
Line(3) = {5, 4};
//+
Line(4) = {4, 3};
//+
Line(5) = {3, 7};
//+
Line(6) = {7, 8};
//+
Line(7) = {8, 2};
//+
Line(8) = {2, 1};
//+
Curve Loop(1) = {3, 4, 5, 6, 7, 8, 1, 2};
//+
Plane Surface(1) = {1};
//+
Physical Curve(1) = {3, 5};
//+
Physical Curve(2) = {4, 1, 8, 7};
//+
Physical Surface(1) = {1};

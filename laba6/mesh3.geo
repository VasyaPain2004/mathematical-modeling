Mesh.MshFileVersion = 2.0;

p = 0.01;

Point(1) = {0, 0, 0, p};
Point(2) = {1, 0, 0, p};
Point(3) = {0, 1, 0, p};
Point(4) = {2, 0, 0, p};
Point(5) = {0, 2, 0, p};

Circle(6) = {2, 1, 3};
Circle(7) = {4, 1, 5};

Line(1) = {2, 4};
Line(2) = {3, 5};
//+
Curve Loop(1) = {2, -7, -1, 6};
//+
Plane Surface(1) = {1};
//+
Physical Curve(1) = {2};
//+
Physical Curve(2) = {7};
//+
Physical Curve(3) = {1};
//+
Physical Curve(4) = {6};

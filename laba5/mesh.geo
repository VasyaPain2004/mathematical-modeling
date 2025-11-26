Mesh.MshFileVersion = 2.0;

LX = 10.0;
LY = 5.0;
p = LY/50.0;

Point(1) = { 0, 0, 0, p};
Point(2) = { 0, LY, 0, p};
Point(3) = { LX, LY, 0, p};
Point(4) = { LX, 0, 0, p};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Line Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Physical Surface(1) = {1};

Physical Line(1) = {1};
Physical Line(2) = {3};
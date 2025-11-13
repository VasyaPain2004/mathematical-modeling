Mesh.MshFileVersion = 2.0;

p = 0.01;

Point(1) = {0, 0, 0, p};
Point(2) = {1, 0, 0, p};
Point(3) = {0, 1, 0, p};
Point(4) = {2, 0, 0, p};
Point(5) = {0, 2, 0, p};


//+
Line(1) = {5, 3};
//+
Line(2) = {2, 4};
//+
Circle(3) = {3, 1, 2};
//+
Circle(4) = {5, 1, 4};
//+
Curve Loop(1) = {1, 3, 2, -4};
//+
Plane Surface(1) = {1};
//+
Physical Curve(1) = {1};
//+
Physical Curve(2) = {4};
//+
Physical Curve(3) = {2};
//+
Physical Curve(4) = {3};

Physical Surface(1) = {1};
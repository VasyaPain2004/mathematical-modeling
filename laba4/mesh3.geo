Mesh.MshFileVersion = 2.0;
SetFactory("OpenCASCADE");

LX = 10.0;
LY = 5.0;
p = LY/150.0;

Point(1) = {0, 0, 0, p};
Point(2) = {0, LY, 0, p};
Point(3) = {(LX / 2) - 2, LY, 0, p};
Point(4) = {(LX / 2) + 2, LY, 0, p};
Point(5) = {LX, LY, 0, p};
Point(6) = {LX, 0, 0, p};

Line(7) = {1, 2};
Line(8) = {2, 3};
Line(9) = {4, 5};
Line(10) = {5, 6};
Line(11) = {6, 1};

Point(12) = {5, 5, 0, p}; // Center
Point(13) = {7, 5, 0, p};  // Right point
Point(14) = {3, 5, 0, p};  // Left point
Point(15) = {5, 3, 0, p}; // Bottom point
Point(16) = {7, 5, 0, p}; 
//+
Circle(12) = {15, 12, 3};
//+
Circle(13) = {15, 12, 4};
//+
Curve Loop(1) = {7, 8, -12, 13, 9, 10, 11};
//+
Plane Surface(1) = {1};
//+
Physical Curve(1) = {7};
//+
Physical Curve(2) = {11};
//+
Physical Curve(3) = {10};
//+
Physical Curve(4) = {9, 13, 12, 8};
//+
Physical Surface(1) = {1};
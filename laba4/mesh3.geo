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
Point(15) = {5, 3.5, 0, p}; // Bottom point
Point(16) = {7, 5, 0, p}; 

Ellipse(17) = {13, 12, 16, 15}; // правая нижняя четверть
Ellipse(18) = {15, 12, 16, 14}; // левая нижняя четверть

Line Loop(1) = {7, 8, -18, -17, 9, 10, 11};

Plane Surface(1) = {1};

Physical Surface(1) = {1};

Physical Line(1) = {7}; // левая граница
Physical Line(2) = {8, 9}; // верхняя граница
Physical Line(3) = {10}; // правая граница
Physical Line(4) = {11}; // нижняя граница
Physical Line(5) = {17, 18}; // эллиптический вырез

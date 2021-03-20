hole_radius = 0.2;
plate_len = 1.0;
Point(1) = {0, 0, 0};
Point(2) = {hole_radius, 0, 0};
Point(3) = {plate_len, 0, 0};
Point(4) = {plate_len, plate_len, 0};
Point(5) = {0, plate_len, 0};
Point(6) = {0, hole_radius, 0};




//+
Circle(1) = {2, 1, 6};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 5};
//+
Line(5) = {5, 6};
//+
Curve Loop(1) = {5, -1, 2, 3, 4};
//+
Plane Surface(1) = {1};

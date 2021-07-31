airfoil_cl = 0.050000;
wall_cl = 0.250000;
Lx_plus = 1.500000;
Lx_minus = -0.500000;
Ly_plus = 0.500000;
Ly_minus = -0.500000;

Point(1) = {0.000000, 0.000000, 0.000000, airfoil_cl};
Point(2) = {0.006819, 0.014185, 0.000000, airfoil_cl};
Point(3) = {0.027091, 0.027121, 0.000000, airfoil_cl};
Point(4) = {0.060263, 0.038445, 0.000000, airfoil_cl};
Point(5) = {0.105430, 0.047719, 0.000000, airfoil_cl};
Point(6) = {0.161359, 0.054541, 0.000000, airfoil_cl};
Point(7) = {0.226526, 0.058654, 0.000000, airfoil_cl};
Point(8) = {0.299152, 0.060007, 0.000000, airfoil_cl};
Point(9) = {0.377257, 0.058770, 0.000000, airfoil_cl};
Point(10) = {0.458710, 0.055295, 0.000000, airfoil_cl};
Point(11) = {0.541290, 0.050047, 0.000000, airfoil_cl};
Point(12) = {0.622743, 0.043533, 0.000000, airfoil_cl};
Point(13) = {0.700848, 0.036253, 0.000000, airfoil_cl};
Point(14) = {0.773474, 0.028672, 0.000000, airfoil_cl};
Point(15) = {0.838641, 0.021228, 0.000000, airfoil_cl};
Point(16) = {0.894570, 0.014344, 0.000000, airfoil_cl};
Point(17) = {0.939737, 0.008430, 0.000000, airfoil_cl};
Point(18) = {0.972909, 0.003870, 0.000000, airfoil_cl};
Point(19) = {0.993181, 0.000987, 0.000000, airfoil_cl};
Point(20) = {1.000000, -0.000000, 0.000000, airfoil_cl};
Point(21) = {0.993181, -0.000987, 0.000000, airfoil_cl};
Point(22) = {0.972909, -0.003870, 0.000000, airfoil_cl};
Point(23) = {0.939737, -0.008430, 0.000000, airfoil_cl};
Point(24) = {0.894570, -0.014344, 0.000000, airfoil_cl};
Point(25) = {0.838641, -0.021228, 0.000000, airfoil_cl};
Point(26) = {0.773474, -0.028672, 0.000000, airfoil_cl};
Point(27) = {0.700848, -0.036253, 0.000000, airfoil_cl};
Point(28) = {0.622743, -0.043533, 0.000000, airfoil_cl};
Point(29) = {0.541290, -0.050047, 0.000000, airfoil_cl};
Point(30) = {0.458710, -0.055295, 0.000000, airfoil_cl};
Point(31) = {0.377257, -0.058770, 0.000000, airfoil_cl};
Point(32) = {0.299152, -0.060007, 0.000000, airfoil_cl};
Point(33) = {0.226526, -0.058654, 0.000000, airfoil_cl};
Point(34) = {0.161359, -0.054541, 0.000000, airfoil_cl};
Point(35) = {0.105430, -0.047719, 0.000000, airfoil_cl};
Point(36) = {0.060263, -0.038445, 0.000000, airfoil_cl};
Point(37) = {0.027091, -0.027121, 0.000000, airfoil_cl};
Point(38) = {0.006819, -0.014185, 0.000000, airfoil_cl};

// Define bounding box corners for fluid domain 
Point(39) = {Lx_minus, Ly_minus, 0.0000, wall_cl};
Point(40) = {Lx_minus, Ly_plus, 0.0000, wall_cl};
Point(41) = {Lx_plus, Ly_plus, 0.0000, wall_cl};
Point(42) = {Lx_plus, Ly_minus, 0.0000, wall_cl};

// Make line loops for the fluid domain
Line(39) = {39, 40};
Line(40) = {40, 41};
Line(41) = {41, 42};
Line(42) = {42, 39};

// Combine the fluid domain lines
Line Loop(101) = {39, 40, 41, 42};

// combine foil points into a bspline 
BSpline(1) = {1:20};
BSpline(2) = {20:38,1};
Line Loop(102) = {1, 2};

Plane Surface(201) = {101, 102};

Extrude {0, 0, 1} {
 Surface{201};
 Layers{1};
 Recombine;
}
// Define the physical surfaces -> values are obtained via the GUI 
Physical Surface("overset") = {212, 216, 220, 224}; 
Physical Surface("symFront") = {233}; 
Physical Surface("symBack") = {201}; 
Physical Surface("wingTop") = {232}; 
Physical Surface("wingBottom") = {228}; 

// Define the internal volume - usually 1 but obtained from GUI
Physical Volume("internal") = {1};

//Define Boundary Layer
Field[1].EdgesList = {5, 227, 208, -226, 206, 6};
Field[1].NodesList = {1,1021, 1022, 160};
Field[1].hfar = 0.01;
Field[1].hwall_n = 0.0005;
Field[1].thickness = 0.02;
Field[1].ratio = 1.11;
Field[1].AnisoMax = 10;
Field[1].Quads = 1;
Field[1].IntersectMetrics = 0;
BoundaryLayer Field = 1;

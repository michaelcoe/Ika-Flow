ss_cl =  0.002500;
ps_cl =  0.002500;
le_cl = 0.001000;
te_cl = 0.001000;
wall_cl =  0.030000;
extrude_length =  1.000000;

Lx_plus = 1.500000;
Lx_minus = -0.500000;
Ly_plus = 0.500000;
Ly_minus = -0.500000;

back_plane = -0.500000;

Point(1) = {0.000000, 0.000000, back_plane, le_cl};
Point(2) = {0.001027, 0.005653, back_plane, le_cl};
Point(3) = {0.004105, 0.011181, back_plane, ss_cl};
Point(4) = {0.009220, 0.016567, back_plane, ss_cl};
Point(5) = {0.016353, 0.021785, back_plane, ss_cl};
Point(6) = {0.025472, 0.026806, back_plane, ss_cl};
Point(7) = {0.036542, 0.031592, back_plane, ss_cl};
Point(8) = {0.049516, 0.036106, back_plane, ss_cl};
Point(9) = {0.064341, 0.040308, back_plane, ss_cl};
Point(10) = {0.080956, 0.044160, back_plane, ss_cl};
Point(11) = {0.099293, 0.047626, back_plane, ss_cl};
Point(12) = {0.119277, 0.050678, back_plane, ss_cl};
Point(13) = {0.140825, 0.053295, back_plane, ss_cl};
Point(14) = {0.163850, 0.055465, back_plane, ss_cl};
Point(15) = {0.188255, 0.057186, back_plane, ss_cl};
Point(16) = {0.213942, 0.058469, back_plane, ss_cl};
Point(17) = {0.240804, 0.059340, back_plane, ss_cl};
Point(18) = {0.268731, 0.059833, back_plane, ss_cl};
Point(19) = {0.297608, 0.059999, back_plane, ss_cl};
Point(20) = {0.327317, 0.059888, back_plane, ss_cl};
Point(21) = {0.357736, 0.059504, back_plane, ss_cl};
Point(22) = {0.388740, 0.058839, back_plane, ss_cl};
Point(23) = {0.420200, 0.057893, back_plane, ss_cl};
Point(24) = {0.451988, 0.056667, back_plane, ss_cl};
Point(25) = {0.483974, 0.055169, back_plane, ss_cl};
Point(26) = {0.505475, 0.054019, back_plane, ss_cl};
Point(27) = {0.526976, 0.052756, back_plane, ss_cl};
Point(28) = {0.548477, 0.051384, back_plane, ss_cl};
Point(29) = {0.569979, 0.049906, back_plane, ss_cl};
Point(30) = {0.591480, 0.048324, back_plane, ss_cl};
Point(31) = {0.612981, 0.046641, back_plane, ss_cl};
Point(32) = {0.634482, 0.044861, back_plane, ss_cl};
Point(33) = {0.655983, 0.042986, back_plane, ss_cl};
Point(34) = {0.677484, 0.041019, back_plane, ss_cl};
Point(35) = {0.698985, 0.038964, back_plane, ss_cl};
Point(36) = {0.720486, 0.036822, back_plane, ss_cl};
Point(37) = {0.741987, 0.034598, back_plane, ss_cl};
Point(38) = {0.763488, 0.032293, back_plane, ss_cl};
Point(39) = {0.784989, 0.029912, back_plane, ss_cl};
Point(40) = {0.806490, 0.027456, back_plane, ss_cl};
Point(41) = {0.827991, 0.024929, back_plane, ss_cl};
Point(42) = {0.849492, 0.022334, back_plane, ss_cl};
Point(43) = {0.870994, 0.019673, back_plane, ss_cl};
Point(44) = {0.892495, 0.016950, back_plane, ss_cl};
Point(45) = {0.913996, 0.014167, back_plane, ss_cl};
Point(46) = {0.935497, 0.011328, back_plane, ss_cl};
Point(47) = {0.956998, 0.008435, back_plane, te_cl};
Point(48) = {0.978499, 0.005491, back_plane, te_cl};
Point(49) = {1.000000, 0.002500, back_plane, te_cl};
Point(50) = {1.000000, 0.000000, back_plane, te_cl};
Point(51) = {1.000000, -0.002500, back_plane, te_cl};
Point(52) = {0.978499, -0.005491, back_plane, te_cl};
Point(53) = {0.956998, -0.008435, back_plane, te_cl};
Point(54) = {0.935497, -0.011328, back_plane, ps_cl};
Point(55) = {0.913996, -0.014167, back_plane, ps_cl};
Point(56) = {0.892495, -0.016950, back_plane, ps_cl};
Point(57) = {0.870994, -0.019673, back_plane, ps_cl};
Point(58) = {0.849492, -0.022334, back_plane, ps_cl};
Point(59) = {0.827991, -0.024929, back_plane, ps_cl};
Point(60) = {0.806490, -0.027456, back_plane, ps_cl};
Point(61) = {0.784989, -0.029912, back_plane, ps_cl};
Point(62) = {0.763488, -0.032293, back_plane, ps_cl};
Point(63) = {0.741987, -0.034598, back_plane, ps_cl};
Point(64) = {0.720486, -0.036822, back_plane, ps_cl};
Point(65) = {0.698985, -0.038964, back_plane, ps_cl};
Point(66) = {0.677484, -0.041019, back_plane, ps_cl};
Point(67) = {0.655983, -0.042986, back_plane, ps_cl};
Point(68) = {0.634482, -0.044861, back_plane, ps_cl};
Point(69) = {0.612981, -0.046641, back_plane, ps_cl};
Point(70) = {0.591480, -0.048324, back_plane, ps_cl};
Point(71) = {0.569979, -0.049906, back_plane, ps_cl};
Point(72) = {0.548477, -0.051384, back_plane, ps_cl};
Point(73) = {0.526976, -0.052756, back_plane, ps_cl};
Point(74) = {0.505475, -0.054019, back_plane, ps_cl};
Point(75) = {0.483974, -0.055169, back_plane, ps_cl};
Point(76) = {0.451988, -0.056667, back_plane, ps_cl};
Point(77) = {0.420200, -0.057893, back_plane, ps_cl};
Point(78) = {0.388740, -0.058839, back_plane, ps_cl};
Point(79) = {0.357736, -0.059504, back_plane, ps_cl};
Point(80) = {0.327317, -0.059888, back_plane, ps_cl};
Point(81) = {0.297608, -0.059999, back_plane, ps_cl};
Point(82) = {0.268731, -0.059833, back_plane, ps_cl};
Point(83) = {0.240804, -0.059340, back_plane, ps_cl};
Point(84) = {0.213942, -0.058469, back_plane, ps_cl};
Point(85) = {0.188255, -0.057186, back_plane, ps_cl};
Point(86) = {0.163850, -0.055465, back_plane, ps_cl};
Point(87) = {0.140825, -0.053295, back_plane, ps_cl};
Point(88) = {0.119277, -0.050678, back_plane, ps_cl};
Point(89) = {0.099293, -0.047626, back_plane, ps_cl};
Point(90) = {0.080956, -0.044160, back_plane, ps_cl};
Point(91) = {0.064341, -0.040308, back_plane, ps_cl};
Point(92) = {0.049516, -0.036106, back_plane, ps_cl};
Point(93) = {0.036542, -0.031592, back_plane, ps_cl};
Point(94) = {0.025472, -0.026806, back_plane, ps_cl};
Point(95) = {0.016353, -0.021785, back_plane, ps_cl};
Point(96) = {0.009220, -0.016567, back_plane, ps_cl};
Point(97) = {0.004105, -0.011181, back_plane, ps_cl};
Point(98) = {0.001027, -0.005653, back_plane, le_cl};

// Define bounding box corners for fluid domain 
Point(99) = {Lx_minus, Ly_minus, back_plane, wall_cl};
Point(100) = {Lx_minus, Ly_plus, back_plane, wall_cl};
Point(101) = {Lx_plus, Ly_plus, back_plane, wall_cl};
Point(102) = {Lx_plus, Ly_minus, back_plane, wall_cl};

// Make line loops for the fluid domain
Line(99) = {99, 100};
Line(100) = {100, 101};
Line(101) = {101, 102};
Line(102) = {102, 99};

// Combine the fluid domain lines
Line Loop(101) = { 99, 100, 101, 102}; 

// combine foil points into a bspline 
BSpline(5) = {1:50}; 
BSpline(6) = {50:98,1}; 
Line Loop(102) = {5, 6}; 

Plane Surface(201) = {101, 102}; 

//Recombine Surface {201}; 

Extrude {0, 0, extrude_length} {
 Surface{201};
 Layers{1};
 Recombine;
}
// Define the physical surfaces -> values are obtained via the GUI 
Physical Surface("inlet") = {212}; 
Physical Surface("outlet") = {220}; 
Physical Surface("symFront") = {233}; 
Physical Surface("symBack") = {201}; 
Physical Surface("topAndBottom") = {216, 224}; 
Physical Surface("wingT") = {232}; 
Physical Surface("wingB") = {228}; 

// Define the internal volume - usually 1 but obtained from GUI 
Physical Volume("internal") = {1};

//Define the boundary layer 
Field[1] = BoundaryLayer; 
Field[1].EdgesList = {5, 6}; 
Field[1].hfar = 0.00177; 
Field[1].hwall_n = 3.75e-05; 
Field[1].thickness = 0.0162; 
Field[1].ratio = 1.12; 
Field[1].AnisoMax = 10; 
Field[1].Quads = 1; 
Field[1].IntersectMetrics = 0; 
BoundaryLayer Field = 1; 

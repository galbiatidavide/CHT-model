// 2D Conjugate Heat Transfer Channel Mesh
// Parameters
L = 0.05;
H_top = 0.004;
H_bot = 0.001;
H_tot = H_top + H_bot;

// ----- Geometry points -----
Point(1) = {0, H_bot, 0, 1.0};
Point(2) = {L, H_bot, 0, 1.0};
Point(3) = {L, H_tot, 0, 1.0};
Point(4) = {0, H_tot, 0, 1.0};
Point(5) = {0, 0,     0, 1.0};
Point(6) = {L, 0,     0, 1.0};

// ----- Interface line (shared) -----
Line(1) = {1, 2};  // interface

// Top rectangle (fluid)
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line Loop(10) = {1,2,3,4};
Plane Surface(200) = {10};

// Bottom rectangle (solid)
Line(5) = {5, 6};
Line(6) = {6, 2};
Line(7) = {1, 5};
Line Loop(20) = {5,6,-1,7}; 
Plane Surface(100) = {20};

// ----- Physical groups -----
Physical Surface(200) = {200};   // Fluid
Physical Surface(100) = {100};   // Solid

Physical Line(5) = {1};         // Interface
Physical Line(2) = {5};         // Bottom
Physical Line(3) = {2,6};       // Right (outlet)
Physical Line(4) = {3};         // Top
Physical Line(1) = {4,7};       // Left (inlet)


// ----- Mesh refinement (separate interface + inlet) -----

//+
Transfinite Curve {3, 1, 5} = 150 Using Progression 1;
//+
Transfinite Curve {4, 2, 6, 7} = 15 Using Progression 1;
//+

//+
Transfinite Surface {100};
//+
Transfinite Surface {200};
//+

//+
Recombine Surface {100, 200};



// Distance from interface
Field[1] = Distance;
Field[1].CurvesList = {1};       
Field[1].NumPointsPerCurve = 100;

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 4.0e-4;  // 0.40 mm
Field[2].SizeMax = 2.0e-3;  // 2.00 mm
Field[2].DistMin = 5.0e-4;  // fine zone 0.5 mm near interface
Field[2].DistMax = 2.0e-3;  // transition until 2 mm

// Distance from inlet
Field[3] = Distance;
Field[3].CurvesList = {4, 7};    // inlet lines (fluid+solid)
Field[3].NumPointsPerCurve = 100;

Field[4] = Threshold;
Field[4].InField = 3;
Field[4].SizeMin = 4.0e-4;  // 0.40 mm
Field[4].SizeMax = 2.0e-3;  // 2.00 mm
Field[4].DistMin = 2.0e-3;  // fine zone extends 2 mm into the channel
Field[4].DistMax = 3.75e-2;  // transition until 10 mm


Field[5] = Min;
Field[5].FieldsList = {2, 4};

Background Field = 5;

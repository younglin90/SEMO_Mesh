// Gmsh project created on Thu Apr 28 12:23:13 2011
Point(1) = {0, 0, 0, 1.0};
Point(2) = {1, 0, 0, 1.0};
Point(3) = {1, 1, 0, 1.0};
Point(4) = {0, 1, 0, 1.0};
Line(1) = {4, 3};
Line(2) = {3, 2};
Line(3) = {2, 1};
Line(4) = {1, 4};
Line Loop(5) = {1, 2, 3, 4};
Plane Surface(6) = {5};
Extrude {0, 0, 0.01} {
  Surface{6};
  Layers{1};
  Recombine;
}
Physical Surface("bottom") = {23};
Physical Surface("left") = {27};
Physical Surface("top") = {15};
Physical Surface("right") = {19};
Physical Surface("front") = {28};
Physical Surface("back") = {6};
Physical Volume("internal") = {1};

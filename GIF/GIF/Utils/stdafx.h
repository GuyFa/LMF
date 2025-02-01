#pragma once

//////  CGAL  //////
#include <CGAL/Timer.h>
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Polyhedron_3.h>
#include <CGAL/Polyhedron_incremental_builder_3.h>
#include <CGAL/Arr_segment_traits_2.h>
#include <CGAL/Arrangement_2.h>
#include <CGAL/Arr_landmarks_point_location.h>
#include <CGAL/Arr_extended_dcel.h>
#include <CGAL/Boolean_set_operations_2.h> //needed for testing intersection of two polygons (do_intersect)
#include <CGAL/Arr_default_overlay_traits.h>
#include <CGAL/Sweep_line_2_algorithms.h>
#include <CGAL/Polygon_2.h>





#undef short2
#undef short3
#undef long2
#undef long3
#undef int2
#undef int3
#undef float2
#undef float3
#undef double2
#undef double3
#undef double4

using namespace std;

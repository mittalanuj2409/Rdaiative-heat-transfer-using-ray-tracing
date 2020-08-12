//#include<bits/stdc++.h>
#include<iostream>
#include<string>
#include<vector>
#include<sstream>

using namespace std;

/*https://wiki.eclipse.org/IGES_file_Specification*/
class point {
public:
	float x, y, z;
};

class cir_arc {
public:
	int de_no, t_matrix, form;
	/*type 100*/
	float z, x_center, y_center, x_start, y_start, x_end, y_end;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~cir_arc() {};
};

class comp_curve {
	/*type 102*/
public:
	int n, de_no, t_matrix, form;
	float* list = (float*)malloc(1e4 * sizeof(float)); //*n
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~comp_curve() {};
};

class conic_arc {
public:
	int de_no, t_matrix, form;
	/*type 104*/
	float a, b, c, d, e, f, x_start, y_start, z_start, x_end, y_end, z_end;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~conic_arc() {};
};

class plane {
public:
	int de_no, t_matrix, form;
	/*type 108*/
	float a, b, c, d, x, y, z, size;
	int pointer;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~plane() {};
};

class line {
public:
	int de_no, t_matrix, form;
	/*type 110*/
	float x_start, y_start, z_start, x_end, y_end, z_end;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~line() {};
};

class par_spline_curve {
public:
	int de_no, t_matrix, form;
	/*type 112*/
	int type, k, dim, n;
	float* t_points = (float*)malloc(sizeof(float) * 1e4); //*n
	class segment {
	public:
		float ax, bx, cx, dx, ay, by, cy, dy, az, bz, cz, dz;
	};
	segment* list = (class segment*) malloc(sizeof(class segment) * 1e4); //*n
	float xn, xn1, xn2, xn3, yn, yn1, yn2, yn3, zn, zn1, zn2, zn3;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~par_spline_curve() {};
};

class par_spline_sur {
public:
	int de_no, t_matrix, form;
	/*type 114*/
	int type, k, m, n;
	float* tu_points = (float*)malloc(sizeof(float) * 1e4); //*m
	float* tv_points = (float*)malloc(sizeof(float) * 1e4); //*n
	float* p_ij = (float*)malloc(sizeof(float) * 1e4); //*m*n
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~par_spline_sur() {};
};

class point_ {
public:
	int de_no, t_matrix, form;
	/*type 116*/
	float x, y, z;
	int pointer;
	~point_() {};
};

class ruled_sur {
public:
	int de_no, t_matrix, form;
	/*type 118*/
	int pointer1, pointer2, dir_flag, dev_flag;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~ruled_sur() {};
};

class sur_revolution {
public:
	int de_no, t_matrix, form;
	/*type 120*/
	int pointer1, pointer2;
	float sa, ea;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~sur_revolution() {};
};

class tab_cylinder {
public:
	int de_no, t_matrix, form;
	/*type 122*/
	int pointer;
	float lx, ly, lz;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~tab_cylinder() {};
};

class direction {
public:
	int de_no, t_matrix, form;
	/*type 123*/
	float x, y, z;
	~direction() {};
};

class trans_matrix {
public:
	int de_no, t_matrix, form;
	/*type 124*/
	float r11, r12, r13, t1, r21, r22, r23, t2, r31, r32, r33, t3;
	~trans_matrix() {};
};

class rat_spline_curve {
public:
	int de_no, t_matrix, form;
	/*type 126*/
	int k, m, f1, f2, f3, f4;
	/*
	k: number of control points minus one / upper index sum
	m: degree of B-spline basis function
	*/
	float* knot = (float*)malloc(sizeof(float) * 1e4); //*(2+k+m)
	float* weight = (float*)malloc(sizeof(float) * 1e4); //*(1+k)
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	point* control_points = (class point*) malloc(sizeof(class point) * 1e4); //*(1+k)
	float v0, v1, xn, yn, zn;
	~rat_spline_curve() {};
};

class rat_spline_sur {
public:
	int de_no, t_matrix, form;
	/*type 128*/
	int k1, k2, m1, m2, f1, f2, f3, f4, f5;
	/*
	k1: number of control points in first direction minus one / upper index of outer sum
	k2: number of control points in second direction minus one / upper index of inner sum
	m1: degree of B-spline basis functions in first direction
	m2: degree of B-spline basis functions in second direction
	*/
	float* knot1 = (float*)malloc(sizeof(float) * 1e4); //(2+k1+m1));
	float* knot2 = (float*)malloc(sizeof(float) * 1e4); //(2+k2+m2));
	float* weight = (float*)malloc(sizeof(float) * 1e4);//((1+k1)*(1+k2)));
	point* control_points = (class point*) calloc(1e4, sizeof(class point)); //*(1+k1)*(1+k2)
	point* actual_points = (class point*) calloc(1e4, sizeof(class point));
	float u0, u1, v0, v1;
	~rat_spline_sur() {};
};

class offset_curve {
public:
	int de_no, t_matrix, form;
	/*type 130*/
	int pointer1, f1, pointer2, dim, f2;
	float d1, td1, d2, td2, x, y, z, tt1, tt2;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~offset_curve() {};
};

class offset_sur {
public:
	int de_no, t_matrix, form;
	/*type 140*/
	float nx, ny, nz, d;
	int pointer;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~offset_sur() {};
};

class boundary {
public:
	int de_no, t_matrix, form;
	/*type 141*/
	int type, pref, pointer1, n;
	class curve {
	public:
		int pointer, f, k;
		float* p = (float*)malloc(sizeof(float) * 1e4); //*k
	};
	curve* curves = (class curve*) malloc(sizeof(class curve) * 1e4); //*n
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~boundary() {};
};

class curve_on_sur {
public:
	int de_no, t_matrix, form;
	/*type 142*/
	int flag, pointer1, pointer2, pointer3, rep;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~curve_on_sur() {};
};

class bounded_sur {
public:
	int de_no, t_matrix, form;
	/*type 143*/
	int type, pointer1, n;
	int* list = (int*)malloc(1e4 * sizeof(int)); //*n
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~bounded_sur() {};
};

class trimmed_sur {
public:
	int de_no, t_matrix, form;
	/*type 144*/
	int pointer1, flag, n, pointer2;
	int* list = (int*)malloc(1e4 * sizeof(int)); //*n
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~trimmed_sur() {};
};

class block {
public:
	int de_no, t_matrix, form;
	/*type 150. CSG block object*/
	float length_x, length_y, length_z, x_corner, y_corner, z_corner, x_i, x_j, x_k, z_i, z_j, z_k;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~block() {};
};

class wedge {
public:
	int de_no, t_matrix, form;
	/*type 152. CSG wedge object*/
	float length_x, length_y, length_z, tlx, x_corner, y_corner, z_corner, x_i, x_j, x_k, z_i, z_j, z_k;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~wedge() {};
};

class cylinder {
public:
	int de_no, t_matrix, form;
	/*type 154. CSG cylinder*/
	float h, r, x_center, y_center, z_center, i, j, k;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~cylinder() {};
};

class cone {
public:
	int de_no, t_matrix, form;
	/*type 156. CSG cone*/
	float h, r_large, r_small, x_center, y_center, z_center, i, j, k;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~cone() {};
};

class sphere {
public:
	int de_no, t_matrix, form;
	/*type 158. CSG sphere*/
	float r, x, y, z;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~sphere() {};
};

class torus {
public:
	int de_no, t_matrix, form;
	/*type 160. CSG torus*/
	float r1, r2, x_center, y_center, z_center, i, j, k;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~torus() {};
};

class solid_rev {
public:
	int de_no, t_matrix, form;
	/*type 162*/
	int pointer;
	float f, x, y, z, i, j, k;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~solid_rev() {};
};

class solid_ext {
public:
	int de_no, t_matrix, form;
	/*type 164*/
	int pointer;
	float l, i, j, k;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~solid_ext() {};
};

class ellipsoid {
public:
	int de_no, t_matrix, form;
	/*type 168. CSG ellipsoid*/
	float length_x, length_y, length_z, x_center, y_center, z_center, x_i, x_j, x_k, z_i, z_j, z_k;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~ellipsoid() {};
};

class bool_tree {
public:
	int de_no, t_matrix, form;
	/*type 180*/
	int n;
	int* pointer = (int*)malloc(1e4 * sizeof(int)); //*n
	~bool_tree() {};
};

class solid {
public:
	int de_no, t_matrix, form;
	/*type 186*/
	int pointer, flag, n;
	class shell {
	public:
		int pointer, vflag;
	};
	shell* list = (class shell*) malloc(1e4 * sizeof(class shell)); //*n
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~solid() {};
};

class plane_sur {
public:
	int de_no, t_matrix, form;
	/*type 190*/
	int pointer1, pointer2, pointer3;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~plane_sur() {};
};

class cyl_sur {
public:
	int de_no, t_matrix, form;
	/*type 192*/
	int point_, axis, ref;
	float r;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~cyl_sur() {};
};

class con_sur {
public:
	int de_no, t_matrix, form;
	/*type 194*/
	int point_, axis, ref;
	float r, ang;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~con_sur() {};
};

class sph_sur {
public:
	int de_no, t_matrix, form;
	/*type 196*/
	int point_, axis, ref;
	float r;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~sph_sur() {};
};

class tor_sur {
public:
	int de_no, t_matrix, form;
	/*type 198*/
	int point_, axis, ref;
	float r1, r2;
	point* actual_points = (class point*) malloc(sizeof(class point) * 1e4);
	~tor_sur() {};
};

class subfigure {
public:
	int de_no, t_matrix, form;
	/*type 308*/
	int depth, n;
	string name;
	int* pointer = (int*)malloc(1e4 * sizeof(int)); //*n
	~subfigure() {};
};

class color {
public:
	int de_no, t_matrix, form;
	/*type 314*/
	float red, blue, green;
	string name;
	~color() {};
};

class subfig_instance {
public:
	int de_no, t_matrix, form;
	/*type 408*/
	int sd;
	float x, y, z, s;
	~subfig_instance() {};
};

class vertex_list {
public:
	int de_no, t_matrix, form;
	/*type 502*/
	int n;
	point* list = (class point*) malloc(1e4 * sizeof(class point)); //*n
	~vertex_list() {};
};

class edge_list {
public:
	int de_no, t_matrix, form;
	/*type 504*/
	int n;
	class edge {
	public:
		int pointer1, start_pointer, start, end_pointer, end;
	};
	edge* list = (class edge*) malloc(1e4 * sizeof(class edge)); //*n
	~edge_list() {};
};

class loop {
public:
	int de_no, t_matrix, form;
	/*type 508*/
	int n;
	class edge_or_vertex {
	public:
		int type, pointer, index, flag, k;
		//		if(k!=0) {
		class space_curve {
		public:
			int iso_flag, pointer;
		};
		space_curve* list1 = (class space_curve*) malloc(1e4 * sizeof(class space_curve)); //*k
//		}
	};
	edge_or_vertex* list2 = (class edge_or_vertex*) malloc(1e4 * sizeof(class edge_or_vertex)); //*n
	~loop() {};
};

class face {
public:
	int de_no, t_matrix, form;
	/*type 510*/
	int pointer, n, flag;
	int* loop_pointer = (int*)malloc(1e4 * sizeof(int)); //*n
	~face() {};
};

class shell_ {
public:
	int de_no, t_matrix, form;
	/*type 514*/
	int n;
	class faces {
	public:
		int pointer, flag;
	};
	faces* list = (class faces*) malloc(1e4 * sizeof(class faces)); //*n
	~shell_() {};
};

vector <cir_arc> cir_arc_vec;
vector <comp_curve> comp_curve_vec;
vector <conic_arc> conic_arc_vec;
vector <plane> plane_vec;
vector <line> line_vec;
vector <par_spline_curve> par_spline_curve_vec;
vector <par_spline_sur> par_spline_sur_vec;
vector <point_> point_vec;
vector <ruled_sur> ruled_sur_vec;
vector <sur_revolution> sur_revolution_vec;
vector <tab_cylinder> tab_cylinder_vec;
vector <direction> direction_vec;
vector <trans_matrix> trans_matrix_vec;
vector <rat_spline_curve> rat_spline_curve_vec;
vector <rat_spline_sur> rat_spline_sur_vec;
vector <offset_curve> offset_curve_vec;
vector <offset_sur> offset_sur_vec;
vector <boundary> boundary_vec;
vector <curve_on_sur> curve_on_sur_vec;
vector <bounded_sur> bounded_sur_vec;
vector <trimmed_sur> trimmed_sur_vec;
vector <block> block_vec;
vector <wedge> wedge_vec;
vector <cylinder> cylinder_vec;
vector <cone> cone_vec;
vector <sphere> sphere_vec;
vector <torus> torus_vec;
vector <solid_rev> solid_rev_vec;
vector <solid_ext> solid_ext_vec;
vector <ellipsoid> ellipsoid_vec;
vector <bool_tree> bool_tree_vec;
vector <solid> solid_vec;
vector <plane_sur> plane_sur_vec;
vector <cyl_sur> cyl_sur_vec;
vector <con_sur> con_sur_vec;
vector <sph_sur> sph_sur_vec;
vector <tor_sur> tor_sur_vec;
vector <subfigure> subfigure_vec;
vector <color> color_vec;
vector <subfig_instance> subfig_instance_vec;
vector <vertex_list> vertex_list_vec;
vector <edge_list> edge_list_vec;
vector <loop> loop_vec;
vector <face> face_vec;
vector <shell_> shell_vec;

vector <int> t_matrix_pointer;
vector <int> form_number;

vector <point*> all; //all the actual points are stored here
vector <point*> source; //source actual points are stored here
vector <point*> receiver; //receiver actual points are stored here



#include<iostream>
#include<string>
#include<vector>
#include<sstream>
#include<math.h>
#include"classes.h"
//#include"functions.cpp"
using namespace std;

template <class T>

bool confirm(vector <T> a){
	int seq;
	char c;
	cout<<"Enter a sequence no. to select an object\n":
	cin>>seq;
	if(seq<1 || seq>a.size()) return FALSE;
	system("CLS");
	cout<<"\n The representative actual points on selected objects are:\n";
	for(int i = 0; i < Np_max; i++){
		cout<<"( "<<a[seq-1].actual_points[i].x<<", "<<a[seq-1].actual_points[i].y<<", "<<a[seq-1].actual_points[i].z<<" )"<<endl;
	}
	cout<<"\nDo you want to continue? (y/n)\n":
	cin>>c;
	if(c != 'y' || c != 'Y') return FALSE;
	source.push_back(a[seq-1].actual_points);
	return TRUE;

}

void type_disp_0() {}

void cir_arc_disp() {
	/*type 100. index 1*/
	cout<<"\n Circular Arc (anti-clock)\n\n";
	for (int i = 0; i < cir_arc_vec.size(), i++){
		cout<<" sequence no: "<<i+1<<endl;
		cout<<" de_no      : "<<cir_arc_vec[i].de_no<<endl;
		cout<<" center     : ("<<cir_arc_vec[i].x_center<<","<<cir_arc_vec[i].y_center<<","<<cir_arc_vec[i].z<<")"<<endl;
		cout<<" Start      : ("<<cir_arc_vec[i].x_start<<","<<cir_arc_vec[i].y_start<<","<<cir_arc_vec[i].z<<")"<<endl;
		cout<<" End        : ("<<cir_arc_vec[i].x_end<<","<<cir_arc_vec[i].y_end<<","<<cir_arc_vec[i].z<<")"<<endl;
		cout<<"\n";
	}
	if(!confirm<cir_arc>(cir_arc_vec)){
		cir_arc_disp();
	}
}

void comp_curve_disp() {
	/*type 102. index 2*/

}

void conic_arc_disp() {
	/*type 104. index 3*/

}

void type_disp_106() {}

void plane_disp() {
	/*type 108. index 5*/

}

void line_disp() {
	/*type 110. index 6*/

}

void par_spline_curve_disp() {
	/*type 112. index 7*/
	
}

void par_spline_sur_disp() {
	/*type 114. index 8*/
	
}

void point_disp() {
	/*type 116. index 9*/
	
}

void ruled_sur_disp() {
	/*type 118. index 10*/
	
}

void sur_revolution_disp() {
	/*type 120. index 11*/
	
}

void tab_cylinder_disp() {
	/*type 122. index 12*/
	
}

void direction_disp() {
	/*type 123. index 13*/
	
}

void trans_matrix_disp() {
	/*typr 124. index 14*/
	
}

void type_disp_125() {}

void rat_spline_curve_disp() {
	/*type 126. index 16*/
	
}

void rat_spline_sur_disp() {
	/*type 128. index 17*/
	
}

void offset_curve_disp() {
	/*type 130. index 18*/
	
}

void type_disp_132() {}

void type_disp_134() {}

void type_disp_136() {}

void type_disp_138() {}

void offset_sur_disp() {
	/*type 140. index 23*/
	
}

void boundary_disp() {
	/*type 141. index 24*/
	
}

void curve_on_sur_disp() {
	/*type 142. index 25*/
	
}

void bounded_sur_disp() {
	/*type 143. index 26*/
	
}

void trimmed_sur_disp() {
	/*type 144. index 27*/
	
}

void type_disp_146() {}

void type_disp_148() {}

void block_disp() {
	/*typr 150. index 30*/
	
}

void wedge_disp() {
	/*typr 152. index 31*/
	
}

void cylinder_disp() {
	/*typr 154. index 32*/
	
}

void cone_disp() {
	/*typr 156. index 33*/
	
}

void sphere_disp() {
	/*type 158. index 34*/
	
}

void torus_disp() {
	/*typr 160. index 35*/
	
}

void solid_rev_disp() {
	/*type 162. index 36*/
	
}

void solid_ext_disp() {
	/*type 164. index 37*/
	
}

void ellipsoid_disp() {
	/*type 168. index 38*/
	
}

void bool_tree_disp() {
	/*type 180. index 39*/
	
}

void type_disp_182() {}

void type_disp_184() {}

void solid_disp() {
	/*type 186. index 42*/
	
}

void plane_sur_disp() {
	/*type 190. index 43*/
	
}


void cyl_sur_disp() {
	/*type 192. index 44*/
	
}

void con_sur_disp() {
	/*type 194. index 45*/
	
}

void sph_sur_disp() {
	/*type 196. index 46*/
	
}

void tor_sur_disp() {
	/*type 198. index 47*/
	
}

void type_disp_202() {}

void type_disp_204() {}

void type_disp_206() {}

void type_disp_208() {}

void type_disp_210() {}

void type_disp_212() {}

void type_disp_213() {}

void type_disp_214() {}

void type_disp_216() {}

void type_disp_218() {}

void type_disp_220() {}

void type_disp_222() {}

void type_disp_228() {}

void type_disp_230() {}

void type_disp_302() {}

void type_disp_304() {}

void type_disp_306() {}

void type_disp_308() {} //dummy

void subfigure_disp() {
	/*type 308. index 64*/
	
}

void type_disp_310() {}

void type_disp_312() {}

void type_disp_314() {} //dummy

void color_disp() {
	/*type 314. index 68*/
	
}

void type_disp_316() {}

void type_disp_320() {}

void type_disp_322() {}

void type_disp_402() {}

void type_disp_404() {}

void type_disp_406() {}

void subfig_instance_disp() {
	/*type 408. index 75*/
	
}

void type_disp_410() {}

void type_disp_412() {}

void type_disp_414() {}

void type_disp_416() {}

void type_disp_418() {}

void type_disp_420() {}

void type_disp_422() {}

void type_disp_430() {}

void vertex_list_disp() {
	/*type 502. index 84*/
	
}

void edge_list_disp() {
	/*type 504. index 85*/
	
}

void loop_disp() {
	/*type 508. index 86*/
	
}

void face_disp() {
	/*type 510. index 87*/
	
}

void shell_disp() {
	/*type 514. index 88*/
	
}


int sur_select_disp(){
	int type;
	contents();
	cout<<"\nEnter type no. for required surface:\n";
	cin>>type;
	system("CLS");
	if (index(type) == -1)	cout << "error: unknown type no";
	else (*disp_list[index(type)])();
	
}

int index(int type) {
	/*provide index function in function pointer array*/
	int i, key = 0;
	for (i = 0; i < 89; i++) {
		if (type == entity[i]) {
			key = 1;
			break;
		}
	}
	if (key == 1)	return i;
	return -1;
}

void (*disp_list[])() =
{ type_disp_0,
cir_arc_disp, comp_curve_disp, conic_arc_disp, type_disp_106, plane_disp,
line_disp, par_spline_curve_disp, par_spline_sur_disp, point_disp, ruled_sur_disp,
sur_revolution_disp, tab_cylinder_disp, direction_disp, trans_matrix_disp, type_disp_125,
rat_spline_curve_disp, rat_spline_sur_disp, offset_curve_disp, type_disp_132, type_disp_134,
type_disp_136, type_disp_138, offset_sur_disp, boundary_disp, curve_on_sur_disp,
bounded_sur_disp, trimmed_sur_disp, type_disp_146, type_disp_148, block_disp,
wedge_disp, cylinder_disp, cone_disp, sphere_disp, torus_disp,
solid_rev_disp, solid_ext_disp, ellipsoid_disp, bool_tree_disp, type_disp_182,
type_disp_184, solid_disp, plane_sur_disp, cyl_sur_disp, con_sur_disp,
sph_sur_disp, tor_sur_disp, type_disp_202, type_disp_204, type_disp_206,
type_disp_208, type_disp_210, type_disp_212, type_disp_213, type_disp_214,
type_disp_216, type_disp_218, type_disp_220, type_disp_222, type_disp_228,
type_disp_230, type_disp_302, type_disp_304, type_disp_306, type_disp_308,
type_disp_310, type_disp_312, type_disp_314, type_disp_316, type_disp_320,
type_disp_322, type_disp_402, type_disp_404, type_disp_406, subfig_instance_disp,
type_disp_410, type_disp_412, type_disp_414, type_disp_416, type_disp_418,
type_disp_420, type_disp_422, type_disp_430, vertex_list_disp, edge_list_disp,
loop_disp, face_disp, shell_disp
};
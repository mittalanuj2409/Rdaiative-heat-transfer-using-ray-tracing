//#include <bits/stdc++.h>
#include<iostream>
#include<string>
#include<vector>
#include<sstream>
#include<math.h>
#include"classes.h"
using namespace std;

#define PI 3.14159265
int Np_max = 100; //max number of actual points on any curve or one direction of surface + 1
int Np_min = 5;   //min number of actual points on any curve or one direction of surface + 1

int span_index(float u, int no_of_control_points, int degree, float* knots) {
	int n = no_of_control_points - 1;
	int p = degree;
	/*'u' should belong to (knots(i), knots(i+1)]*/
	int i;
	/*special case: u is equal to last p+1 knots. here we take (knots(n), knots(n+1)]*/
	if (u == knots[n + 1]) return n;
	/*binary search for i*/
	int low = p;
	int high = n + 1;
	i = (low + high) / 2;
	while (u < knots[i] || u >= knots[i + 1]) {
		if (u < knots[i]) high = i;
		else low = i;
		i = (low + high) / 2;
	}
	return i;
}

vector <float> b_spline_func(float u, int no_of_control_points, int degree, float* knots) {
	/*return B-spline functions N(i,degree) of required degree */
	vector <float> b_spline;
	int n = no_of_control_points - 1;
	int p = degree;
	/*'u' should belong to (knots(i), knots(i+1)]*/
	int i = span_index(u, n + 1, p, knots);

	/*Now,
		B-Spline functions N(i-p,p) to N(i,p) only are non-vanishing.
		Thus only these are calculated and else are equated to zero
	*/

	float* arr = (float*)malloc(sizeof(float) * 1e4);
	arr[0] = 1;
	for (int j = 1; j <= p; j++) arr[j] = 0;
	float down, m1, m2;
	for (int a = 0; a < p; a++) {
		down = arr[0];
		for (int b = 1; b <= p - a; b++) {
			m1 = knots[i + a] - knots[i - b] == 0 ? 1 : (u - knots[i - b]) / (knots[i + a] - knots[i - b]);
			m2 = knots[i + a + 1] - knots[i - b + 1] == 0 ? 1 : (knots[i + a + 1] - u) / (knots[i + a + 1] - knots[i - b + 1]);
			arr[b] = arr[b] * m1 + down * m2;
			down = arr[b];
		}
		arr[0] = arr[0] * (u - knots[i]) / (knots[i + a + 1] - knots[i]);
	}
	for (int j = 1; j <= i - p; j++) b_spline.push_back(0);
	for (int j = 0; j < p + 1; j++) b_spline.push_back(arr[j]);
	for (int j = 1; j <= n - i; j++) b_spline.push_back(0);
	free(arr);
	return b_spline;
}

int de_no(string l) {
	/*Extracting the Data Entry pointer number from a line.*/
	if (l.at(72) == 'P') return stoi(l.substr(64, 8));
	else return -1;
}

vector <string> parameters(string l, string de_no, char del1, char del2) {
	/*each parameter is seperated as string and returned in a
	  vector.*/
	vector <string> par;
	string data = "";
	int i = 0;
	while (l.at(i) != del2) {
		while (l.at(i) != del1 & l.at(i) != del2) {
			data.append(l.substr(i, 1));
			i++;
		}
		if (l.at(i) != del2) i++;
		par.push_back(data);
		data.clear();
	}
	par.push_back(de_no);
	/*de_no is used as pointer*/
	return par;
}

int entity[90] = { 0,
				100, 102, 104, 106, 108, //5
				110, 112, 114, 116, 118, //10
				120, 122, 123, 124, 125, //15
				126, 128, 130, 132, 134, //20
				136, 138, 140, 141, 142, //25
				143, 144, 146, 148, 150, //30
				152, 154, 156, 158, 160, //35
				162, 164, 168, 180, 182, //40
				184, 186, 190, 192, 194, //45
				196, 198, 202, 204, 206, //50
				208, 210, 212, 213, 214, //55
				216, 218, 220, 222, 228, //60
				230, 302, 304, 306, 308, //65
				310, 312, 314, 316, 320, //70
				322, 402, 404, 406, 408, //75
				410, 412, 414, 416, 418, //80
				420, 422, 430, 502, 504, //85
				508, 510, 514 };
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

float theta(float xc, float yc, float x, float y){
	float t = atan2((y-yc),(x-xc));
	t = t<0?t+2*PI:t;
	return t;
}

void type_0(vector <float> par_num) {}

void cir_arc_func(vector <float> par_num) {
	/*type 100. index 1*/
	cir_arc obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.z = par_num.at(1);
	obj.x_center = par_num.at(2);
	obj.y_center = par_num.at(3);
	obj.x_start = par_num.at(4);
	obj.y_start = par_num.at(5);
	obj.x_end = par_num.at(6);
	obj.y_end = par_num.at(7);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	float theta1, theta2, rad;
	rad = sqrt(pow((obj.y_start - obj.y_center),2) + pow((obj.x_start - obj.x_center),2));
	theta1 = theta(obj.x_center, obj.y_center, obj.x_start, obj.y_start);
	theta2 = theta(obj.x_center, obj.y_center, obj.x_end, obj.y_end);
	theta2 = theta2<theta1?theta2+2*PI:theta2;
	for(int i = 0; i<=Np_max; i++){
		obj.actual_points[i].x = obj.x_center + rad*cos(theta1 + (theta2 - theta1)/Np_max*i);
		obj.actual_points[i].y = obj.y_center + rad*sin(theta1 + (theta2 - theta1)/Np_max*i);
		obj.actual_points[i].z = obj.z;
	}

	cir_arc_vec.push_back(obj);
	
}

void comp_curve_func(vector <float> par_num) {
	/*type 102. index 2*/
	comp_curve obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.n = par_num.at(1);
	int i;
	for (i = 0; i < obj.n; i++) {
		obj.list[i] = par_num.at(1 + i + 1);
	}
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	comp_curve_vec.push_back(obj);
	
}

void conic_arc_func(vector <float> par_num) {
	/*type 104. index 3*/
	conic_arc obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.a = par_num.at(1);
	obj.b = par_num.at(2);
	obj.c = par_num.at(3);
	obj.d = par_num.at(4);
	obj.e = par_num.at(5);
	obj.f = par_num.at(6);
	obj.x_start = par_num.at(7);
	obj.y_start = par_num.at(8);
	obj.z_start = par_num.at(9);
	obj.x_end = par_num.at(10);
	obj.y_end = par_num.at(11);
	obj.z_end = par_num.at(12);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	float q1 = obj.a*(obj.c*obj.f - obj.e*obj.e/4) - (obj.b/2)*(obj.b*obj.f/2 - obj.d*obj.e/4) + (obj.d/2)*(obj.b*obj.e/4 - obj.c*obj.d/2);
	float q2 = obj.a*pbj.c - obj.b*obj.b/4;
	float q3 = obj.a + obj.c;
	if(q2>0 && q1*q3 <0){

	}
	else if(q2<0 && q1!=0){

	}
	else if(q2=0 && q1!=0){
	
		}
	conic_arc_vec.push_back(obj);
	
}

void type_106(vector <float> par_num) {}

void plane_func(vector <float> par_num) {
	/*type 108. index 5*/
	plane obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.a = par_num.at(1);
	obj.b = par_num.at(2);
	obj.c = par_num.at(3);
	obj.d = par_num.at(4);
	obj.x = par_num.at(5);
	obj.y = par_num.at(6);
	obj.z = par_num.at(7);
	obj.size = par_num.at(8);
	obj.pointer = par_num.at(9);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	plane_vec.push_back(obj);
	
}

void line_func(vector <float> par_num) {
	/*type 110. index 6*/
	line obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.x_start = par_num.at(1);
	obj.y_start = par_num.at(2);
	obj.z_start = par_num.at(3);
	obj.x_end = par_num.at(4);
	obj.y_end = par_num.at(5);
	obj.z_end = par_num.at(6);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	line_vec.push_back(obj);
	
}

void par_spline_curve_func(vector <float> par_num) {
	/*type 112. index 7*/
	par_spline_curve obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.type = par_num.at(1);
	obj.k = par_num.at(2);
	obj.dim = par_num.at(3);
	obj.n = par_num.at(4);
	for (int i = 0; i < obj.n; i++) {
		obj.t_points[i] = par_num.at(4 + i + 1);
		obj.list[i].ax = par_num.at(4 + obj.n + i * 4 + 1);
		obj.list[i].bx = par_num.at(4 + obj.n + i * 4 + 2);
		obj.list[i].cx = par_num.at(4 + obj.n + i * 4 + 3);
		obj.list[i].dx = par_num.at(4 + obj.n + i * 4 + 4);
		obj.list[i].ay = par_num.at(4 + obj.n * 5 + i * 4 + 1);
		obj.list[i].by = par_num.at(4 + obj.n * 5 + i * 4 + 2);
		obj.list[i].cy = par_num.at(4 + obj.n * 5 + i * 4 + 3);
		obj.list[i].dy = par_num.at(4 + obj.n * 5 + i * 4 + 4);
		obj.list[i].az = par_num.at(4 + obj.n * 9 + i * 4 + 1);
		obj.list[i].bz = par_num.at(4 + obj.n * 9 + i * 4 + 2);
		obj.list[i].cz = par_num.at(4 + obj.n * 9 + i * 4 + 3);
		obj.list[i].dz = par_num.at(4 + obj.n * 9 + i * 4 + 4);
	}
	obj.xn = par_num.at(4 + 13 * obj.n + 1);
	obj.xn1 = par_num.at(4 + 13 * obj.n + 2);
	obj.xn2 = par_num.at(4 + 13 * obj.n + 3);
	obj.xn3 = par_num.at(4 + 13 * obj.n + 4);
	obj.yn = par_num.at(4 + 13 * obj.n + 5);
	obj.yn1 = par_num.at(4 + 13 * obj.n + 6);
	obj.yn2 = par_num.at(4 + 13 * obj.n + 7);
	obj.yn3 = par_num.at(4 + 13 * obj.n + 8);
	obj.zn = par_num.at(4 + 13 * obj.n + 9);
	obj.zn1 = par_num.at(4 + 13 * obj.n + 10);
	obj.zn2 = par_num.at(4 + 13 * obj.n + 11);
	obj.zn3 = par_num.at(4 + 13 * obj.n + 12);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	par_spline_curve_vec.push_back(obj);
	
}

void par_spline_sur_func(vector <float> par_num) {
	/*type 114. index 8*/
	par_spline_sur obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.type = par_num.at(1);
	obj.k = par_num.at(2);
	obj.m = par_num.at(3);
	obj.n = par_num.at(4);
	for (int i = 0; i < obj.m; i++) {
		obj.tu_points[i] = par_num.at(4 + i + 1);
	}
	for (int j = 0; j < obj.n; j++) {
		obj.tv_points[j] = par_num.at(4 + j + obj.m + 1);
	}

	for (int l = 0; l < obj.m * obj.n; l++) {
		obj.p_ij[l] = par_num.at(4 + l + obj.m + obj.n + 1);
	}
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	par_spline_sur_vec.push_back(obj);
	
}

void point_func(vector <float> par_num) {
	/*type 116. index 9*/
	point_ obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.x = par_num.at(1);
	obj.y = par_num.at(2);
	obj.z = par_num.at(3);
	obj.pointer = int(par_num.at(4));
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	point_vec.push_back(obj);
	
}

void ruled_sur_func(vector <float> par_num) {
	/*type 118. index 10*/
	ruled_sur obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.pointer1 = par_num.at(1);
	obj.pointer2 = par_num.at(2);
	obj.dir_flag = par_num.at(3);
	obj.dev_flag = int(par_num.at(4));
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	ruled_sur_vec.push_back(obj);
	
}

void sur_revolution_func(vector <float> par_num) {
	/*type 120. index 11*/
	sur_revolution obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.pointer1 = par_num.at(1);
	obj.pointer2 = par_num.at(2);
	obj.sa = par_num.at(3);
	obj.ea = int(par_num.at(4));
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	sur_revolution_vec.push_back(obj);
	
}

void tab_cylinder_func(vector <float> par_num) {
	/*type 122. index 12*/
	tab_cylinder obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.pointer = par_num.at(1);
	obj.lx = par_num.at(2);
	obj.ly = par_num.at(3);
	obj.lz = par_num.at(4);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	tab_cylinder_vec.push_back(obj);
	
}

void direction_func(vector <float> par_num) {
	/*type 123. index 13*/
	direction obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.x = par_num.at(1);
	obj.y = par_num.at(2);
	obj.z = par_num.at(3);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	direction_vec.push_back(obj);
	
}

void trans_matrix_func(vector <float> par_num) {
	/*typr 124. index 14*/
	trans_matrix obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.r11 = par_num.at(1);
	obj.r12 = par_num.at(2);
	obj.r13 = par_num.at(3);
	obj.t1 = par_num.at(4);
	obj.r21 = par_num.at(5);
	obj.r22 = par_num.at(6);
	obj.r23 = par_num.at(7);
	obj.t2 = par_num.at(8);
	obj.r31 = par_num.at(9);
	obj.r32 = par_num.at(10);
	obj.r33 = par_num.at(11);
	obj.t3 = par_num.at(12);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	trans_matrix_vec.push_back(obj);
	
}

void type_125(vector <float> par_num) {}

void rat_spline_curve_func(vector <float> par_num) {
	/*type 126. index 16*/
	rat_spline_curve obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.k = par_num.at(1);
	obj.m = par_num.at(2);
	obj.f1 = par_num.at(3);
	obj.f2 = par_num.at(4);
	obj.f3 = par_num.at(5);
	obj.f4 = par_num.at(6);
	for (int a = 0; a < 2 + obj.k + obj.m; a++) {
		obj.knot[a] = par_num.at(6 + a + 1);
	}
	/*knot vector should be normalized
	  It is non-decreasing order
	  Thus, each term is divided by the last term
	*/
	for (int a = 0; a < 2 + obj.k + obj.m; a++) {
		obj.knot[a] = obj.knot[a] / obj.knot[1 + obj.k + obj.m];
	}

	float max_x, max_y, max_z, min_x, min_y, min_z;
	max_x = par_num.at(6 + (2 + obj.k + obj.m) + (1 + obj.k) + 1);
	max_y = par_num.at(6 + (2 + obj.k + obj.m) + (1 + obj.k) + 2);
	max_z = par_num.at(6 + (2 + obj.k + obj.m) + (1 + obj.k) + 3);
	min_x = max_x;
	min_y = max_y;
	min_z = max_z;

	for (int b = 0; b < 1 + obj.k; b++) {
		obj.weight[b] = par_num.at(6 + (2 + obj.k + obj.m) + b + 1);
		obj.control_points[b].x = par_num.at(6 + (2 + obj.k + obj.m) + (1 + obj.k) + 3 * b + 1);
		obj.control_points[b].z = par_num.at(6 + (2 + obj.k + obj.m) + (1 + obj.k) + 3 * b + 3);
		obj.control_points[b].y = par_num.at(6 + (2 + obj.k + obj.m) + (1 + obj.k) + 3 * b + 2);
		if (max_x < obj.control_points[b].x) max_x = obj.control_points[b].x;
		if (min_x > obj.control_points[b].x) min_x = obj.control_points[b].x;
		if (max_y < obj.control_points[b].y) max_y = obj.control_points[b].y;
		if (min_y > obj.control_points[b].y) min_y = obj.control_points[b].y;
		if (max_z < obj.control_points[b].z) max_z = obj.control_points[b].z;
		if (min_z > obj.control_points[b].z) min_z = obj.control_points[b].z;
	}
	obj.v0 = par_num.at(6 + 5 * (1 + obj.k) + 1 + obj.m + 1);
	obj.v1 = par_num.at(6 + 5 * (1 + obj.k) + 1 + obj.m + 2);
	obj.xn = par_num.at(6 + 5 * (1 + obj.k) + 1 + obj.m + 3);
	obj.yn = par_num.at(6 + 5 * (1 + obj.k) + 1 + obj.m + 4);
	obj.zn = par_num.at(6 + 5 * (1 + obj.k) + 1 + obj.m + 5);

	int Np_ = (max_x - min_x) >= (max_y - min_y) ?
		((max_x - min_x) >= (max_z - min_z) ? (max_x - min_x) : (max_z - min_z)) :
		((max_y - min_y) >= (max_z - min_z) ? (max_y - min_y) : (max_z - min_z));
	Np_ = Np_ > Np_max ? Np_max : (Np_ < Np_min ? Np_min : Np_);
	vector <float> b_spline;
	int s, i;
	for (int j = 0; j < Np_; j++) {
		b_spline = b_spline_func(float(j) / float(Np_), 1 + obj.k, obj.m, obj.knot);
		s = span_index(float(j) / float(Np_), 1 + obj.k, obj.m, obj.knot);
		float x_i = 0, y_i = 0, z_i = 0, w_i = 0;
		for (int l = 0; l < 1 + obj.m; l++) {
			i = s - obj.m + l;
			x_i += b_spline.at(i) * obj.control_points[i].x * obj.weight[i];
			y_i += b_spline.at(i) * obj.control_points[i].y * obj.weight[i];
			z_i += b_spline.at(i) * obj.control_points[i].z * obj.weight[i];
			w_i += b_spline.at(i) * obj.weight[i];
		}
		obj.actual_points[j].x = x_i / w_i;
		obj.actual_points[j].y = y_i / w_i;
		obj.actual_points[j].z = z_i / w_i;
	}
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	rat_spline_curve_vec.push_back(obj);
	
}

void rat_spline_sur_func(vector <float> par_num) {
	/*type 128. index 17*/
	rat_spline_sur obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.k1 = par_num.at(1);
	obj.k2 = par_num.at(2);
	obj.m1 = par_num.at(3);
	obj.m2 = par_num.at(4);
	obj.f1 = par_num.at(5);
	obj.f2 = par_num.at(6);
	obj.f3 = par_num.at(7);
	obj.f4 = par_num.at(8);
	obj.f5 = par_num.at(9);
	for (int a = 0; a < 2 + obj.k1 + obj.m1; a++) {
		obj.knot1[a] = par_num.at(9 + a + 1);
	}
	for (int a = 0; a < 2 + obj.k2 + obj.m2; a++) {
		obj.knot2[a] = par_num.at(9 + (2 + obj.k1 + obj.m1) + a + 1);
	}
	/*knot vector should be normalized
	  It is non-decreasing order
	  Thus, each term is divided by the last term
	*/
	for (int a = 0; a < 2 + obj.k1 + obj.m1; a++) {
		obj.knot1[a] = obj.knot1[a] / obj.knot1[1 + obj.k1 + obj.m1];
	}
	for (int a = 0; a < 2 + obj.k2 + obj.m2; a++) {
		obj.knot2[a] = obj.knot2[a] / obj.knot2[1 + obj.k2 + obj.m2];
	}
	float max_x, max_y, max_z, min_x, min_y, min_z, Np;
	int Np_1 = 0; //Number of points in direction 1
	int Np_2 = 0; //Number of points in direction 2
	max_x = par_num.at(9 + (4 + obj.k1 + obj.k2 + obj.m1 + obj.m2) + (1 + obj.k1) * (1 + obj.k2) + 1);
	max_y = par_num.at(9 + (4 + obj.k1 + obj.k2 + obj.m1 + obj.m2) + (1 + obj.k1) * (1 + obj.k2) + 2);
	max_z = par_num.at(9 + (4 + obj.k1 + obj.k2 + obj.m1 + obj.m2) + (1 + obj.k1) * (1 + obj.k2) + 3);
	min_x = max_x;
	min_y = max_y;
	min_z = max_z;
	for (int b = 0; b < (1 + obj.k2); b++) {
		for (int c = 0; c < (1 + obj.k1); c++) {
			obj.weight[(1 + obj.k1) * b + c] = par_num.at(9 + (4 + obj.k1 + obj.k2 + obj.m1 + obj.m2) + (1 + obj.k1) * b + c + 1);
			obj.control_points[(1 + obj.k1) * b + c].x = par_num.at(9 + (4 + obj.k1 + obj.k2 + obj.m1 + obj.m2) + (1 + obj.k1) * (1 + obj.k2) + 3 * ((1 + obj.k1) * b + c) + 1);
			obj.control_points[(1 + obj.k1) * b + c].y = par_num.at(9 + (4 + obj.k1 + obj.k2 + obj.m1 + obj.m2) + (1 + obj.k1) * (1 + obj.k2) + 3 * ((1 + obj.k1) * b + c) + 2);
			obj.control_points[(1 + obj.k1) * b + c].z = par_num.at(9 + (4 + obj.k1 + obj.k2 + obj.m1 + obj.m2) + (1 + obj.k1) * (1 + obj.k2) + 3 * ((1 + obj.k1) * b + c) + 3);
			if (max_x < obj.control_points[(1 + obj.k1) * b + c].x) max_x = obj.control_points[(1 + obj.k1) * b + c].x;
			if (min_x > obj.control_points[(1 + obj.k1) * b + c].x) min_x = obj.control_points[(1 + obj.k1) * b + c].x;
			if (max_y < obj.control_points[(1 + obj.k1) * b + c].y) max_y = obj.control_points[(1 + obj.k1) * b + c].y;
			if (min_y > obj.control_points[(1 + obj.k1) * b + c].y) min_y = obj.control_points[(1 + obj.k1) * b + c].y;
			if (max_z < obj.control_points[(1 + obj.k1) * b + c].z) max_z = obj.control_points[(1 + obj.k1) * b + c].z;
			if (min_z > obj.control_points[(1 + obj.k1) * b + c].z) min_z = obj.control_points[(1 + obj.k1) * b + c].z;
		}
		Np = (max_x - min_x) >= (max_y - min_y) ?
			((max_x - min_x) >= (max_z - min_z) ? (max_x - min_x) : (max_z - min_z)) :
			((max_y - min_y) >= (max_z - min_z) ? (max_y - min_y) : (max_z - min_z));
		Np_1 = Np > Np_1 ? Np : Np_1;
	}
	for (int b = 0; b < (1 + obj.k1); b++) {
		for (int c = 0; c < (1 + obj.k2); c++) {
			if (max_x < obj.control_points[(1 + obj.k2) * b + c].x) max_x = obj.control_points[(1 + obj.k2) * b + c].x;
			if (min_x > obj.control_points[(1 + obj.k2) * b + c].x) min_x = obj.control_points[(1 + obj.k2) * b + c].x;
			if (max_y < obj.control_points[(1 + obj.k2) * b + c].y) max_y = obj.control_points[(1 + obj.k2) * b + c].y;
			if (min_y > obj.control_points[(1 + obj.k2) * b + c].y) min_y = obj.control_points[(1 + obj.k2) * b + c].y;
			if (max_z < obj.control_points[(1 + obj.k2) * b + c].z) max_z = obj.control_points[(1 + obj.k2) * b + c].z;
			if (min_z > obj.control_points[(1 + obj.k2) * b + c].z) min_z = obj.control_points[(1 + obj.k2) * b + c].z;
		}
		Np = (max_x - min_x) >= (max_y - min_y) ?
			((max_x - min_x) >= (max_z - min_z) ? (max_x - min_x) : (max_z - min_z)) :
			((max_y - min_y) >= (max_z - min_z) ? (max_y - min_y) : (max_z - min_z));
		Np_2 = Np > Np_2 ? Np : Np_2;
	}
	Np_1 = Np_1 > Np_max ? Np_max : (Np_1 < Np_min ? Np_min : Np_1);
	Np_2 = Np_2 > Np_max ? Np_max : (Np_2 < Np_min ? Np_min : Np_2);
	vector <float> b_spline_1, b_spline_2;
	int s1, s2;
	for (int j = 0; j < Np_1; j++) {
		for (int m = 0; m < Np_2; m++) {
			b_spline_1 = b_spline_func(float(j) / float(Np_1), 1 + obj.k1, obj.m1, obj.knot1);
			s1 = span_index(float(j) / float(Np_1), 1 + obj.k1, obj.m1, obj.knot1);
			b_spline_2 = b_spline_func(float(m) / float(Np_2), 1 + obj.k2, obj.m2, obj.knot2);
			s2 = span_index(float(m) / float(Np_2), 1 + obj.k2, obj.m2, obj.knot2);
			float xt = 0, yt = 0, zt = 0, wt = 0;
			for (int l = 0; l < 1 + obj.m1; l++) {
				float x = 0, y = 0, z = 0, w = 0;
				for (int n = 0; n < 1 + obj.m2; n++) {
					int i = s1 - obj.m1 + l + (n + s2 - obj.m2) * (1 + obj.k1);
					x += obj.control_points[i].x * b_spline_2.at(n + s2 - obj.m2);
					y += obj.control_points[i].y * b_spline_2.at(n + s2 - obj.m2);
					z += obj.control_points[i].z * b_spline_2.at(n + s2 - obj.m2);
					w += obj.weight[i] * b_spline_2.at(n + s2 - obj.m2);
				}
				xt += x * b_spline_1.at(s1 - obj.m1 + l);
				yt += y * b_spline_1.at(s1 - obj.m1 + l);
				zt += z * b_spline_1.at(s1 - obj.m1 + l);
				wt += w * b_spline_1.at(s1 - obj.m1 + l);
			}
			obj.actual_points[Np_1 * (m)].x = xt / wt;
			obj.actual_points[Np_1 * (m)].y = yt / wt;
			obj.actual_points[Np_1 * (m)].z = zt / wt;
		}
	}
	obj.u0 = par_num.at(9 + (4 + obj.k1 + obj.k2 + obj.m1 + obj.m2) + 4 * (1 + obj.k1) * (1 + obj.k2) + 1);
	obj.u1 = par_num.at(9 + (4 + obj.k1 + obj.k2 + obj.m1 + obj.m2) + 4 * (1 + obj.k1) * (1 + obj.k2) + 2);
	obj.v0 = par_num.at(9 + (4 + obj.k1 + obj.k2 + obj.m1 + obj.m2) + 4 * (1 + obj.k1) * (1 + obj.k2) + 3);
	obj.v1 = par_num.at(9 + (4 + obj.k1 + obj.k2 + obj.m1 + obj.m2) + 4 * (1 + obj.k1) * (1 + obj.k2) + 4);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	rat_spline_sur_vec.push_back(obj);
	
}

void offset_curve_func(vector <float> par_num) {
	/*type 130. index 18*/
	offset_curve obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.pointer1 = par_num.at(1);
	obj.f1 = par_num.at(2);
	obj.pointer2 = par_num.at(3);
	obj.dim = par_num.at(4);
	obj.f2 = par_num.at(5);
	obj.d1 = par_num.at(6);
	obj.td1 = par_num.at(7);
	obj.d2 = par_num.at(8);
	obj.td2 = par_num.at(9);
	obj.x = par_num.at(10);
	obj.y = par_num.at(11);
	obj.z = par_num.at(12);
	obj.tt1 = par_num.at(13);
	obj.tt2 = par_num.at(14);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	offset_curve_vec.push_back(obj);
	
}

void type_132(vector <float> par_num) {}

void type_134(vector <float> par_num) {}

void type_136(vector <float> par_num) {}

void type_138(vector <float> par_num) {}

void offset_sur_func(vector <float> par_num) {
	/*type 140. index 23*/
	offset_sur obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.nx = par_num.at(1);
	obj.ny = par_num.at(2);
	obj.nz = par_num.at(3);
	obj.d = par_num.at(4);
	obj.pointer = par_num.at(5);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	offset_sur_vec.push_back(obj);
	
}

void boundary_func(vector <float> par_num) {
	/*type 141. index 24*/
	int j, k = 0;
	boundary obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.type = par_num.at(1);
	obj.pref = par_num.at(2);
	obj.pointer1 = par_num.at(3);
	obj.n = par_num.at(4);
	//k = obj.curves[0].k;
	for (int i = 0; i < obj.n; i++) {
		obj.curves[i].pointer = par_num.at(4 + (k)+1);
		obj.curves[i].f = par_num.at(4 + (k)+2);
		obj.curves[i].k = par_num.at(4 + (k)+3);
		for (j = 0; j < k; j++) {
			obj.curves[i].p[j] = par_num.at(4 + (k)+j + 4);
		}
		k += 3 + obj.curves[i].k;
	}
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	boundary_vec.push_back(obj);
	
}

void curve_on_sur_func(vector <float> par_num) {
	/*type 142. index 25*/
	curve_on_sur obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.flag = par_num.at(1);
	obj.pointer1 = par_num.at(2);
	obj.pointer2 = par_num.at(3);
	obj.pointer3 = par_num.at(4);
	obj.rep = par_num.at(5);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	curve_on_sur_vec.push_back(obj);
	
}

void bounded_sur_func(vector <float> par_num) {
	/*type 143. index 26*/
	bounded_sur obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.type = par_num.at(1);
	obj.pointer1 = par_num.at(2);
	obj.n = par_num.at(3);
	for (int i = 0; i < obj.n; i++) {
		obj.list[i] = par_num.at(3 + i + 1);
	}
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	bounded_sur_vec.push_back(obj);
	
}

void trimmed_sur_func(vector <float> par_num) {
	/*type 144. index 27*/
	trimmed_sur obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.pointer1 = par_num.at(1);
	obj.flag = par_num.at(2);
	obj.n = par_num.at(3);
	obj.pointer2 = par_num.at(4);
	for (int i = 0; i < obj.n; i++) {
		obj.list[i] = par_num.at(4 + i + 1);
	}
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	trimmed_sur_vec.push_back(obj);
	
}

void type_146(vector <float> par_num) {}

void type_148(vector <float> par_num) {}

void block_func(vector <float> par_num) {
	/*typr 150. index 30*/
	block obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.length_x = par_num.at(1);
	obj.length_y = par_num.at(2);
	obj.length_z = par_num.at(3);
	obj.x_corner = par_num.at(4);
	obj.y_corner = par_num.at(5);
	obj.z_corner = par_num.at(6);
	obj.x_i = par_num.at(7);
	obj.x_j = par_num.at(8);
	obj.x_k = par_num.at(9);
	obj.z_i = par_num.at(10);
	obj.z_j = par_num.at(11);
	obj.z_k = par_num.at(12);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	block_vec.push_back(obj);
	
}

void wedge_func(vector <float> par_num) {
	/*typr 152. index 31*/
	wedge obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.length_x = par_num.at(1);
	obj.length_y = par_num.at(2);
	obj.length_z = par_num.at(3);
	obj.tlx = par_num.at(4);
	obj.x_corner = par_num.at(5);
	obj.y_corner = par_num.at(6);
	obj.z_corner = par_num.at(7);
	obj.x_i = par_num.at(8);
	obj.x_j = par_num.at(9);
	obj.x_k = par_num.at(10);
	obj.z_i = par_num.at(11);
	obj.z_j = par_num.at(12);
	obj.z_k = par_num.at(13);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	wedge_vec.push_back(obj);
	
}

void cylinder_func(vector <float> par_num) {
	/*typr 154. index 32*/
	cylinder obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.h = par_num.at(1);
	obj.r = par_num.at(2);
	obj.x_center = par_num.at(3);
	obj.y_center = par_num.at(4);
	obj.z_center = par_num.at(5);
	obj.i = par_num.at(6);
	obj.j = par_num.at(7);
	obj.k = par_num.at(8);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	cylinder_vec.push_back(obj);
	
}

void cone_func(vector <float> par_num) {
	/*typr 156. index 33*/
	cone obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.h = par_num.at(1);
	obj.r_large = par_num.at(2);
	obj.r_small = par_num.at(3);
	obj.x_center = par_num.at(4);
	obj.y_center = par_num.at(5);
	obj.z_center = par_num.at(6);
	obj.i = par_num.at(7);
	obj.j = par_num.at(8);
	obj.k = par_num.at(9);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	cone_vec.push_back(obj);
	
}

void sphere_func(vector <float> par_num) {
	/*type 158. index 34*/
	sphere obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.r = par_num.at(1);
	obj.x = par_num.at(2);
	obj.y = par_num.at(3);
	obj.z = par_num.at(4);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	sphere_vec.push_back(obj);
	
}

void torus_func(vector <float> par_num) {
	/*typr 160. index 35*/
	torus obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.r1 = par_num.at(1);
	obj.r2 = par_num.at(2);
	obj.x_center = par_num.at(3);
	obj.y_center = par_num.at(4);
	obj.z_center = par_num.at(5);
	obj.i = par_num.at(6);
	obj.j = par_num.at(7);
	obj.k = par_num.at(8);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	torus_vec.push_back(obj);
	
}

void solid_rev_func(vector <float> par_num) {
	/*type 162. index 36*/
	solid_rev obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.pointer = par_num.at(1);
	obj.f = par_num.at(2);
	obj.x = par_num.at(3);
	obj.y = par_num.at(4);
	obj.z = par_num.at(5);
	obj.i = par_num.at(6);
	obj.j = par_num.at(7);
	obj.k = par_num.at(8);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	solid_rev_vec.push_back(obj);
	
}

void solid_ext_func(vector <float> par_num) {
	/*type 164. index 37*/
	solid_ext obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.pointer = par_num.at(1);
	obj.l = par_num.at(2);
	obj.i = par_num.at(3);
	obj.j = par_num.at(4);
	obj.k = par_num.at(5);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	solid_ext_vec.push_back(obj);
	
}

void ellipsoid_func(vector <float> par_num) {
	/*type 168. index 38*/
	ellipsoid obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.length_x = par_num.at(1);
	obj.length_y = par_num.at(2);
	obj.length_z = par_num.at(3);
	obj.x_center = par_num.at(4);
	obj.y_center = par_num.at(5);
	obj.z_center = par_num.at(6);
	obj.x_i = par_num.at(7);
	obj.x_j = par_num.at(8);
	obj.x_k = par_num.at(9);
	obj.z_i = par_num.at(10);
	obj.z_j = par_num.at(11);
	obj.z_k = par_num.at(12);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	ellipsoid_vec.push_back(obj);
	
}

void bool_tree_func(vector <float> par_num) {
	/*type 180. index 39*/
	bool_tree obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.n = par_num.at(1);
	for (int i = 0; i < obj.n; i++) {
		obj.pointer[i] = par_num.at(1 + i + 1);
	}
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	bool_tree_vec.push_back(obj);
}

void type_182(vector <float> par_num) {}

void type_184(vector <float> par_num) {}

void solid_func(vector <float> par_num) {
	/*type 186. index 42*/
	solid obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.pointer = par_num.at(1);
	obj.flag = par_num.at(2);
	obj.n = par_num.at(3);
	int j;
	for (j = 0; j < obj.n; j++) {
		obj.list[j].pointer = par_num.at(3 + 2 * j + 1);
		obj.list[j].vflag = par_num.at(3 + 2 * j + 2);
	}
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	solid_vec.push_back(obj);
	
}

void plane_sur_func(vector <float> par_num) {
	/*type 190. index 43*/
	plane_sur obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.pointer1 = par_num.at(1);
	obj.pointer2 = par_num.at(2);
	obj.pointer3 = par_num.at(3);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	plane_sur_vec.push_back(obj);
	
}


void cyl_sur_func(vector <float> par_num) {
	/*type 192. index 44*/
	cyl_sur obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.point_ = par_num.at(1);
	obj.axis = par_num.at(2);
	obj.ref = par_num.at(3);
	obj.r = par_num.at(4);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	cyl_sur_vec.push_back(obj);
	
}

void con_sur_func(vector <float> par_num) {
	/*type 194. index 45*/
	con_sur obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.point_ = par_num.at(1);
	obj.axis = par_num.at(2);
	obj.ref = par_num.at(3);
	obj.r = par_num.at(4);
	obj.ang = par_num.at(5);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	con_sur_vec.push_back(obj);
	
}

void sph_sur_func(vector <float> par_num) {
	/*type 196. index 46*/
	sph_sur obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.point_ = par_num.at(1);
	obj.axis = par_num.at(2);
	obj.ref = par_num.at(3);
	obj.r = par_num.at(4);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	sph_sur_vec.push_back(obj);
	
}

void tor_sur_func(vector <float> par_num) {
	/*type 198. index 47*/
	tor_sur obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.point_ = par_num.at(1);
	obj.axis = par_num.at(2);
	obj.ref = par_num.at(3);
	obj.r1 = par_num.at(4);
	obj.r2 = par_num.at(5);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	tor_sur_vec.push_back(obj);
	
}

void type_202(vector <float> par_num) {}

void type_204(vector <float> par_num) {}

void type_206(vector <float> par_num) {}

void type_208(vector <float> par_num) {}

void type_210(vector <float> par_num) {}

void type_212(vector <float> par_num) {}

void type_213(vector <float> par_num) {}

void type_214(vector <float> par_num) {}

void type_216(vector <float> par_num) {}

void type_218(vector <float> par_num) {}

void type_220(vector <float> par_num) {}

void type_222(vector <float> par_num) {}

void type_228(vector <float> par_num) {}

void type_230(vector <float> par_num) {}

void type_302(vector <float> par_num) {}

void type_304(vector <float> par_num) {}

void type_306(vector <float> par_num) {}

void type_308(vector <float> par_num) {} //dummy

void subfigure_func(vector <string> par) {
	/*type 308. index 64*/
	subfigure obj;
	obj.de_no = stoi(*(par.rbegin()));
	obj.depth = stoi(par.at(1));
	obj.name = string(par.at(2));
	obj.n = stoi(par.at(3));
	for (int i = 0; i < obj.n; i++) {
		obj.pointer[i] = stoi(par.at(3 + i + 1));
	}
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	subfigure_vec.push_back(obj);
}

void type_310(vector <float> par_num) {}

void type_312(vector <float> par_num) {}

void type_314(vector <float> par_num) {} //dummy

void color_func(vector <string> par) {
	/*type 314. index 68*/
	color obj;
	stringstream a(par.at(1)), b(par.at(2)), c(par.at(3));
	obj.de_no = stoi(*(par.rbegin()));
	a >> obj.red;
	b >> obj.blue;
	c >> obj.green;
	obj.name = string(par.at(4));
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	color_vec.push_back(obj);
}

void type_316(vector <float> par_num) {}

void type_320(vector <float> par_num) {}

void type_322(vector <float> par_num) {}

void type_402(vector <float> par_num) {}

void type_404(vector <float> par_num) {}

void type_406(vector <float> par_num) {}

void subfig_instance_func(vector <float> par_num) {
	/*type 408. index 75*/
	subfig_instance obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.sd = par_num.at(1);
	obj.x = par_num.at(2);
	obj.y = par_num.at(3);
	obj.z = par_num.at(4);
	obj.s = par_num.at(5);
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	subfig_instance_vec.push_back(obj);
}

void type_410(vector <float> par_num) {}

void type_412(vector <float> par_num) {}

void type_414(vector <float> par_num) {}

void type_416(vector <float> par_num) {}

void type_418(vector <float> par_num) {}

void type_420(vector <float> par_num) {}

void type_422(vector <float> par_num) {}

void type_430(vector <float> par_num) {}

void vertex_list_func(vector <float> par_num) {
	/*type 502. index 84*/
	vertex_list obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.n = par_num.at(1);
	int j;
	for (j = 0; j < obj.n; j++) {
		obj.list[j].x = par_num.at(1 + 3 * j + 1);
		obj.list[j].y = par_num.at(1 + 3 * j + 2);
		obj.list[j].z = par_num.at(1 + 3 * j + 3);
	}
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	vertex_list_vec.push_back(obj);
}

void edge_list_func(vector <float> par_num) {
	/*type 504. index 85*/
	edge_list obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.n = par_num.at(1);
	int j;
	for (j = 0; j < obj.n; j++) {
		obj.list[j].pointer1 = par_num.at(1 + 5 * j + 1);
		obj.list[j].start_pointer = par_num.at(1 + 5 * j + 2);
		obj.list[j].start = par_num.at(1 + 5 * j + 3);
		obj.list[j].end_pointer = par_num.at(1 + 5 * j + 4);
		obj.list[j].end = par_num.at(1 + 5 * j + 5);
	}
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	edge_list_vec.push_back(obj);
}

void loop_func(vector <float> par_num) {
	/*type 508. index 86*/
	loop obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.n = par_num.at(1);
	int j;
	for (j = 0; j < obj.n; j++) {
		obj.list2[j].type = par_num.at(1 + 5 * j + 1);
		obj.list2[j].pointer = par_num.at(1 + 5 * j + 2);
		obj.list2[j].index = par_num.at(1 + 5 * j + 3);
		obj.list2[j].flag = par_num.at(1 + 5 * j + 4);
		obj.list2[j].k = par_num.at(1 + 5 * j + 5);
		//if(obj.list2[j].k == 0) continue;
		for (int i = 0; i < obj.list2[j].k; i++) {
			obj.list2[j].list1[i].iso_flag = par_num.at(1 + 5 * j + 5 + 2 * i + 1);
			obj.list2[j].list1[i].pointer = par_num.at(1 + 5 * j + 5 + 2 * i + 2);
		}
	}
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	loop_vec.push_back(obj);
}

void face_func(vector <float> par_num) {
	/*type 510. index 87*/
	face obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.pointer = par_num.at(1);
	obj.n = par_num.at(2);
	obj.flag = par_num.at(3);
	int i;
	for (i = 0; i < obj.n; i++) {
		obj.loop_pointer[i] = par_num.at(3 + i + 1);
	}
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	face_vec.push_back(obj);
}

void shell_func(vector <float> par_num) {
	/*type 514. index 88*/
	shell_ obj;
	obj.de_no = int(*(par_num.rbegin()));
	obj.n = par_num.at(1);
	int j;
	for (j = 0; j < obj.n; j++) {
		obj.list[j].pointer = par_num.at(1 + 2 * j + 1);
		obj.list[j].flag = par_num.at(1 + 2 * j + 2);
	}
	obj.t_matrix = t_matrix_pointer.at((obj.de_no - 1) / 2);
	obj.form = form_number.at((obj.de_no - 1) / 2);
	shell_vec.push_back(obj);
}


void (*function_list[])(vector <float>) =
{ type_0,
cir_arc_func, comp_curve_func, conic_arc_func, type_106, plane_func,
line_func, par_spline_curve_func, par_spline_sur_func, point_func, ruled_sur_func,
sur_revolution_func, tab_cylinder_func, direction_func, trans_matrix_func, type_125,
rat_spline_curve_func, rat_spline_sur_func, offset_curve_func, type_132, type_134,
type_136, type_138, offset_sur_func, boundary_func, curve_on_sur_func,
bounded_sur_func, trimmed_sur_func, type_146, type_148, block_func,
wedge_func, cylinder_func, cone_func, sphere_func, torus_func,
solid_rev_func, solid_ext_func, ellipsoid_func, bool_tree_func, type_182,
type_184, solid_func, plane_sur_func, cyl_sur_func, con_sur_func,
sph_sur_func, tor_sur_func, type_202, type_204, type_206,
type_208, type_210, type_212, type_213, type_214,
type_216, type_218, type_220, type_222, type_228,
type_230, type_302, type_304, type_306, type_308,
type_310, type_312, type_314, type_316, type_320,
type_322, type_402, type_404, type_406, subfig_instance_func,
type_410, type_412, type_414, type_416, type_418,
type_420, type_422, type_430, vertex_list_func, edge_list_func,
loop_func, face_func, shell_func
};



void contents() {
	system("CLS");
	cout << "type 100" << "\t" << "circular arc __________" << "\t" << cir_arc_vec.size() << endl;
	cout << "type 102" << "\t" << "composite arc _________" << "\t" << comp_curve_vec.size() << endl;
	cout << "type 104" << "\t" << "conic arc _____________" << "\t" << conic_arc_vec.size() << endl;
	cout << "type 108" << "\t" << "plane _________________" << "\t" << plane_vec.size() << endl;
	cout << "type 110" << "\t" << "line __________________" << "\t" << line_vec.size() << endl;
	cout << "type 112" << "\t" << "parametric spline curve" << "\t" << par_spline_curve_vec.size() << endl;
	cout << "type 116" << "\t" << "point _________________" << "\t" << point_vec.size() << endl;
	cout << "type 118" << "\t" << "ruled surface _________" << "\t" << ruled_sur_vec.size() << endl;
	cout << "type 120" << "\t" << "surface of revolution  " << "\t" << sur_revolution_vec.size() << endl;
	cout << "type 122" << "\t" << "tabulated cylinder ____" << "\t" << tab_cylinder_vec.size() << endl;
	cout << "type 123" << "\t" << "direction _____________" << "\t" << direction_vec.size() << endl;
	cout << "type 124" << "\t" << "transformation matrix  " << "\t" << trans_matrix_vec.size() << endl;
	cout << "type 126" << "\t" << "rational spline curve  " << "\t" << rat_spline_curve_vec.size() << endl;
	cout << "type 128" << "\t" << "rational spline surface" << "\t" << rat_spline_sur_vec.size() << endl;
	cout << "type 130" << "\t" << "offset curve __________" << "\t" << offset_curve_vec.size() << endl;
	cout << "type 140" << "\t" << "offset surface ________" << "\t" << offset_sur_vec.size() << endl;
	cout << "type 141" << "\t" << "boundary ______________" << "\t" << boundary_vec.size() << endl;
	cout << "type 142" << "\t" << "curve on surface ______" << "\t" << curve_on_sur_vec.size() << endl;
	cout << "type 143" << "\t" << "bounded surface _______" << "\t" << bounded_sur_vec.size() << endl;
	cout << "type 144" << "\t" << "trimmed surface _______" << "\t" << trimmed_sur_vec.size() << endl;
	cout << "type 150" << "\t" << "block _________________" << "\t" << block_vec.size() << endl;
	cout << "type 152" << "\t" << "wedge _________________" << "\t" << wedge_vec.size() << endl;
	cout << "type 154" << "\t" << "cylinder ______________" << "\t" << cylinder_vec.size() << endl;
	cout << "type 156" << "\t" << "cone __________________" << "\t" << cone_vec.size() << endl;
	cout << "type 158" << "\t" << "sphere ________________" << "\t" << sphere_vec.size() << endl;
	cout << "type 160" << "\t" << "torus _________________" << "\t" << torus_vec.size() << endl;
	cout << "type 162" << "\t" << "solid of revolution ___" << "\t" << solid_rev_vec.size() << endl;
	cout << "type 164" << "\t" << "solid extrusion _______" << "\t" << solid_ext_vec.size() << endl;
	cout << "type 168" << "\t" << "ellipsoid _____________" << "\t" << ellipsoid_vec.size() << endl;
	cout << "type 180" << "\t" << "boolean tree __________" << "\t" << bool_tree_vec.size() << endl;
	cout << "type 186" << "\t" << "solid _________________" << "\t" << solid_vec.size() << endl;
	cout << "type 190" << "\t" << "plane surface _________" << "\t" << plane_sur_vec.size() << endl;
	cout << "type 192" << "\t" << "cylindrical surface ___" << "\t" << cyl_sur_vec.size() << endl;
	cout << "type 194" << "\t" << "conical surface _______" << "\t" << con_sur_vec.size() << endl;
	cout << "type 196" << "\t" << "spherical surface _____" << "\t" << sph_sur_vec.size() << endl;
	cout << "type 198" << "\t" << "toroidal surface ______" << "\t" << tor_sur_vec.size() << endl;
	cout << "type 308" << "\t" << "subfigure _____________" << "\t" << subfigure_vec.size() << endl;
	cout << "type 314" << "\t" << "color _________________" << "\t" << color_vec.size() << endl;
	cout << "type 408" << "\t" << "subfigure instance ____" << "\t" << subfig_instance_vec.size() << endl;
	cout << "type 502" << "\t" << "vertex list ___________" << "\t" << vertex_list_vec.size() << endl;
	cout << "type 504" << "\t" << "edge list _____________" << "\t" << edge_list_vec.size() << endl;
	cout << "type 508" << "\t" << "loop __________________" << "\t" << loop_vec.size() << endl;
	cout << "type 510" << "\t" << "face __________________" << "\t" << face_vec.size() << endl;
	cout << "type 514" << "\t" << "shell _________________" << "\t" << shell_vec.size() << endl;
}

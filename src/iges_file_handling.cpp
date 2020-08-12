//#include<bits/stdc++.h>

#include<iostream>
#include<string>
#include<fstream>
#include<vector>
#include<sstream>
#include <buid_entity.h>
#include"sur_funcs.cpp"
using namespace std;

/*https://wiki.eclipse.org/IGES_file_Specification
Refer this documentation for iges files*/

/*Information about the delimiters is obtained from first two
  words of Global section*/
char del1 = ',';
char del2 = ';';

int main() {
	string filename, exit_ = "";
	cout << "enter filename with extension" << endl;
	cin >> filename;
	ifstream f1;
	f1.open(filename);
	if (!f1) {
		cerr << "unable to open " << filename << endl << "Please restart" << endl << "type \"exit\" to exit" << endl;;

		while (exit_ != "exit" && exit_ != "EXIT" && exit_ != "Exit") cin >> exit_;
		exit(1);
	}

	string line1, line2;
	getline(f1, line1);
	while (line1.at(72) == 'S' || line1.at(72) == 'G') {
		/*All the required information about our model can be obtained
		  from Parameter Data (P) and Data Entry sections. Thus, Start(S) and Global(G)
		  sections are skipped.*/
		getline(f1, line1);

	}

	while (line1.at(72) == 'D') {
		if (line1.substr(48, 8) != "        ") t_matrix_pointer.push_back(stoi(line1.substr(48, 8)));
		else t_matrix_pointer.push_back(0);
		getline(f1, line1);
		if (line1.substr(32, 8) != "        ") form_number.push_back(stoi(line1.substr(32, 8)));
		else form_number.push_back(0);
		getline(f1, line1);
	}

	string next_line, line_n = "";
	string formatted_line = "";
	/*formatted_line will contain only the columns 1-64 of each line
	  as these columns contain the parameter data required.*/

	while (line1.at(72) == 'P') {
		getline(f1, next_line);
		formatted_line.append(line1.substr(0, 64));

		while (de_no(line1) == de_no(next_line)) {
			/*Some enities have paramters written in multiple lines.
			  But, their DE pointer number in same. Thus, parameters
			  of such entities are joined in a single formatted_string.*/
			formatted_line.append(next_line.substr(0, 64));
			line1 = next_line;
			getline(f1, next_line);
		}
		int pos = formatted_line.find(";");
		line_n = formatted_line.substr(0, pos + 1);
		vector <string> par = parameters(line_n, line1.substr(64, 8), del1, del2);
		//pass par to entity builder
		build_entity(par);
		/*cout << formatted_line << endl;*/
		formatted_line.clear();
		line1 = next_line;
	}
	sur_select_disp();

	if (line1.at(72) == 'T') cout << "END" << endl;
	//cout<<line_vec.size();

	cout << endl << endl << "type \"exit\" to exit" << endl;
	while (exit_ != "exit" && exit_ != "EXIT" && exit_ != "Exit") cin >> exit_;
}


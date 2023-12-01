#include <bits/stdc++.h>
#include "inc/Matrix.h"
using namespace std; 

typedef long long ll;
typedef vector<long long> vll;

/*
*	All constants details
*/
ll n_elem_points;
ll n_node_points;
double a;
double b;
double c1;

/*
*	All data structures details
*/
vector<vector<ll>> elem_points;
vector<vector<double>> node_points;
vector<vector<double>> ke_even;
vector<vector<double>> ke_odd; 
vector<vector<double>> k_global;
vector<vector<double>> k_global_with_rows_deleted;
vector<vector<double>> k_global_with_cols_also_deleted;
vector<ll> boundary_at_0;
vector<ll> boundary_at_spec_temp;
set<ll> nodes_to_be_deleted;
set<ll> nodes_with_spec_temp;
vector<vector<double>> cols_copied_from_k;
vector<double> rhs;
vector<double> nodes_with_spec_temp_vec;

/*
*	Taking input for fine.elem
*/
void input_fine_elem(){
	cin >> n_elem_points;
	elem_points.clear();
	elem_points.resize(n_elem_points, vector<ll>(4,0));

	for(int i = 0; i < n_elem_points; i++){
		for(int j = 0; j < 4; j++){
			cin >> elem_points[i][j];
		}
	}
}

/*
*	Taking input for fine.node
*/
void input_fine_node(){
	cin >> n_node_points;
	node_points.clear();
	node_points.resize(n_node_points, vector<double>(3,0));

	for(int i = 0; i < n_node_points; i++){
		for(int j = 0; j < 3; j++){
			cin >> node_points[i][j];
		}
	}
}

/*
*	Initialize a and b
*/
void init_a_and_b(){
	a = 0.06;
	b = 0.04;
	c1 = 1/(2*a*b);
}

/*
*	Initialize smaller(3x3) stiffness matrices Ke_even and Ke_odd
*/
void init_ke_even_and_ke_odd(){
	// clear and reshape
	ke_odd.clear();
	ke_even.clear();
	ke_odd.resize(3, vector<double>(3,0));
	ke_even.resize(3, vector<double>(3,0));

	// fill ke_odd
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			if((i == 0 || i == 2) && (j == 0 || j == 2)){
				ke_odd[i][j] += b*b;
			}
			if(i >= 1 && j >= 1){
				ke_odd[i][j] += a*a;
			}
			if((j == 2 && (i == 0 || i == 1)) || (i == 2 && (j == 0 || j == 1)))
			{
				ke_odd[i][j] *= -1;
			}
		}
	}

	// fill ke_even
	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			if(i <= 1 && j <= 1){
				ke_even[i][j] += a*a;
			}
			if(i >= 1 && j >= 1)
			{
				ke_even[i][j] += b*b;
			}
			if((i == 1 && (j == 0 || j == 2)) || (j == 1 && (i == 0 || i == 2))){
				ke_even[i][j] *= -1;
			}
		}
	}

	for(int i=0; i<3; i++){
		for(int j=0; j<3; j++){
			ke_even[i][j] *= c1;
			ke_odd[i][j] *= c1;
		}
	}
}

/*
*	Initialize k_global(121x121)
*/
void init_k_global(){
	//clear and resize
	k_global.clear();
	k_global.resize(122, vector<double>(122,0));

	// fill in values
	for(auto elem : elem_points){
		vector<ll> numbering = {elem[1], elem[2], elem[3]};

		if(elem[0]%2 == 0){
			// even elem, use ke_even
			for(int i=0; i<3; i++){
				for(int j=0; j<3; j++){
					k_global[numbering[i]][numbering[j]] += ke_even[i][j];
				}
			}

		} else {
			// odd elem, use ke_odd
			for(int i=0; i<3; i++){
				for(int j=0; j<3; j++){
					k_global[numbering[i]][numbering[j]] += ke_odd[i][j];
				}
			}
		}
	}
}

/*
*	Initialize boundary vectors, for T=0, and T=500cos(..)
*/	
void init_boundary_condition_vectors(){
	boundary_at_0.clear();
	boundary_at_spec_temp.clear();
	nodes_to_be_deleted.clear();
	nodes_with_spec_temp.clear();

	for(auto node : node_points){
		if(node[1] == 0.60){
			boundary_at_0.push_back(node[0]);
			nodes_to_be_deleted.insert(node[0]);
		}
		if(node[2] == 0.40){
			boundary_at_spec_temp.push_back(node[0]);
			nodes_to_be_deleted.insert(node[0]);
			nodes_with_spec_temp.insert(node[0]);
		}
	}
}

/*
*	Any rows in k_global related to dirichlet boundary cond.
*	has to be deleted.
*/	
void delete_rows_in_k_global(){
	k_global_with_rows_deleted.clear();
	k_global_with_rows_deleted.resize(0, vector<double>(122));
	ll row_no = 0;

	for(auto row : k_global){
		if(nodes_to_be_deleted.count(row_no)){
			row_no++;
			continue;
		}

		k_global_with_rows_deleted.push_back(row);
		row_no++;
	}
}

/*
*	Any cols with T=500cos() has to be copied as their values
*	will be used later
*/
void copy_spec_temp_cols(){
	ll rows = k_global_with_rows_deleted.size();
	ll cols = k_global_with_rows_deleted[0].size();

	for(int i=0; i<cols; i++){
		if(nodes_with_spec_temp.count(i)){
			vector<double> col_elems;
			col_elems.clear();
			for(int j=0; j<rows; j++){
				col_elems.push_back(k_global_with_rows_deleted[j][i]);
			}
			cols_copied_from_k.push_back(col_elems);
		}
	}
}

/*
*	Deleting the cols with dirichlet boundary cond.
*/
void delete_cols_also_from_k_global(){
	ll new_mat_size = k_global_with_rows_deleted.size();
	ll total_cols = k_global.size();
	k_global_with_cols_also_deleted.clear();
	k_global_with_cols_also_deleted.resize(new_mat_size);

	for(int i=0; i<new_mat_size; i++){
		int col_no = 0;
		while(col_no < total_cols){
			if(nodes_to_be_deleted.count(col_no)){
				col_no++;
				continue;
			}

			k_global_with_cols_also_deleted[i].push_back(k_global_with_rows_deleted[i][col_no]);
			col_no++;
		}
	}
}

double temp_at_boundary(double x){
	return 500*cos(2*acos(0.0)*x/1.2);
}

/*
*	Finding the RHS i.e. B in eq. Ax = B
*/
void find_rhs_vector(){
	vector<double> copy_of_nodes_with_spec_temp;
	for(auto si : nodes_with_spec_temp){
		copy_of_nodes_with_spec_temp.push_back(si);
	}

	for(int i=0; i<cols_copied_from_k.size(); i++){
		for(int j=0; j<cols_copied_from_k[0].size(); j++){
			cols_copied_from_k[i][j] *= -temp_at_boundary(node_points[copy_of_nodes_with_spec_temp[i]-1][1]);
		}
	}

	rhs.clear();
	rhs.resize(cols_copied_from_k[0].size(), 0);

	for(int i=0; i<cols_copied_from_k.size(); i++){
		for(int j=0; j<cols_copied_from_k[0].size(); j++){
			rhs[j] += cols_copied_from_k[i][j];
		}
	}
}


int main(){

	#ifndef ONLINE_JUDGE
    freopen("input.txt", "r", stdin);
    freopen("output.txt", "w", stdout);
	#endif

	input_fine_elem();
	input_fine_node();

	init_a_and_b();
	init_ke_even_and_ke_odd();
	init_k_global();
	init_boundary_condition_vectors();

	delete_rows_in_k_global();
	copy_spec_temp_cols();
	delete_cols_also_from_k_global();

	find_rhs_vector();

	vector <vector <double>> temp(rhs.size()-1, vector <double>(1, -1));
	for(int i=0;i<rhs.size()-1;i++){
		temp[i][0] = rhs[i+1];
	}

	int size = k_global_with_cols_also_deleted.size()-1;
	vector <vector <double>> mat(size, vector <double>(size, 0));

	for(int i=0;i<size;i++)	{
		for(int j=0;j<size;j++)		{
			mat[i][j] = k_global_with_cols_also_deleted[i+1][j+1];
		}
	}

	// For Ax = b equation our b is line 301 and A is line 304 
	Matrix b(temp, temp.capacity(), temp[0].capacity());
	Matrix matrix(mat, mat.capacity(), mat.capacity());
	matrix.printMatrix();
	Matrix inv(matrix.getInverse(), mat.capacity(), mat.capacity());
	Matrix finaltemp = inv*b;
	
	vector<double> temperature_values = finaltemp.copy_temperature_matrix();
	// finaltemp.printMatrix();
	vector<double> final_temperature_vector;

	int index = 0;
	for(int i=1; i<=n_node_points; i++){
		if(nodes_to_be_deleted.count(i)){
			final_temperature_vector.push_back(0);
		} else if(nodes_with_spec_temp.count(i)){
			final_temperature_vector.push_back(temp_at_boundary(node_points[boundary_at_spec_temp[i]-1][1]));
		} else {
			final_temperature_vector.push_back(temperature_values[index]);
			index++;
		}
	}
	for(auto i:ke_even)
	{
		for(auto j:i)
		{
			cout<<j<<" ";
		}
		cout<<endl;
	}
	for(auto i:ke_odd)
	{
		for(auto j:i)
		{
			cout<<j<<" ";
		}
		cout<<endl;
	}
	for(auto i:final_temperature_vector)
	{
		cout<<i<<" ";
	}
	cout<<endl;
	// final temperature vector is the final thing that should be printed in the output file
}
//acknowledgement: cs294-73 colella@eecs.berkeley.edu
#include "FEGrid.h"
#include "Matrix.h"

#include <iostream>
#include <set>

#define CONDUCTIVITY 1
typedef long long ll;
using namespace std;

vector<vector<double>> k_global_with_rows_deleted;
vector<vector<double>> k_global_with_cols_also_deleted;
vector<ll> boundary_at_0;
vector<ll> boundary_at_spec_temp;
set<ll> nodes_to_be_deleted;
set<ll> nodes_with_spec_temp;
vector<vector<double>> cols_copied_from_k;
vector<double> rhs;
vector<double> nodes_with_spec_temp_vec;

double temp_at_boundary(double x)
{
	return 500 * cos(2 * acos(0.0) * x / 1.2);
}

void init_boundary_condition_vectors(vector<vector<double>> node_points)
{
	boundary_at_0.clear();
	boundary_at_spec_temp.clear();
	nodes_to_be_deleted.clear();
	nodes_with_spec_temp.clear();

	for (auto node : node_points)
	{
		if (node[1] == 0.60)
		{
			boundary_at_0.push_back(node[0]);
			nodes_to_be_deleted.insert(node[0]);
		}
		if (node[2] == 0.40)
		{
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
void delete_rows_in_k_global(vector<vector<double>> k_global)
{
	k_global_with_rows_deleted.clear();
	k_global_with_rows_deleted.resize(0, vector<double>(122));
	ll row_no = 0;

	for (auto row : k_global)
	{
		if (nodes_to_be_deleted.count(row_no))
		{
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
void copy_spec_temp_cols()
{
	ll rows = k_global_with_rows_deleted.size();
	ll cols = k_global_with_rows_deleted[0].size();

	for (int i = 0; i < cols; i++)
	{
		if (nodes_with_spec_temp.count(i))
		{
			vector<double> col_elems;
			col_elems.clear();
			for (int j = 0; j < rows; j++)
			{
				col_elems.push_back(k_global_with_rows_deleted[j][i]);
			}
			cols_copied_from_k.push_back(col_elems);
		}
	}
}

/*
*	Deleting the cols with dirichlet boundary cond.
*/
void delete_cols_also_from_k_global(vector<vector<double>> k_global)
{
	ll new_mat_size = k_global_with_rows_deleted.size();
	ll total_cols = k_global.size();
	k_global_with_cols_also_deleted.clear();
	k_global_with_cols_also_deleted.resize(new_mat_size);

	for (int i = 0; i < new_mat_size; i++)
	{
		int col_no = 0;
		while (col_no < total_cols)
		{
			if (nodes_to_be_deleted.count(col_no))
			{
				col_no++;
				continue;
			}

			k_global_with_cols_also_deleted[i].push_back(k_global_with_rows_deleted[i][col_no]);
			col_no++;
		}
	}
}

/*
*	Finding the RHS i.e. B in eq. Ax = B
*/
void find_rhs_vector(vector<vector<double>> node_points)
{
	vector<double> copy_of_nodes_with_spec_temp;
	for (auto si : nodes_with_spec_temp)
	{
		copy_of_nodes_with_spec_temp.push_back(si);
	}

	for (int i = 0; i < cols_copied_from_k.size(); i++)
	{
		for (int j = 0; j < cols_copied_from_k[0].size(); j++)
		{
			cols_copied_from_k[i][j] *= -temp_at_boundary(node_points[copy_of_nodes_with_spec_temp[i] - 1][1]);
		}
	}

	rhs.clear();
	rhs.resize(cols_copied_from_k[0].size(), 0);

	for (int i = 0; i < cols_copied_from_k.size(); i++)
	{
		for (int j = 0; j < cols_copied_from_k[0].size(); j++)
		{
			rhs[j] += cols_copied_from_k[i][j];
		}
	}
}

int main(int argc, char **argv)
{
	if (argc != 2)
	{
		cout << "this program takes one argument that is the common name prefix of .node and .elem files ";
		cout << "file prefix" << endl;
		return 1;
	}
	string prefix(argv[1]);
	string nodeFile = prefix + ".node";
	string eleFile = prefix + ".elem";
	freopen("solution.txt", "w", stdout);
	FEGrid grid(nodeFile, eleFile);
	//insert code here to form the global stiffness matrix and the N.C. vector. Also, compute temperature at each node of the grid based on the stiffness matrix and NC vector. Refer to 5_10_Notes.pdf (for global stiffness matrix), 6_10_Notes.pdf, and 7_10_Notes.pdf (for 2D triangular elements)

	//write temp at each node into a file. The file should contain <number of nodes> number of lines.
	//On each line print two values in exactly this format: <node number> <blank_space> <temperature>.

	int interior_node_count = grid.getNumInteriorNodes();
	vector<vector<double>> k_global(interior_node_count, vector<double>(interior_node_count, 0));
	float grad[DIM];
	grid.gradient(grad, 200, 1);

	for (int i = 0; i < grid.getNumElts(); i++)
	{
		vector<vector<double>> b(VERTICES, vector<double>(DIM, 0));
		vector<vector<double>> b_transpose(DIM, vector<double>(VERTICES, 0));
		for (int localNode1 = 0; localNode1 < 3; localNode1++)
		{
			for (int localNode2 = 0; localNode2 < 3; localNode2++)
			{
				Node elemNode1 = grid.getNode(i, localNode1);
				Node elemNode2 = grid.getNode(i, localNode2);
				float localGradient1[DIM];
				float localGradient2[DIM];
				if (elemNode1.isInterior() && elemNode2.isInterior())
				{
					grid.gradient(localGradient1, i, localNode1);
					grid.gradient(localGradient2, i, localNode2);
					float gradientDot = 0;
					for (int i = 0; i < DIM; i++)
					{
						gradientDot += localGradient1[i] * localGradient2[i];
					}
					k_global[elemNode1.getInteriorNodeID()][elemNode2.getInteriorNodeID()] += grid.elementArea(i) * gradientDot;
				}
			}
		}
	}

	Matrix K_GLOBAL(k_global, interior_node_count, interior_node_count);
	vector<vector<double>> k_glo = K_GLOBAL.getMatrix();

	delete_rows_in_k_global(k_glo);
	delete_cols_also_from_k_global(k_glo);

	return 0;
}

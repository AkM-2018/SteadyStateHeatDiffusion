#pragma once
#include <iostream>
#include <vector>
#include <iomanip>
#include <cmath>
#define EPSILON 0.000001
class Matrix
{
private:
	std::vector<std::vector<double>> m_matrix;
	std::vector<std::vector<double>> m_identity; //to get the inverse matrix.
	//we use identity matrix to get inverse in gaussJordan.
	std::vector<std::vector<double>> m_copy; //to not change the orignal matrix..
	int m_rownum;
	int m_colnum;
	int m_rank;
	bool isInvertible;
	int m_iteration = 0;

public:
	Matrix(std::vector<std::vector<double>> matrix, const int &rownum, const int &colnum)
		: m_matrix(matrix), m_rownum(rownum), m_colnum(colnum)
	{
		m_copy.assign(m_matrix.begin(), m_matrix.end());
		m_identity.assign(rownum, std::vector<double>(rownum, 0));
		for (int i = 0; i < rownum; i++)
			m_identity[i][i] = 1;
		isInvertible = true;
		m_rank = m_rownum; //initialization of rank.
	}
	void printMatrix()
	{
		for (const auto &i : m_matrix)
		{
			for (const auto &j : i)
				std::cout << std::setprecision(8) << j << " ";
			std::cout << std::endl;
		}
	}

	int getM()
	{
		return m_rownum;
	}
	int getN()
	{
		return m_colnum;
	}
	void printCurrentMatrix()
	{
		for (int i = 0; i < m_rownum; i++)
		{
			for (const auto &e : m_copy.at(i))
				std::cout << std::setprecision(4) << e << " ";
			std::cout << " ";
			for (const auto &e2 : m_identity.at(i))
				std::cout << std::setprecision(4) << e2 << " ";
			std::cout << "\n";
		}
		//Shows How the matrix changes as the iteration goes, you may comment if you dont need
	}
	void GaussJordan()
	{
		rowEchleon();	  // First, changes the matrix in to row Echleon Form
		reducedEchleon(); //Second, using the matrix used in RowEchleon, changes it into reduced Echelon.
		if (!isInvertible)
			std::cout << "Can't Invert matrix.."
					  << "your rank of the matrix is " << m_rank << std::endl;
	}
	void rowEchleon()
	{
		for (int i = 0; i < m_rownum; i++)
		{
			int count = i + 1;
			//checking if swapping is available..
			//i.e if matrix[0][0]is 0 then it needs to be swapped to nonzero row..
			while ((std::abs(m_copy[i][i]) < EPSILON) && count < m_rownum)
			{
				if (std::abs(m_copy[count][i] < EPSILON))
				{
					count++;

					if (count == m_rownum)
						isInvertible = false;
				}
				else
				{
					rowSwap(i, count);
					count++;
					break;
				}
			}
			//after checking if it could not be swapped, isInvertibel is false.
			if (std::abs(m_copy[i][i]) < 0.000001)
			{
				isInvertible = false;
				return;
			}
			//after swapping change the diagonals to 1
			if (std::abs(m_copy[i][i] - 1) > 0.00001)
			{
				double den = m_copy[i][i];
				for (int j = 0; j < m_rownum; j++)
				{
					//making Diagonal elements to 1..
					m_identity[i][j] = m_identity[i][j] / den;
				}
				for (int j = 0; j < m_rownum; j++)
				{
					//making Diagonal elements to 1..
					m_copy[i][j] = m_copy[i][j] / den;
				}
			}
			//for every column, make every row except the down-half the diagnoal ones 0 by doing row element operation.
			for (int j = i + 1; j < m_rownum; j++)
			{
				double div = -(m_copy[j][i] / m_copy[i][i]);
				for (int k = i; k < m_rownum; k++)
				{
					m_copy[j][k] = m_copy[i][k] * div + m_copy[j][k];
				}
				for (int l = (m_rownum - 1); l >= 0; l--)
				{
					m_identity[j][l] = m_identity[i][l] * div + m_identity[j][l];
				}
			} //This were the statement for changes to row-echeleon form..
			  //std::cout << "Iteration " << m_iteration++ << std::endl;
			  //printCurrentMatrix();
			  //std::cout << std::endl;
		}
	}
	void reducedEchleon()
	{
		// Even though the matrix is not invertible, to get the rank, do reduced Echleon for every matrix.
		for (int i = m_rownum - 1; i >= 0; i--)
		{
			if (abs(m_copy[i][i]) < EPSILON)
			{
				m_rank--;
				for (int r = 0; r < m_rownum; r++)
				{
					m_copy[r][i] = 0;
				}
			} ///reduces the rank  one by one..

			else
			{
				for (int j = i - 1; j >= 0; j--)
				{
					double div = -(m_copy[j][i] / m_copy[i][i]);
					updateRow(m_copy, j, i, div);
					if (isInvertible)
					{
						updateRow(m_identity, j, i, div);
						//Make every row except the diagonal to zero;
					}
				}
			}
			//std::cout << "Iteration : " << m_iteration++ << std::endl;
			//printCurrentMatrix();
			//std::cout << std::endl;
		}
	}
	void printInverse()
	{
		GaussJordan();
		if (isInvertible)
		{
			std::cout << "The Matrix is Invertible" << std::endl;
			std::cout << "Your Inversed matrix is " << std::endl;
			for (const auto &i : m_identity)
			{
				for (const auto &j : i)
					std::cout << std::setprecision(3) << j << " ";
				std::cout << std::endl;
			}
		}
		else
		{
			std::cout << "Inverse matrix is not availalbe." << std::endl;
		}
	}
	void rowSwap(const int &row1, const int &row2)
	{
		//swaps rows of the matrix.
		std::vector<double> temp;
		std::vector<double> temp2;
		temp.assign(m_copy.at(row1).begin(), m_copy.at(row1).end());
		m_copy.at(row1) = m_copy.at(row2);
		m_copy.at(row2) = temp;
		temp2.assign(m_identity.at(row1).begin(), m_identity.at(row1).end());
		m_identity.at(row1) = m_identity.at(row2);
		m_identity.at(row2) = temp2;
	}

	void updateRow(std::vector<std::vector<double>> &vec, const int &cur_rownum,
				   const int &chg_rownum, const double &div)
	{
		for (int l = 0; l < m_rownum; l++)
			vec[cur_rownum][l] = vec[chg_rownum][l] * div + vec[cur_rownum][l];
	}
	int getRank() { return m_rank; }
	std::vector<std::vector<double>> getMatrix()
	{
		return m_matrix;
	}
	std::vector<std::vector<double>> getInverse()
	{
		GaussJordan();
		return m_identity;
	}
	Matrix operator*(Matrix &b)
	{
		int n = b.getN();
		int m = m_rownum;
		std::vector<std::vector<double>> result(m, std::vector<double>(n, 0));
		for (int i = 0; i < m; i++)
		{
			for (int k = 0; k < m_colnum; k++)
			{
				for (int j = 0; j < n; j++)
				{
					result[i][j] += m_matrix[i][k] * b.getMatrix()[k][j];
				}
			}
		}
		Matrix c(result, m, n);
		return c;
	}
};
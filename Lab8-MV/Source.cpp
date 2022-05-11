#include <iostream>
#include "windows.h"
#include <vector>
#include "string"
#include <cmath>

using namespace std;

class Matrix {
private:
	vector<vector<double>> matrix;
	vector<int> width_form;
	int m = 0, n = 0;

	void Format() {
		width_form.clear();
		int width, buf_width, width_null;
		string str_width;
		width = 1;
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < m; j++) {
				str_width = to_string(matrix[j][i]);
				for (int i = str_width.size() - 1; i >= 0; i--) {
					if (str_width[i] == '0' || str_width[i] == ',') {
						str_width.erase(i);
					}
					else break;
				}
				buf_width = str_width.size();
				if (width < buf_width) width = buf_width;
			}
			width_form.push_back(width);
		}
	}
	void Size() {
		m = matrix.size();
		n = matrix[0].size();
		for (int i = 1; i < m; i++) {
			if (matrix[i].size() != n) {
				n = 0;
				break;
			}
		}
	}
public:
	void Input_Matrix(int m, int n)
	{
		double number;
		vector<double> str_matrix;
		cout << "Enter the [" << m << "," << n << "] matrix:" << endl;
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cin >> number;
				str_matrix.push_back(number);
			}
			matrix.push_back(str_matrix);
			str_matrix.clear();
		}
	}
	void Output_Matrix() {
		Size();
		if (n == 0 || m == 0) {
			return;
		}
		Format();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cout.width(width_form[j]);
				cout << matrix[i][j] << " ";
			}
			cout << endl;
		}
		cout << endl;
	}
	void Output_SLAU() {
		Size();
		if (n == 0 || m == 0) {
			return;
		}
		Format();
		for (int i = 0; i < m; i++) {
			for (int j = 0; j < n; j++) {
				cout.width(width_form[j]);
				cout << matrix[i][j] << " ";
				if (j == n - 2) cout << "| ";
			}
			cout << endl;
		}
		cout << endl;
	}
	vector<vector<double>> get_Matrix() {
		return matrix;
	}
	void set_SLAU(vector<vector<double>> slau) {
		matrix = slau;
	}
};

class Determinant {
private:
	vector<vector<double>> matrix_det;
	double det = 0;
	Matrix M;

	double Accur(double n) {
		double buf = abs(n);
		if (buf != 0) {
			if (fmod(buf, 0.001) < 1) {
				n = round(n * 1000) / 1000;
			}
		}
		return n;
	}
	double Det_2x(vector<vector<double>> matrix_2x) {
		double right = matrix_2x[0][0] * matrix_2x[1][1];
		double left = matrix_2x[0][1] * matrix_2x[1][0];
		return Accur(right - left);
	}
	double Det_3x(vector<vector<double>> matrix_3x) {
		double right = 0;
		double buf = 1, buf1 = 1;
		for (int i = 0; i < 3; i++) {
			buf *= matrix_3x[i][i];
		}
		right += buf;
		buf = 1;
		for (int i = 0; i < 3; i++) {
			if (i + 1 < 3) {
				buf *= matrix_3x[i][i + 1];
				buf1 *= matrix_3x[i + 1][i];
			}
			else {
				buf *= matrix_3x[i][0];
				buf1 *= matrix_3x[0][i];
			}
		}
		right += buf;
		right += buf1;
		double left = 0;
		buf = matrix_3x[0][2] * matrix_3x[2][0] * matrix_3x[1][1];
		left += buf;
		buf = matrix_3x[0][1] * matrix_3x[1][0] * matrix_3x[2][2];
		left += buf;
		buf = matrix_3x[1][2] * matrix_3x[2][1] * matrix_3x[0][0];
		left += buf;
		return Accur(right - left);
	}
	double Det_mx(vector<vector<double>> matrix_mx) {
		if (matrix_mx.size() == 2) return Det_2x(matrix_mx);
		else if (matrix_mx.size() == 3) return Det_3x(matrix_mx);
		else {
			double det_mx = 0;
			for (int i = 0; i < matrix_mx.size(); i++) {
				vector<vector<double>> buf_m = matrix_mx;
				buf_m.erase(buf_m.begin());
				int buf_size = buf_m.size();
				for (int buf_index = 0; buf_index < buf_size; buf_index++) {
					buf_m[buf_index].erase(buf_m[buf_index].begin() + i);
				}
				if (i % 2 == 0) {
					det_mx += matrix_mx[0][i] * Det_mx(buf_m);
				}
				else det_mx -= matrix_mx[0][i] * Det_mx(buf_m);
			}
			return Accur(det_mx);
		}
	}
public:
	void Input_Det_Matrix(vector<vector<double>> matrix) {
		matrix_det = matrix;
	}
	void Decision() {
		det = Det_mx(matrix_det);
		cout << "det(A) = " << det << endl;
	}
	double get_det() {
		return det;
	}
};

class SLAU {
private:
	vector<vector<double>> matrix_slau, matrix_a;
	vector<double> X;
	Matrix M;
	int m, k = 0;
	void Accur_str(int str_index) {
		double buf;
		vector<double> str_slau = matrix_slau[str_index];
		for (int i = 0; i < str_slau.size(); i++) {
			buf = abs(str_slau[i]);
			if (buf != 0) {
				if (fmod(buf, 0.001) < 1) {
					str_slau[i] = round(str_slau[i] * 1000) / 1000;
				}
			}
		}
		matrix_slau[str_index] = str_slau;
	}
	double Accur(double n) {
		double buf = abs(n);
		if (buf != 0) {
			if (fmod(buf, 0.001) < 1) {
				n = round(n * 1000) / 1000;
			}
		}
		return n;
	}
	void Output_SLAU() {
		M.set_SLAU(matrix_slau);
		M.Output_SLAU();
	}
	void Output_a() {
		M.set_SLAU(matrix_a);
		M.Output_SLAU();
	}
	void OutputX() {
		cout << "k = " << k;
		for (int i = 0; i < m; i++) {
			cout << " x" << i + 1 << " = " << X[i];
		}
		cout << endl;
	}
	bool test_Accur(vector<double> x_buf) {
		double accur_max = 0;
		double buf = 0;
		for (int i = 0; i < m; i++) {
			buf = abs(X[i] - x_buf[i]);
			if (accur_max < buf) accur_max = buf;
		}
		if (accur_max < 0.001) return true;
		else return false;
	}
	void X_Record() {
		OutputX();
		double buf_x = 0;
		vector<double> buf_;
		for (int i = 0; i < m; i++) {
			buf_.push_back(0);
			for (int j = 0; j < i; j++) {
				buf_[i] += matrix_a[i][j] * buf_[j];
			}
			for (int j = i; j < m; j++) {
				buf_[i] += matrix_a[i][j] * X[j];
			}
			buf_[i] += matrix_a[i][m];
		}
		k++;
		if (test_Accur(buf_)) {
			X = buf_;
			for (int i = 0; i < m; i++) {
				X[i] = Accur(X[i]);
			}
			return;
		}
		else {
			X = buf_;
			X_Record();
		}
	}
	bool str_diagonal(int index_i) {
		double a_i = abs(matrix_slau[index_i][index_i]);
		for (int i = 0; i < m; i++) {
			if (i != index_i) a_i -= abs(matrix_slau[index_i][i]);
		}
		if (a_i < 0) return false;
		return true;
	}
	bool str_a(vector<double> a_buf) {
		double a_ = 0;
		for (int i = 0; i < m; i++) {
			a_ += abs(a_buf[i]);
		}
		if (a_ >= 1) return false;
		return true;
	}

	bool Test1() {
		bool test = true;
		for (int i = 0; i < m; i++) {
			if (!str_diagonal(i)) {
				test = false;
				break;
			}
		}
		return test;
	}
	bool Create_a() {
		bool test = true;
		double aii;
		for (int i = 0; i < m; i++) {
			vector<double> buf;
			aii = matrix_slau[i][i];
			for (int j = 0; j < m; j++) {
				if (j == i) buf.push_back(0);
				else buf.push_back(-(matrix_slau[i][j] / aii));
			}
			buf.push_back((matrix_slau[i][m] / aii));
			if (!str_a(buf)) test = false;
			matrix_a.push_back(buf);
		}
		return test;
	}
public:
	bool Input_SLAU(vector<vector<double>> matrix_A, vector<vector<double>> matrix_B) {
		if (matrix_A.size() != matrix_B.size() || matrix_B[0].size() != 1) {
			return false;
		}
		matrix_slau = matrix_A;
		for (int i = 0; i < matrix_A.size(); i++) {
			matrix_slau[i].push_back(matrix_B[i][0]);
		}
		m = matrix_A.size();
		if (!Test1()) {
			cout << "Diagonals do not prevail." << endl;
			return false;
		}
		Output_SLAU();
		if (!Create_a()) {
			cout << "The eigenvalues of the matrix a >= 1" << endl;
			Output_a();
			return false;
		}
		Output_a();
		for (int i = 0; i < m; i++) {
			X.push_back(matrix_a[i][m]);
		}
		return true;
	}
	void Decision() {
		X_Record();
	}
	void Rezult() {
		cout << "Answer: ";
		for (int i = 0; i < X.size(); i++) {
			cout << "x" << i + 1 << " = " << X[i] << ", ";
		}
		cout << endl;
	}

};

class Task {
private:
	Matrix M;
	Matrix M1;
	SLAU Slau;
	Determinant Det;
	int m, n;
	void Create_SLAU() {
		M.Input_Matrix(m, m);
		Det.Input_Det_Matrix(M.get_Matrix());
		n = 1;
		M1.Input_Matrix(m, n);
	}
public:
	void SLAU_Task() {
		cout << "Input the size of the matrix." << endl << "m = ";
		cin >> m;
		if (m < 3 || m > 49) {
			cout << "Incorrect size of the matrix." << endl;
			return;
		}
		Create_SLAU();
		Det.Decision();
		if (Det.get_det() == 0) {
			cout << "Error! The matrix det = 0." << endl;
			return;
		}
		bool flag = Slau.Input_SLAU(M.get_Matrix(), M1.get_Matrix());
		if (!flag) {
			cout << "Incorrect matrix." << endl;
			return;
		}
		Slau.Decision();
		Slau.Rezult();
	}
};

int main() {
	setlocale(LC_ALL, "Rus");
	SetConsoleCP(1251);
	int ans, exit = 1;
	while (exit == 1) {
		Task T;
		cout << "1.SLAU" << endl << "2.Exit" << endl << "Choose a way:" << endl;
		cin >> ans;
		switch (ans)
		{
		case 1:
			T.SLAU_Task();
			break;
		case 2:
			exit = 0;
			break;
		default:
			cout << "This task does not exist" << endl;
			break;
		}
	}
	system("pause");
	return 0;
}
/*

4.6 1.8 -1.7
2.7 -5.6 1.9
1.5 0.5 3.3
3.8
0.4
-1.6

*/
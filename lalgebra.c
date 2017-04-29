/* Linear Algebra Library
 * Developed by Joao Francisco B. S. Martins <joaofbsm@dcc.ufmg.br>
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

#include "lalgebra.h"

void parse_input(FILE* input, double** matrix, int m, int n) {
	int i, j, input_size;

	i = 0;

    input_size = (2 + (4 * m) + (m * (8 * n))); // Approximated size of input
    char* input_matrix = malloc(input_size * sizeof(char)); // Unformated LP
    fgets(input_matrix, input_size, input); 

    char* row = strtok(input_matrix, "{"); // Parsed first row

    // Parsing of rows to get elements and fill matrix
    while(row != NULL) { 
    	char* element = row;
		while (*element) {
			if ((element[0] == '-' && isdigit(element[1])) || isdigit(element[0])) { // Check if char is a number
		    	double val = strtol(element, &element, 10); // Convert char to number in base 10
		    	matrix[i][j] = val; // Add element to the corresponding position in the matrix
		    	j++;
//		    	printf("%g ", val); 
			} 
			else {
		    	element++;
			}
		}
//		printf("\n");

    	row = strtok(NULL, "{");
    	i++;
    	j = 0;
    }
//    printf("\n");
}

// Allocate matrix of dimensions m x n in the pointer of pointers **matrix.
double** allocate_matrix(int m, int n) {
	double** matrix;
	int i;

	matrix = malloc(m * sizeof(double*));
	for(i = 0; i < m; i++) {
		matrix[i] = malloc(n * sizeof(double)); 
	}

	return matrix;
}

// Creates identity matrix of dimensions m x m.
double** identity(int m) {
	double** identity_matrix;
	int i, j;

	identity_matrix = allocate_matrix(m, m);
	for(i = 0; i < m; i++) {
		for(j = 0; j < m; j++) {
			if(i == j) {
				identity_matrix[i][j] = 1; 
			}
			else {
				identity_matrix[i][j] = 0; 
			}
		}
	}

	return identity_matrix;
}

void print_matrix(double** matrix, int m, int n) {
	int i, j;

	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
			if(matrix[i][j] < 0 || matrix[i][j] > 9) {
				printf("%g ", round(matrix[i][j] * 100000) / 100000);
			}
			else {
				printf(" %g ", round(matrix[i][j] * 100000) / 100000);
			}
		}
		printf("\n");
	}
	printf("\n");
}

// Copy the value of every element in the original matrix to the new one.
void copy_matrix(double** original, double** copy, int m, int n) {
	int i, j;

	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
			copy[i][j] = original[i][j];
		}
	}
}

/* Insert source matrix inside target matrix in the submatrix comprehended by the integer offsets(from_row, to_row, etc). 
 * The dimensions must match. The received indexes start at 1 so we convert them. */
void insert_matrix(double** source, double** target, int from_row, int to_row, int from_column, int to_column) {
	int i, j, m, n;

	m = 0;
	n = 0;

	for(i = (from_row - 1); i < to_row; i++) {
		for(j = (from_column - 1); j < to_column; j++) {
//			printf("Source: [%d][%d] Target: [%d][%d] Value: %f\n", m, n, i, j, source[m][n]);
			target[i][j] = source[m][n];
			n++;
		}
		m++;
		n = 0;
	}
//	printf("END\n");
}

/* Offers a way to make linear operations. If sum_to is -1, replaces the actual line. m stands for the row 
 * to operate on and n stands for the dimension of columns. The index for sum_to begins at 0 */
void operate_on_rows(double** matrix, int m, int n, double multiply_by, int sum_to) {
	double* new_row;
	int j;

	new_row = malloc(n * sizeof(double)); 

	//Solves problem for really small negative numbers causing -0 to be printed and for really big ratios to appear
	if(fabs(multiply_by) < 0.00001) {
		multiply_by = 0;
	}

//	printf("Multiply row %d by %g and sum to row %d\n", m, multiply_by, sum_to);

//	printf("The new row is:\n");
	for(j = 0; j < n; j++) {
		if(matrix[m][j] != 0) {
			new_row[j] = matrix[m][j] * multiply_by;
		}
//		printf("%g ", new_row[j]);
	}
//	printf("\n\n");

	if(sum_to == -1) {
		for(j = 0; j < n; j++) {
			matrix[m][j] = new_row[j]; // TODO - Check if matrix[m] = new_row, or its variants, works
			if(fabs(matrix[m][j]) < 0.00001) {
				matrix[m][j] = 0;
			}
		}	
	}
	else {
		for(j = 0; j < n; j++) {
			matrix[sum_to][j] += new_row[j];
			if(fabs(matrix[sum_to][j]) < 0.00001) {
				matrix[sum_to][j] = 0;
			}
		}	
	}
}

void operate_on_columns(double** matrix, int m, int n, double multiply_by, int sum_to) {
	double* new_column;
	int i;

	new_column = malloc(n * sizeof(double)); 

	//Solves problem for really small negative numbers causing -0 to be printed and for really big ratios to appear
	if(fabs(multiply_by) < 0.00001) {
		multiply_by = 0;
	}

//	printf("Multiply col %d by %g and sum to row %d\n", m, multiply_by, sum_to);

	for(i = 0; i < m; i++) {
		if(matrix[i][n] != 0) {
			new_column[i] = matrix[i][n] * multiply_by;
		}
	}

	if(sum_to == -1) {
		for(i = 0; i < m; i++) {
			matrix[i][n] = new_column[i];
			if(fabs(matrix[i][n]) < 0.00001) {
				matrix[i][n] = 0;
			}
		}	
	}
	else {
		for(i = 0; i < m; i++) {
			matrix[i][sum_to] += new_column[i];
			if(fabs(matrix[i][sum_to]) < 0.00001) {
				matrix[i][sum_to] = 0;
			}
		}	
	}
}

// Format the LP to the Standard Equality Form
double** format_sef(double** original_matrix, int m, int n) {
	double** new_matrix;
	int i, j, new_columns, new_n;
	
	new_columns = m - 1;
	new_n = n + new_columns;

	new_matrix = allocate_matrix(m, (n + new_columns));

	// Set the values in the submatrix that will register the operations
	for(j = 0; j < new_columns; j++) {
		new_matrix[0][j] = 0; // Set c row to 0
	}

	// Insert the original matrix in the beginning of the new one without the last column
	insert_matrix(original_matrix, new_matrix, 1, m, 1, (n - 1));

	// Adds the identity matrix in the correct position
	insert_matrix(identity(new_columns), new_matrix, 2, m, n, (n + new_columns)); 

	// Adds the last column of the original LP to the new one
	for(i = 0; i < m; i++) {
		new_matrix[i][new_n - 1] = original_matrix[i][n - 1];
	}

	free(original_matrix);

	return new_matrix;
}

// Negates the entries in c for the tableau
void format_tableau(double** matrix, int m, int n) {
	int i;
	for(i = 0; i < n; i++) {
		if(matrix[0][i] != 0) {
			matrix[0][i] = matrix[0][i] * -1;
		}
	}
}

// Create a new matrix that consists of the operation register submatrix added to the original one.
double** add_operations_register(double** original_matrix, int m, int n) {
	double** new_matrix;
	int j, new_columns;

	new_columns = m - 1; // Number of rows added

	new_matrix = allocate_matrix(m, (n + new_columns));

	// Set the values in the submatrix that will register the operations
	for(j = 0; j < new_columns; j++) {
		new_matrix[0][j] = 0; // Set c row to 0
	}
	insert_matrix(identity(new_columns), new_matrix, 2, m, 1, new_columns); // Adds the identity matrix in the correct position

	// Insert the original matrix as a submatrix of the new one
	insert_matrix(original_matrix, new_matrix, 1, m, m, (n + new_columns));

	free(original_matrix);

	return new_matrix;
}

void make_b_non_negative(double** matrix, int m, int n) {
	int i;

	// First check if b > 0
	for(i = 1; i < m; i++) {
		if(matrix[i][n - 1] < 0) {
			operate_on_rows(matrix, i, n, -1, -1);
		}
	}
}

void make_c_non_negative(double** matrix, int m, int n) {
	int j;

	// First check if c > 0
	for(j = 0; j < (n - 1); j++) {
		if(matrix[0][j] < 0) {
			operate_on_columns(matrix, m, j, -1, -1);
		}
	}
}

double** create_auxiliar_lp(double** matrix, int m, int n) {
	double** auxiliar_lp;
	double** copied_matrix;
	int i, j, auxiliar_n;

	auxiliar_n = n + m - 1;
	auxiliar_lp = allocate_matrix(m, auxiliar_n);

	copied_matrix = allocate_matrix(m, n); 
	copy_matrix(matrix, copied_matrix, m, n); // Holds the values of the original matrix but doesn't mess with the original data in any sense

	make_b_non_negative(copied_matrix, m, n);
	
	// Fulfill the auxiliar lp without the last column(thats why to_column equals n - 2)
	insert_matrix(copied_matrix, auxiliar_lp, 1, m, 1, (n - 1));	

	// Adds the last column to the auxiliar lp
	for(i = 0; i < m; i++) {
		auxiliar_lp[i][auxiliar_n - 1] = copied_matrix[i][n - 1];
	}

	// Creates the first row of the auxiliar LP in the correct form
	for(j = 0; j < auxiliar_n; j++) { 
		if(j >= (auxiliar_n - m) && j < (auxiliar_n - 1)) {
			auxiliar_lp[0][j] = 1;
		} 
		else {
			auxiliar_lp[0][j] = 0; 
		}
	} 

	insert_matrix(identity(m - 1), auxiliar_lp, 2, m, (auxiliar_n - m + 1), (auxiliar_n - 1)); // Adds the identity matrix in the correct position

	/*
	for(i = 0; i < m; i++) {
		for(j = 0; j < auxiliar_n; j++) {
			printf("%g ", round(auxiliar_lp[i][j] * 100000) / 100000);
		}
		printf("\n");
	}*/



	return auxiliar_lp;
}

void set_initial_base(double** matrix, int m, int n, int* base) {
	int b, i;
	b = (n - 1 - (m - 1));
    for(i = 0; i < (m - 1); i++) {
    	base[i] = b;
    	b++;
    }
}

// Format the LP to the Canonical Form
void format_canonical(double** matrix, int m, int n, int* base) {
	/* Function flow *
	 * - Check if matrix[i + 1][base[i]] == 1
	 *   - If not, make it be(operate_on_rows)
	 * - Check if all the other rows for this collumn are 0
	 *   - If not, make it be
	 */
	int i, j;

	for(i = 0; i < (m - 1); i++) { // Goes through all the basic columns
		if(matrix[i + 1][base[i]] != 1) { 
			if(matrix[i + 1][base[i]] == 0) {
				// TODO - Find row(beginning from 1) where this column is not zero and add it to the current one
			}
			operate_on_rows(matrix, (i + 1), n, (1 / matrix[i + 1][base[i]]), -1);
		}
		for(j = 0; j < m; j++) {
			if(j != (i + 1)) {
				if(matrix[j][base[i]] != 0) {
					operate_on_rows(matrix, (i + 1), n, (-1 * (matrix[j][base[i]] / matrix[i + 1][base[i]])), j);
				}
			} 
		}
	}
}

// Needs to return the column, the row and the base that is leaving
int primal_next_base(double** matrix, int m, int n, int* base_row, int* base_column) {
	int i, j;
	double min_ratio, row_ratio;

	min_ratio = 999999;

//	printf("\n");

	for(j = (m - 1); j < (n - 1); j++) {
		if(matrix[0][j] < 0) {
			*base_column = j;
//			printf("Chosen column to enter the base was %d\n", j);
			for(i = 1; i < m; i++) {
				// b will never be negative after primal simplex starts to run, so, for a valid ratio, we need a positive number that is not zero
				if(matrix[i][j] > 0) {
					row_ratio = matrix[i][n - 1] / matrix[i][j];
//					printf("Ratio of row %d and column %d is %g\n", i, j, row_ratio);
					if(row_ratio < min_ratio) {
						min_ratio = row_ratio;
						*base_row = i;
					}
				}
			}
			if(min_ratio == 999999) {
				return j; // LP is unbounded
			}
			else {
				return 0; // Goes to the next round of simplex
			}
		}
	}

	return -1; // LP is optimal
}

// If return = -1 the LP is optimal and if return > 0 it's unbounded and the return value = the column where we can get the certificate
int primal_simplex(double** matrix, int m, int n, int* base, int print_output) {	
	int i, result, new_base_row, new_base_column;

	i = 0; // REMOVE THIS IN THE FINAL VERSION

	make_b_non_negative(matrix, m, n);

//	printf("Starting Primal Simplex with LP\n\n");
//	print_matrix(matrix, m, n);
//	printf("____________________________________________________\n\n");

	while(1) {

		new_base_row = 0;
		new_base_column = 0;

		// First we need to present the LP in the canonical form
		format_canonical(matrix, m, n, base);

//		printf("Canonical form for bases ");
//		for(i = 0; i < (m - 1); i++) {
//			printf("%d ", base[i]);
//		}
//		printf("\n\n");
//		print_matrix(matrix, m, n);

		if(print_output) {
			print_output_matrix(matrix, m, n);
		}

		// Find the next base for the primal simplex
		result = primal_next_base(matrix, m, n, &new_base_row, &new_base_column);

		if(result != 0) {
//			printf("Result %d\n\n", result);
			return result;
		}

//		printf("The new base is column %d and row %d\n\n", new_base_column, new_base_row);

		base[new_base_row - 1] = new_base_column;

//		printf("____________________________________________________\n\n");
//		getchar();
	}
}

int dual_next_base(double** matrix, int m, int n, int* base_row, int* base_column) {
	int i, j;
	double min_ratio, row_ratio;

	min_ratio = 999999;

//	printf("\n");

	for(i = 1; i < m; i++) {
		if(matrix[i][n - 1] < 0) {
			*base_row = i;
//			printf("Chosen row to enter the base was %d\n", i);
			for(j = (m - 1); j < (n - 1); j++) {
				// c will never be negative after primal simplex starts to run, so, for a valid ratio, we need a negative number that is not zero
				if(matrix[i][j] < 0) {
					row_ratio = matrix[0][j] / (-1 * matrix[i][j]);
//					printf("Ratio of row %d and column %d is %g\n", i, j, row_ratio);
					if(row_ratio < min_ratio) {
						min_ratio = row_ratio;
						*base_column = j;
					}
				}
			}
			if(min_ratio == 999999) {
				return j; // LP is unbounded
			}
			else {
				return 0; // Goes to the next round of simplex
			}
		}
	}

	return -1; // LP is optimal
}

int dual_simplex(double** matrix, int m, int n, int* base, int print_output) {
	int i, result, new_base_row, new_base_column;

	i = 0; // REMOVE THIS IN THE FINAL VERSION

	make_c_non_negative(matrix, m, n);

//	printf("Starting Dual Simplex with LP\n\n");
//	print_matrix(matrix, m, n);
//	printf("____________________________________________________\n\n");

	while(1) {

		new_base_row = 0;
		new_base_column = 0;

		// First we need to present the LP in the canonical form
		format_canonical(matrix, m, n, base);

//		printf("Canonical form for bases ");
//		for(i = 0; i < (m - 1); i++) {
//			printf("%d ", base[i]);
//		}
//		printf("\n\n");
//		print_matrix(matrix, m, n);

		if(print_output) {
			print_output_matrix(matrix, m, n);
		}

		// Find the next base for the primal simplex
		result = dual_next_base(matrix, m, n, &new_base_row, &new_base_column);

		if(result != 0) {
//			printf("Result %d\n\n", result);
			return result;
		}

//		printf("The new base is column %d and row %d\n\n", new_base_column, new_base_row);

		base[new_base_row - 1] = new_base_column;

//		printf("____________________________________________________\n\n");
//		getchar();
	}
}

double* get_primal_optimal_solution(double** matrix, int m, int n, int* base) {
	double* vector;
	int i;

	vector = malloc((n - 1 - (m - 1)) * sizeof(double));

	for(i = 0; i < (n - 1 - (m - 1)); i++) {
		vector[i] = 0;
	}

	// Assigns the value of b to the columns in the solution that correspond to the columns in the base
	for(i = 0; i < (m - 1); i++) {
		vector[base[i] - (m - 1)]  = matrix[i + 1][n - 1];
	}

	return vector;
}

// Generate both infeasibility and optimality certificates
double* get_dual_optimal_solution(double** matrix, int m) {
	double* vector;
	int i;

	vector = malloc((m - 1) * sizeof(double));

	for(i = 0; i < (m - 1); i++) {
		vector[i] = matrix[0][i];
	}

	return vector;
}

double* generate_unboundedness_certificate(double** matrix, int m, int n, int column, int* base) {
	double* vector;
	int i;

	vector = malloc((n - 1 - (m - 1)) * sizeof(double));

	for(i = 0; i < (n - 1 - (m - 1)); i++) {
		vector[i] = 0;
	}

	// Assigns the value of b to the columns in the solution that correspond to the columns in the base
	for(i = 0; i < (m - 1); i++) {
		if(matrix[i + 1][column] != 0) {
			vector[base[i] - (m - 1)] = -1 * matrix[i + 1][column];
		}
	}

	// Column that shows unboundedness will be 1 to make it easier to create the rest of the certificate
	vector[column - (m - 1)] = 1;

	return vector;
}

// TODO - Change this to fprintf
void print_output_vector(double* vector, int n) {
	int i;

	printf("{");
	for(i = 0; i < n; i++) {
		printf("%g", round(vector[i] * 100000) / 100000);
		if(i != (n - 1)) {
			printf(", ");
		}
	}
	printf("}");
}

void print_output_matrix(double** matrix, int m, int n) {
	int i;

	printf("{");
	for(i = 0; i < m; i++) {
		print_output_vector(matrix[i], n);
		if(i != (m - 1)) {
			printf(", ");
		}
	}
	printf("}\n");
}
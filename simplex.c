#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/*
CERTIFICATES

 - INVALID: Solution of the auxiliar PL's dual(Use add_operations_register in the auxiliar PL 
   and the corresponding c collumns are going to be the solution and the certificate as well)

 - UNBOUNDED: Find vector z such that:
	
		 - z >= 0	
         - Az = 0
         - cTz > 0

 - OPTIMAL: Solution of the dual -> it's in the first row and the collumns correspondent to the operation_register
*/

/*

- In the beginning God created the heavens and the earth. Then we need to put the LP in Standard Equality Form. 

- First call of solve_lp is ALWAYS to solve an auxiliar LP. The tableau for the auxiliar LP will always be in the final form after that.
  The optimal value will be the first element of the last collum and the optimal solution for the dual of the auxiliar LP will be the 
  first m indexes of the first line.

- We need to save the solution of the dual if the LP is unfeasible. Else we extract the base from the vector given
*/

/* PROGRAM FLOW *
 *
 * 1- Format LP to SEF
 * 2- Add operation register submatrix
 * 3- Create Auxiliar LP
 * 4- Solve Auxiliar LP
 * 5- If Auxiliar LP optimal value is:
 *    5.1- < 0, the LP is unfeasible, so we should save the optimal solution for the dual of the 
 *         auxiliar LP(vector in [0->1][0->m - 2]) which is the unfeasibility certificate
 *    5.2- = 0, the LP is feasible and so we need to extract the starting base from the resulting vector
 * 6- Now that we have a base to start with, we run the respective simplex. If mode is 1, we run the primal simplex. 
 *    If mode is 2, we run the chosen simplex. 
 *    6.1- Not regarding the chosen mode, we need to check in the simplex if we cannot proceed due to a collumn being 
           all negative in the primal simplex or a row being all positive in the dual simplex, the LP is unbounded. Using 
           previously created mechanisms you should be able to create an algorithm that finds the certificate.
   7- If we get to the end of the simplexes(canonical form, c and b >= 0) we can show the optimal solution and its value
      which is the first element of the last collumn.
      7.1- If mode equals 1 we need to print the optimal solution for the dual as well
      7.2- If mode equals 2 we need to print the final form of the tableau
 */

// TODO: Function that adds matrix to already created ones(So we can substitute format_sef and add_operations_register for only one function)

double** allocate_matrix(int m, int n);
double** identity(int m);
void print_matrix(double** matrix, int m, int n);
void copy_matrix(double** original, double** copy, int m, int n);
void insert_matrix(double** source, double** target, int from_row, int to_row, int from_collumn, int to_collumn);
void operate_on_rows(double** matrix, int m, int n, int multiply_b, int sum_to);
double** format_sef(double** matrix, int m, int n);
double** add_operations_register(double** matrix, int m, int n);
double** create_auxiliar_lp(double** matrix, int m, int n);
void pivoting(double** matrix, int m, int n, int row, int collumn);
void format_canonical(double** matrix, int m, int n, int* base);
void primal_simplex(double** matrix, int m, int n);
void dual_simplex(double** matrix, int m, int n);

int main(int argc, char* argv[]) {
    double** matrix; // LP
    double** auxiliar_lp; 
    double** tmp;
    int i, j, m, n, auxiliar_n, input_size, mode; // Matrix dimensions  
    char simplex_type; // Primal or Dual Simplex

    FILE* input; // Input file
    if(argc >= 2) { // Input file name has been given
		 input = fopen(argv[1], "r+");
    }
    else { // Default input file name
    	input = fopen("input.txt", "r+");
    }

    FILE* output = fopen("output.txt", "w+"); // Output file

    fscanf(input, "%d %d ", &m, &n);

    input_size = (2 + (4 * m) + (m * (3 * n))); // Approximated size of input
    char* input_matrix = malloc(input_size * sizeof(char)); // Unformated LP
    fgets(input_matrix, input_size, input); 

    char* row = strtok(input_matrix, "{"); // Parsed first row

    // Allocate memory for matrix of dimensions m x n
    matrix = allocate_matrix(m, n);

    i = 0;
    j = 0;

    // Parsing of rows to get elements and fill matrix
    while(row != NULL) { 
    	char* element = row;
		while (*element) {
			if ((element[0] == '-' && isdigit(element[1])) || isdigit(element[0])) { // Check if char is a number
		    	double val = strtol(element, &element, 10); // Convert char to number in base 10
		    	matrix[i][j] = val; // Add element to the corresponding position in the matrix
		    	j++;
		    	printf("%g ", val); 
			} 
			else {
		    	element++;
			}
		}
		printf("\n");

    	row = strtok(NULL, "{");
    	i++;
    	j = 0;
    }

    // User chooses the modus operandi
    printf("Escolha o modo de saída: ");
    //scanf(" %d", &mode);
    mode = 1;
    printf("1\n");

    // Adds the slack variables for the problem, by formating it to the standard equalities form
    tmp = matrix;
    matrix = format_sef(matrix, m, n);
    free(tmp);
    n += m - 1;

    printf("Standard Equality Form:\n");
    print_matrix(matrix, m, n);

   	// Now add the register of operations
   	tmp = matrix; // Holds the position in the memory pointed by auxiliar_lp so we can free it later
	matrix = add_operations_register(matrix, m, n);
	free(tmp); // Finally frees it
	n += m - 1; // Sums the new rows added in the previous operation to the total number of rows

	printf("Operations Register:\n");
    print_matrix(matrix, m, n);

    switch(mode) {
	    case 1:
	    auxiliar_lp = create_auxiliar_lp(matrix, m, n);

	    auxiliar_n = n + m - 1; // Value of n for the auxiliar PL with the operation register matrix on its side



	    break;

	    case 2:
	    	printf("Você gostaria de resolver pelo simplex (P)rimal ou (D)ual? ");
	    	scanf(" %c", &simplex_type);

	    	switch(simplex_type) {
	    		case 'P':
	    		printf("P.");
	    		break;

	    		case 'D':
	    		printf("D.");
	    		break;

	    		default:
				printf("Erro: Opção Inválida.\n");
	    	}
	    break;

		default:
		printf("Erro: Opção Inválida.\n");
	}

    for(i = 0; i < m; i++) { // Prints output matrix in file
		for(j = 0; j < n; j++) {
		    fprintf(output, "%g ", matrix[i][j]);
		}
		fprintf(output, "\n");
	}

	fclose(input);
	fclose(output);

/*	
	free(you);
	free(tu);
	free(a porra toda);
*/
    return 0;
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
			printf("%g ", matrix[i][j]);
		}
		printf("\n");
	}
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
void insert_matrix(double** source, double** target, int from_row, int to_row, int from_collumn, int to_collumn) {
	int i, j, m, n;

	m = 0;
	n = 0;

	for(i = (from_row - 1); i < to_row; i++) {
		for(j = (from_collumn - 1); j < to_collumn; j++) {
//			printf("Source: [%d][%d] Target: [%d][%d] Value: %f\n", m, n, i, j, source[m][n]);
			target[i][j] = source[m][n];
			n++;
		}
		m++;
		n = 0;
	}
//	printf("END\n");
}

/* Offers a way to make linear operations. If sum_to is NULL, replaces the actual line. m stands for the row 
 * to operate on and n stands for the dimension of collumns. The index for sum_to begins at 0 */
void operate_on_rows(double** matrix, int m, int n, int multiply_by, int sum_to) {
	int* new_row;
	int j;

	new_row = malloc(n * sizeof(int)); 

	for(j = 0; j < n; j++) {
		new_row[j] = matrix[m][j] * multiply_by;
	}

	if(sum_to == (int)NULL) {
		printf("É nulo\n");
		for(j = 0; j < n; j++) {
			matrix[m][j] = new_row[j]; // TODO - Check if matrix[m] = new_row, or its variants, works
			printf("%d ", new_row[j]);
		}	
		printf("\n");
	}
	else {
		for(j = 0; j < n; j++) {
			matrix[sum_to][j] = new_row[j]; // TODO - Check if matrix[m] = new_row, or its variants, works
		}	
	}
}

// Format the LP to the Standard Equality Form
double** format_sef(double** original_matrix, int m, int n) {
	double** new_matrix;
	int i, j, new_collumns, new_n;
	
	new_collumns = m - 1;
	new_n = n + new_collumns;

	new_matrix = allocate_matrix(m, (n + new_collumns));

	// Set the values in the submatrix that will register the operations
	for(j = 0; j < new_collumns; j++) {
		new_matrix[0][j] = 0; // Set c row to 0
	}

	// Insert the original matrix in the beginning of the new one without the last collumn
	insert_matrix(original_matrix, new_matrix, 1, m, 1, (n - 1));

	// Adds the identity matrix in the correct position
	insert_matrix(identity(new_collumns), new_matrix, 2, m, n, (n + new_collumns)); 

	// Adds the last collumn of the original LP to the new one
	for(i = 0; i < m; i++) {
		new_matrix[i][new_n - 1] = original_matrix[i][n - 1];
	}

	return new_matrix;
}



// Create a new matrix that consists of the operation register submatrix added to the original one.
double** add_operations_register(double** original_matrix, int m, int n) {
	double** new_matrix;
	int i, j, new_collumns;

	new_collumns = m - 1; // Number of rows added

	new_matrix = allocate_matrix(m, (n + new_collumns));

	// Set the values in the submatrix that will register the operations
	for(j = 0; j < new_collumns; j++) {
		new_matrix[0][j] = 0; // Set c row to 0
	}
	insert_matrix(identity(new_collumns), new_matrix, 2, m, 1, new_collumns); // Adds the identity matrix in the correct position

	// Insert the original matrix as a submatrix of the new one
	insert_matrix(original_matrix, new_matrix, 1, m, m, (n + new_collumns));

	return new_matrix;
}

double** create_auxiliar_lp(double** matrix, int m, int n) {
	double** auxiliar_lp;
	double** copied_matrix;
	int i, j, auxiliar_n;

	auxiliar_n = n + m - 1;
	auxiliar_lp = allocate_matrix(m, auxiliar_n);

	copied_matrix = allocate_matrix(m, n); 
	copy_matrix(matrix, copied_matrix, m, n); // Holds the values of the original matrix but doesn't mess with the original data in any sense

	// First check if b > 0
	for(i = 1; i < m; i++) {
		if(copied_matrix[i][(n - 1)] < 0) {
			operate_on_rows(copied_matrix, m, n, -1, (int)NULL);
		}
	}
	
	// Fulfill the auxiliar lp without the last collumn(thats why to_collumn equals n - 2)
	insert_matrix(copied_matrix, auxiliar_lp, 1, m, 1, (n - 1));	

	// Adds the last collumn to the auxiliar lp
	for(i = 0; i < m; i++) {
		auxiliar_lp[i][auxiliar_n - 1] = copied_matrix[i][n - 1];
	}

	// Creates the first row of the auxiliar LP in the correct form
	for(j = 0; j < auxiliar_n; j++) { 
		if(j >= (auxiliar_n - m) && j < (auxiliar_n - 1)) {
			auxiliar_lp[0][j] = -1;
		} 
		else {
			auxiliar_lp[0][j] = 0; 
		}
	} 

	insert_matrix(identity(m - 1), auxiliar_lp, 2, m, (auxiliar_n - m + 1), (auxiliar_n - 1)); // Adds the identity matrix in the correct position

	printf("Auxiliar LP:\n");
	for(i = 0; i < m; i++) {
		for(j = 0; j < auxiliar_n; j++) {
			printf("%g ", round(auxiliar_lp[i][j] * 100000) / 100000);
		}
		printf("\n");
	}

	return auxiliar_lp;
}

// Format the LP to the Canonical Form
void format_canonical(double** matrix, int m, int n, int* base) {

}

// If we cant generalize this procedure for auxiliar and normal LP's, we should separate this into two different functions
void primal_simplex(double** matrix, int m, int n) {	
	int* base; // Current base used in simplex
	int i;

	base = malloc((m - 1) * sizeof(int));

	for(i = 0; i < (m - 1); i++) {
		base[i] = n - m - 1 - 1;
	}

	// First we need to present the LP in the canonical form
	format_canonical(matrix, m, n, base);

}

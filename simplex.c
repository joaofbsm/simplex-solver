#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <ctype.h>

/*
CERTIFICATES

 - INVALID: Solution of the auxiliar PL's dual(Use add_operations_register in the auxiliar PL 
   and the corresponding c collumns are going to be the solution and the certificate as well)

 - UNBOUNDED: Find vector z such that: (HOW CAN WE FIND THIS ONLY LOOKING AT THE TABLEAU? MAYBE SOMETHING WITH THE OPERATIONS REGISTER?)
	
		 - z >= 0	
         - Az = 0
         - cTz > 0

 - OPTIMAL: Solution of the dual -> it's in the first row and the collumns correspondent to the operation_register
*/

double** allocate_matrix(int m, int n);
double** identity(int m);
void print_matrix(double** matrix, int m, int n);
void copy_matrix(double **original, double **copy, int m, int n);
void insert_matrix(double **source, double **target, int from_row, int to_row, int from_collumn, int to_collumn);
void operate_on_rows(double **matrix, int m, int n, int multiply_b, int sum_to);
void pivoting(double **matrix, int m, int n, int row, int collumn);
double** add_operations_register(double **matrix, int m, int n);
double** create_auxiliar_lp(double **matrix, int m, int n);
void solve_tableau(double **matrix, int m, int n, char simplex_type);
void print_tableau(double **matrix, int m, int n);

int main(int argc, char* argv[]) {
    double** matrix; // LP
    double** auxiliar_lp; 
    int* base; // Current base used in simplex
    int i, j, m, n, input_size, mode; // Matrix dimensions  
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

    base = malloc(m * sizeof(int));

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

    switch(mode) {
	    case 1:
	    auxiliar_lp = create_auxiliar_lp(matrix, m, n);

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

// Copy the value of every element in the original matrix to the new one.
void copy_matrix(double **original, double **copy, int m, int n) {
	int i, j;

	for(i = 0; i < m; i++) {
		for(j = 0; j < n; j++) {
			copy[i][j] = original[i][j];
		}
	}
}

/* Insert source matrix inside target matrix in the submatrix comprehended by the integer offsets(from_row, to_row, etc). 
 * The dimensions must match. The indexes start at 0. */
void insert_matrix(double **source, double **target, int from_row, int to_row, int from_collumn, int to_collumn) {
	int i, j, m, n;

	m = 0;
	n = 0;

	for(i = from_row; i <= to_row; i++) {
		for(j = from_collumn; j <= to_collumn; j++) {
//			printf("Source: [%d][%d] Target: [%d][%d] Value: %f\n", m, n, i, j, source[m][n]);
			target[i][j] = source[m][n];
			n++;
		}
		m++;
		n = 0;
	}
//	printf("END\n");
}

/* Offers a way to make linear operations. If sum_to is NULL, replaces the actual line. 
 * m stands for the row to operate on and n stands for the dimension of collumns */
void operate_on_rows(double **matrix, int m, int n, int multiply_by, int sum_to) {
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


// Create a new matrix that consists of the operation register submatrix added to the original one.
double** add_operations_register(double **original_matrix, int m, int n) {
	double** new_matrix;
	int i, j, new_rows;

	new_rows = m - 1; // Number of rows added

	new_matrix = allocate_matrix(m, (n + new_rows));

	// Set the values in the submatrix that will register the operations
	for(j = 0; j < new_rows; j++) {
		new_matrix[0][j] = 0; // Set c row to 0
	}
	insert_matrix(identity(new_rows), new_matrix, 1, (m - 1), 0, (new_rows - 1)); // Adds the identity matrix in the correct position

	// Insert the original matrix as a submatrix of the new one
	insert_matrix(original_matrix, new_matrix, 0, (m - 1), (m - 1), (n + new_rows - 1));

	return new_matrix;
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

double** create_auxiliar_lp(double **matrix, int m, int n) {
	double** auxiliar_lp;
	double** copied_matrix;
	double** tmp;
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
	insert_matrix(copied_matrix, auxiliar_lp, 0, (m - 1), 0, (n - 2));	

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

	// ESSA LINHA ESTA ADICIONANDO A IDENTIDADE NO INICIO!!!
	insert_matrix(identity(m - 1), auxiliar_lp, 1, (m - 1), (auxiliar_n - m), (auxiliar_n - 2)); // Adds the identity matrix in the correct position

	tmp = auxiliar_lp; // Holds the position in the memory pointed by auxiliar_lp so we can free it later

	// Now add the register of operations
	auxiliar_lp = add_operations_register(auxiliar_lp, m, auxiliar_n);
	auxiliar_n += m - 1; // Sums the new rows added in the previous operation to the total number of rows
	
	free(tmp); // Finally frees it

	printf("Auxiliar LP:\n");
	for(i = 0; i < m; i++) {
		for(j = 0; j < auxiliar_n; j++) {
			printf("%g ", round(auxiliar_lp[i][j] * 100000) / 100000);
		}
		printf("\n");
	}

	return auxiliar_lp;
}

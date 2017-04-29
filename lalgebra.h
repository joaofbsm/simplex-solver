/* Linear Algebra Library
 * Developed by Joao Francisco B. S. Martins <joaofbsm@dcc.ufmg.br>
 */

#ifndef __LALGEBRA_HEADER__
#define __LALGEBRA_HEADER__

void parse_input(FILE* input, double** matrix, int m, int n);

double** allocate_matrix(int m, int n);
double** identity(int m);
void print_matrix(double** matrix, int m, int n);
void copy_matrix(double** original, double** copy, int m, int n);
void insert_matrix(double** source, double** target, int from_row, int to_row, int from_column, int to_column);
void operate_on_rows(double** matrix, int m, int n, double multiply_by, int sum_to);
void operate_on_columns(double** matrix, int m, int n, double multiply_by, int sum_to);

double** format_sef(double** matrix, int m, int n);
void format_tableau(double** matrix, int m, int n);
double** add_operations_register(double** matrix, int m, int n);
int is_b_negative(double** matrix, int m, int n);
void make_b_non_negative(double** matrix, int m, int n);
int is_c_negative(double** matrix, int m, int n);
void make_c_non_negative(double** matrix, int m, int n);
double** create_auxiliar_lp(double** matrix, int m, int n);
void set_initial_base(double** matrix, int m, int n, int* base);
int find_non_zero_element(double** matrix, int m, int column);
void format_canonical(double** matrix, int m, int n, int* base);
int primal_next_base(double** matrix, int m, int n, int* base_row, int* base_column);
int primal_simplex(double** matrix, int m, int n, int* base, int print_output);
int dual_next_base(double** matrix, int m, int n, int* base_row, int* base_column);
int dual_simplex(double** matrix, int m, int n, int* base, int print_output);

double* get_primal_optimal_solution(double** matrix, int m, int n, int* base);
double* get_dual_optimal_solution(double** matrix, int m);
double* generate_unboundedness_certificate(double** matrix, int m, int n, int column, int* base);

void print_output_vector(double* vector, int n);
void print_output_matrix(double** matrix, int m, int n);

#endif
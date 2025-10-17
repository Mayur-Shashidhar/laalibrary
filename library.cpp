#include <iostream>
#include<math.h>
using namespace std;

class Matrix
{
	public:
		//Data Members
		int r; //number of rows
		int c; //number of columns
		float m[100][100]; //matrix

		//Constructor
		Matrix(int row, int column, bool fill=true);
		/*
		 * boolean flag fill defaults to true to prevent this change from
		 * being a breaking change. If you do not want to be asked to input a
		 * value for every single element in the matrix, set fill to false.
		 */

		//Utility functions
		void printMatrix();

		//Member functions
		Matrix addition(); //Matrix addition
		Matrix subtraction(); //Matrix subtraction
		Matrix multiplication(); //Matrix multiplication
		int isIdentity(); //Returns 1 if it is an identity matrix else returns 0
		int isSquare(); //Checks if matrix is a square matrix
		int trace(); //Find the trace of a given matrix
		int *dimensions(); //Find the dimensions of a matrix
		int *gaussElimination();
		double determinant(int , float [100][100]);
		Matrix columnSpace();
		Matrix transpose();
		Matrix SilentTranspose();
		Matrix nullSpace();
		Matrix inverse();
		int *eigenValues();
		Matrix eigenVectors();
		int *graph(); //If the matrix represents a graph return number of edges and vertices else return NULL
		Matrix *LUdecomposition();
		Matrix *LDUdecomposition();
		Matrix *SVD(); //Single-value decomposition
		int determinant(); //Finds the determinant of the given matrix
		Matrix adjoint(); //Finds the adjoint of the given matrix
		int isInvertible(); //Returns 1 if the matrix is invertible else returns 0
		Matrix scalarMul(); //Multiplies the matrix by a scalar value
		int isIdempotent(); //Returns 1 if the matrix is idempotent else returns 0
		int isInvolutory(); //Returns 1 if the matrix is involutory else returns 0
		int isNilpotent(); //Returns 1 if the matrix is nilpotent else returns 0
		int duplicate(); //returns the number of duplicate numbers in the matrix
		Matrix additiveInv(); //finds the additive inverse of the matrix
		Matrix cofactor(float [100][100],int ,int ,int );

		// since this function computes two matrices, we return a pointer to a matrix array
		// TODO check memory safety of this
		Matrix* symmskew(); //expresses the matrix as a sum of a symmetric and skew symmetric matrix

		// operator overloading for addition and subtraction of like matrices
		Matrix operator + (const Matrix& other)
		{
			if (r == other.r && c == other.c)
			{
				// we can perform operations on them
				Matrix result(r, c, false);

				for (int i = 0; i < r; i++)
				{
					for (int j = 0; j < c; j++)
					{
						result.m[i][j] = m[i][j] + other.m[i][j];
					}
				}
				cout << "Added Matrices successfully." << endl;
				return result;
			}
			else
			{
				cout << "Cannot add two unlike matrices" << endl;
				return Matrix(1,1, false);
			}
		};

		Matrix operator - (const Matrix& other)
		{
			if (r == other.r && c == other.c)
			{
				Matrix result(r, c, false);

				for (int i = 0; i < r; i++)
				{
					for (int j = 0; j < c; j++)
					{
						result.m[i][j] = m[i][j] - other.m[i][j];
					}
				}
				cout << "Subtracted Matrices successfully." << endl;
				return result;
			}
			else
			{
				cout << "Cannot subtract two unlike matrices" << endl;
				return Matrix(1,1, false);
			}
		};

		// overload scalar multiplication
		Matrix operator * (const int& scalar)
		{
			Matrix result(r, c, false);

			for (int i = 0; i < r; i++)
			{
				for (int j = 0; j < c; j++)
				{
					result.m[i][j] = m[i][j] * scalar;
				}
			}
			cout << "Successfully multiplied scalar" << endl;
			return result;
		};

		// overload scalar division
		Matrix operator / (const int& scalar)
		{
			Matrix result(r, c, false);

			for (int i = 0; i < r; i++)
			{
				for (int j = 0; j < c; j++)
				{
					result.m[i][j] = m[i][j] / scalar;
				}
			}
			cout << "Subtracted Matrices successfully." << endl;
			return result;
		};
};

Matrix::Matrix(int row, int column, bool fill)
{
	r=row;
	c=column;
	if (fill == true)
	{
		for(int i=0;i<row;i++)
			for(int j=0;j<column;j++)
			{
				float x;
				cin>>x;
				m[i][j]=x/1.0;
			}
	}
	else
	{
		for(int i=0;i<row;i++)
		{
			for(int j=0;j<column;j++)
			{
				m[i][j]=0;
			}
		}
	}
}

void Matrix::printMatrix()
{
	for(int i=0;i<r;i++)
		for(int j=0;j<c;j++)
			// Prints ' ' if j != n-1 else prints '\n'
          		cout << m[i][j] << " \n"[j == c-1];
}

Matrix Matrix::addition() {
    Matrix m1(r,c);
    Matrix m2(r,c);
    Matrix res(r,c);
    for(int i=0;i<r;i++) {
        for(int j=0;j<c;j++) {
            res.m[i][j] = m1.m[i][j] + m2.m[i][j];
        }
    }
   return res;
}

Matrix Matrix::subtraction() { //PR for subtraction
    Matrix m1(r,c);
    Matrix m2(r,c);
    Matrix res(r,c);
    for(int i=0;i<r;i++) {
        for(int j=0;j<c;j++) {
            res.m[i][j] = m1.m[i][j] - m2.m[i][j];
        }
    }
   return res;
}

Matrix Matrix::multiplication() {
    Matrix m1(r,c);
    Matrix m2(r,c);
    Matrix res(r,c);
    int i, j, k;
        for (i = 0; i < r; i++) {
            for (j = 0; j < r; j++) {
                res.m[i][j] = 0;
                for (k = 0; k < r; k++)
                    res.m[i][j] += m1.m[i][k]
                                 * m2.m[k][j];
            }
        }
   return res;
}

int Matrix::isIdentity() { //pr for isidentity
    int flag=0;
    Matrix m1(r,c);
    for(int i=0;i<r;i++)
       {
           if(m1.m[i][i] == 1)
            flag ++;;
       }
       if(flag == m1.r)
       return 1;
       else
        return 0;
}

int Matrix::trace() {
    Matrix m1(r,c);
    int sum=0;
    for(int i=0;i<r;i++)
        sum  = sum + m1.m[i][i];
    return sum;
}

int Matrix::isSquare() {
    Matrix m1(r,c);
    if(r==c)
        return 1;
    else
        return 0;
}

int *Matrix::dimensions() {
    Matrix m1(r,c);
    int arr[2];
    arr[0]=r;
    arr[1]=c;
    return arr;

}

int *Matrix::gaussElimination()
{
	int x_arr[100];
	int *ptr=x_arr;
	for(int i=0;i<c;i++)
		for(int j=i+1;j<r;j++)
		{
			float x=m[j][i]/m[i][i];
			*ptr=x;
			ptr++;
			for(int k=0;k<c;k++)
				m[j][k]=m[j][k]-(x*m[i][k]);
			printMatrix();
			cout<<"\n";
		}
	return x_arr;
}

Matrix Matrix::columnSpace()
{
	gaussElimination();
	int k;
	for(int i=0;i<r;i++)
		if(m[i][i]==0)
		{
			k=i;
			break;
		}
	Matrix out(r,k);
	for(int i=0;i<r;i++)
		for(int j=0;j<k;j++)
		out.m[i][j]=m[i][j];
	return out;
}

Matrix Matrix::transpose()
{
	Matrix out(c,r);
	for(int i=0;i<r;i++)
		for(int j=0;j<c;j++)
			out.m[j][i]=m[i][j];
	return out;
}

Matrix Matrix::SilentTranspose()
{
	Matrix out(c, r, false);
	for(int i=0;i<r;i++)
		for(int j=0;j<c;j++)
			out.m[j][i]=m[i][j];
	return out;
}

double Matrix::determinant(int n,float mat[100][100])
{
    double D = 0;

    //  Base case : if matrix contains single element
    if (n == 1)
        return mat[0][0];

     Matrix temp(10,10); // To store cofactors

    int sign = 1;  // To store sign multiplier


    for (int f = 0; f < n; f++)
    {
        // Getting Cofactor of A[0][f]
        temp = cofactor(mat, 0, f, n);
        D = D + sign * mat[0][f] * determinant(n-1,temp.m);

        // terms are to be added with alternate sign
        sign = -sign;
    }

    return D;
}

Matrix *Matrix::LUdecomposition(){
	 Matrix upper(r,c);
         Matrix lower(r,c);
	 Matrix mat(r,c);
    // Decomposing matrix into Upper and Lower triangular matrix
    for (int i = 0; i < r; i++) {

	// Upper Triangular
        for (int k = i; k < c; k++) {
            //summation of lower[i][j]*upper[j][k]
            int sum = 0;
            for (int j = 0; j < i; j++)
                sum += (lower.m[i][j] * upper.m[j][k]);
            upper.m[i][k] = mat.m[i][k] - sum;
        }

        // Lower Triangular
        for (int k = i; k < c; k++) {
            if (i == k)
                lower.m[i][i] = 1;   // Diagonal as 1
            else {
		//summation of lower[k][j]*upper[j][i]
                int sum = 0;
                for (int j = 0; j < i; j++)
                    sum += (lower.m[k][j] * upper.m[j][i]);
		// Evaluating lower matrix
                lower.m[k][i] = (mat.m[k][i] - sum) / upper.m[i][i];
            }
        }
    }
	upper.printMatrix();
	lower.printMatrix();
	return 0;
}


int Matrix::isIdempotent(){
	Matrix mat(r,c);
	Matrix result(r,c);
	int i, j, k,f=0;
        for (i = 0; i < r; i++) {
            for (j = 0; j < r; j++) {
                result.m[i][j] = 0;
                for (k = 0; k < r; k++)
                    result.m[i][j] += mat.m[i][k]
                                 * mat.m[k][j];
            }
        }
	for (i = 0; i < r; i++) {
            for (j = 0; j < r;) {
		    if(result.m[i][j] == mat.m[i][j])
			    j++;
		    else{
			    f = 1;
			    break;
		    }
	    }
	}
	if(f==0)
		return 1;  //indicates that the matrix is idempotent
	else
		return 0;  //indicates that the matrix is not idempotent
}

int Matrix::isInvolutory(){
	Matrix mat(r,c);
	Matrix result(r,c);
	int i, j, k,f=0;
        for (i = 0; i < r; i++) {
            for (j = 0; j < r; j++) {
                result.m[i][j] = 0;
                for (k = 0; k < r; k++)
                    result.m[i][j] += mat.m[i][k]
                                 * mat.m[k][j];
            }
        }
	for (i = 0; i < r; i++) {
            for (j = 0; j < r;j++) {
		    if(i == j){
			    if(result.m[i][j] == 1)
				    continue;
			    else{
				    f=1;
				    break;
			    }
		    }
		    else if(result.m[i][j]==0)
		            continue;
	            else{
			    f = 1;
			    break;
		    }
	    }
	}
	if(f==0)
		return 1;  //indicates that the matrix is involutory
	else
		return 0;  //indicates that the matrix is not involutory
}


Matrix Matrix::additiveInv() {
    Matrix m1(c,r);
    Matrix out(c,r);
    for(int i=0;i<m1.r;i++)
    {
        for(int j=0;j<m1.c;j++)
        {
            out.m[i][j] = (-1)*m1.m[i][j];
        }
    }
    return out;
}

Matrix* Matrix::symmskew()
{
	if (r == c)
	{
		// get copies of both matrices needed
		Matrix transpose = SilentTranspose();
		Matrix current = (*this);

		// make empty placeholder matrices
		Matrix symmetric(r, c, false);
		Matrix skew_symmetric(r, c, false);

		// perform operations
		symmetric = (current + transpose)/2;
		skew_symmetric = (current - transpose)/2;

		// initialize pointer to array of matrices
		Matrix* result_pointer = NULL;

		// assign positions
		result_pointer[0] = symmetric;
		result_pointer[1] = skew_symmetric;

		// return
		return result_pointer;
	}
	else
	{
		cout << "Matrix must be square \n";
		return NULL;
	}
}


//takes in matrix , row no,column no, and size of matrix
Matrix Matrix::cofactor(float a[100][100],int rw,int cl,int N)
{
    int i = 0, j = 0;
    Matrix temp(c,r);

    for (int row = 0; row < N; row++)
    {
        for (int col = 0; col < N; col++)
        {
            //  Copying into temporary matrix only those element
            //  which are not in given row and column
            if (row != rw && col != cl)
            {
                temp.m[i][j++] = a[row][col];


                if (j == N - 1)
                {
                    j = 0;
                    i++;
                }
            }
        }
    }
    return temp;
}

/*
int Matrix::isInvertible()
{
    Matrix m2(r,c);
    float m3[10][10];
    m3[10][10]= m2.m[10][10];
    float d1=determinant(r,m3);
    if(isSquare() == 1)
    {
        if(d1 != 0)
            return 1;
    }
    else
        return 0;
}
*/


void displayMenu() {
    cout << "\n========== Linear Algebra Library ==========\n";
    cout << "1.  Create and Display Matrix\n";
    cout << "2.  Matrix Addition\n";
    cout << "3.  Matrix Subtraction\n";
    cout << "4.  Matrix Multiplication\n";
    cout << "5.  Scalar Multiplication\n";
    cout << "6.  Matrix Transpose\n";
    cout << "7.  Calculate Determinant\n";
    cout << "8.  Check if Identity Matrix\n";
    cout << "9.  Check if Square Matrix\n";
    cout << "10. Calculate Trace\n";
    cout << "11. Gauss Elimination\n";
    cout << "12. Column Space\n";
    cout << "13. LU Decomposition\n";
    cout << "14. Check if Idempotent\n";
    cout << "15. Check if Involutory\n";
    cout << "16. Additive Inverse\n";
    cout << "17. Symmetric/Skew-Symmetric Decomposition\n";
    cout << "0.  Exit\n";
    cout << "============================================\n";
    cout << "Enter your choice: ";
}

Matrix createMatrix() {
    int rows, cols;
    cout << "Enter number of rows: ";
    cin >> rows;
    cout << "Enter number of columns: ";
    cin >> cols;
    
    if (rows <= 0 || cols <= 0 || rows > 100 || cols > 100) {
        cout << "Invalid matrix dimensions! Please enter values between 1 and 100.\n";
        cout << "Enter number of rows: ";
        cin >> rows;
        cout << "Enter number of columns: ";
        cin >> cols;
        // Clamp values to valid range
        if (rows <= 0) rows = 1;
        if (cols <= 0) cols = 1;
        if (rows > 100) rows = 100;
        if (cols > 100) cols = 100;
    }
    
    cout << "Enter matrix elements:\n";
    Matrix m(rows, cols, true);
    return m;
}

Matrix createMatrixSilent(int rows, int cols) {
    if (rows <= 0 || cols <= 0 || rows > 100 || cols > 100) {
        cout << "Invalid matrix dimensions! Using default 3x3 matrix instead.\n";
        rows = cols = 3;
    }
    
    cout << "Enter matrix elements:\n";
    Matrix m(rows, cols, true);
    return m;
}

// driver code
int main()
{
    int choice;
    Matrix* matrix1 = nullptr;
    Matrix* matrix2 = nullptr;
    
    cout << "Welcome to Linear Algebra Library!\n";
    
    while (true) {
        displayMenu();
        cin >> choice;
        
        switch (choice) {
            case 0:
                cout << "Thank you for using Linear Algebra Library!\n";
                if (matrix1) delete matrix1;
                if (matrix2) delete matrix2;
                return 0;
                
            case 1: {
                cout << "\n--- Create and Display Matrix ---\n";
                if (matrix1) delete matrix1;
                matrix1 = new Matrix(createMatrix());
                cout << "\nMatrix created successfully:\n";
                matrix1->printMatrix();
                break;
            }
            
            case 2: {
                cout << "\n--- Matrix Addition ---\n";
                cout << "Enter first matrix:\n";
                Matrix m1 = createMatrix();
                cout << "Enter second matrix (same dimensions):\n";
                Matrix m2 = createMatrixSilent(m1.r, m1.c);
                
                if (m1.r == m2.r && m1.c == m2.c) {
                    Matrix result = m1 + m2;
                    cout << "\nResult of addition:\n";
                    result.printMatrix();
                } else {
                    cout << "Error: Matrix dimensions don't match for addition!\n";
                }
                break;
            }
            
            case 3: {
                cout << "\n--- Matrix Subtraction ---\n";
                cout << "Enter first matrix:\n";
                Matrix m1 = createMatrix();
                cout << "Enter second matrix (same dimensions):\n";
                Matrix m2 = createMatrixSilent(m1.r, m1.c);
                
                if (m1.r == m2.r && m1.c == m2.c) {
                    Matrix result = m1 - m2;
                    cout << "\nResult of subtraction:\n";
                    result.printMatrix();
                } else {
                    cout << "Error: Matrix dimensions don't match for subtraction!\n";
                }
                break;
            }
            
            case 4: {
                cout << "\n--- Matrix Multiplication ---\n";
                cout << "Enter first matrix:\n";
                Matrix m1 = createMatrix();
                cout << "Enter second matrix:\n";
                Matrix m2 = createMatrix();
                
                if (m1.c == m2.r) {
                    Matrix result = m1.multiplication();
                    cout << "\nResult of multiplication:\n";
                    result.printMatrix();
                } else {
                    cout << "Error: Cannot multiply matrices! First matrix columns must equal second matrix rows.\n";
                }
                break;
            }
            
            case 5: {
                cout << "\n--- Scalar Multiplication ---\n";
                if (!matrix1) {
                    cout << "Please create a matrix first (option 1).\n";
                    break;
                }
                int scalar;
                cout << "Enter scalar value: ";
                cin >> scalar;
                
                Matrix result = (*matrix1) * scalar;
                cout << "\nResult of scalar multiplication:\n";
                result.printMatrix();
                break;
            }
            
            case 6: {
                cout << "\n--- Matrix Transpose ---\n";
                if (!matrix1) {
                    cout << "Please create a matrix first (option 1).\n";
                    break;
                }
                
                Matrix result = matrix1->transpose();
                cout << "\nTranspose of the matrix:\n";
                result.printMatrix();
                break;
            }
            
            case 7: {
                cout << "\n--- Calculate Determinant ---\n";
                if (!matrix1) {
                    cout << "Please create a matrix first (option 1).\n";
                    break;
                }
                
                if (matrix1->isSquare()) {
                    double det = matrix1->determinant(matrix1->r, matrix1->m);
                    cout << "Determinant: " << det << endl;
                } else {
                    cout << "Error: Determinant can only be calculated for square matrices!\n";
                }
                break;
            }
            
            case 8: {
                cout << "\n--- Check if Identity Matrix ---\n";
                if (!matrix1) {
                    cout << "Please create a matrix first (option 1).\n";
                    break;
                }
                
                if (matrix1->isIdentity()) {
                    cout << "The matrix is an identity matrix.\n";
                } else {
                    cout << "The matrix is not an identity matrix.\n";
                }
                break;
            }
            
            case 9: {
                cout << "\n--- Check if Square Matrix ---\n";
                if (!matrix1) {
                    cout << "Please create a matrix first (option 1).\n";
                    break;
                }
                
                if (matrix1->isSquare()) {
                    cout << "The matrix is a square matrix.\n";
                } else {
                    cout << "The matrix is not a square matrix.\n";
                }
                break;
            }
            
            case 10: {
                cout << "\n--- Calculate Trace ---\n";
                if (!matrix1) {
                    cout << "Please create a matrix first (option 1).\n";
                    break;
                }
                
                if (matrix1->isSquare()) {
                    int tr = matrix1->trace();
                    cout << "Trace of the matrix: " << tr << endl;
                } else {
                    cout << "Error: Trace can only be calculated for square matrices!\n";
                }
                break;
            }
            
            case 11: {
                cout << "\n--- Gauss Elimination ---\n";
                if (!matrix1) {
                    cout << "Please create a matrix first (option 1).\n";
                    break;
                }
                
                cout << "Performing Gauss Elimination:\n";
                matrix1->gaussElimination();
                cout << "\nMatrix after Gauss Elimination:\n";
                matrix1->printMatrix();
                break;
            }
            
            case 12: {
                cout << "\n--- Column Space ---\n";
                if (!matrix1) {
                    cout << "Please create a matrix first (option 1).\n";
                    break;
                }
                
                Matrix result = matrix1->columnSpace();
                cout << "Column space of the matrix:\n";
                result.printMatrix();
                break;
            }
            
            case 13: {
                cout << "\n--- LU Decomposition ---\n";
                if (!matrix1) {
                    cout << "Please create a matrix first (option 1).\n";
                    break;
                }
                
                if (matrix1->isSquare()) {
                    cout << "Performing LU Decomposition:\n";
                    matrix1->LUdecomposition();
                } else {
                    cout << "Error: LU decomposition requires a square matrix!\n";
                }
                break;
            }
            
            case 14: {
                cout << "\n--- Check if Idempotent ---\n";
                if (!matrix1) {
                    cout << "Please create a matrix first (option 1).\n";
                    break;
                }
                
                if (matrix1->isSquare()) {
                    if (matrix1->isIdempotent()) {
                        cout << "The matrix is idempotent (A^2 = A).\n";
                    } else {
                        cout << "The matrix is not idempotent.\n";
                    }
                } else {
                    cout << "Error: Idempotent check requires a square matrix!\n";
                }
                break;
            }
            
            case 15: {
                cout << "\n--- Check if Involutory ---\n";
                if (!matrix1) {
                    cout << "Please create a matrix first (option 1).\n";
                    break;
                }
                
                if (matrix1->isSquare()) {
                    if (matrix1->isInvolutory()) {
                        cout << "The matrix is involutory (A^2 = I).\n";
                    } else {
                        cout << "The matrix is not involutory.\n";
                    }
                } else {
                    cout << "Error: Involutory check requires a square matrix!\n";
                }
                break;
            }
            
            case 16: {
                cout << "\n--- Additive Inverse ---\n";
                if (!matrix1) {
                    cout << "Please create a matrix first (option 1).\n";
                    break;
                }
                
                Matrix result = matrix1->additiveInv();
                cout << "Additive inverse of the matrix:\n";
                result.printMatrix();
                break;
            }
            
            case 17: {
                cout << "\n--- Symmetric/Skew-Symmetric Decomposition ---\n";
                if (!matrix1) {
                    cout << "Please create a matrix first (option 1).\n";
                    break;
                }
                
                if (matrix1->isSquare()) {
                    Matrix* result = matrix1->symmskew();
                    if (result != NULL) {
                        cout << "Symmetric part:\n";
                        result[0].printMatrix();
                        cout << "\nSkew-symmetric part:\n";
                        result[1].printMatrix();
                    }
                } else {
                    cout << "Error: Symmetric/Skew-symmetric decomposition requires a square matrix!\n";
                }
                break;
            }
            
            default:
                cout << "Invalid choice! Please enter a number between 0 and 17.\n";
                break;
        }
        
        cout << "\nPress Enter to continue...";
        cin.ignore();
        cin.get();
    }
    
    return 0;
}

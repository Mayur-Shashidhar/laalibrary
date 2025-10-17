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

// Global array to store matrices for menu operations
Matrix* storedMatrices[10];
int matrixCount = 0;

// Function to clear input buffer
void clearInputBuffer() {
    cin.clear();
    cin.ignore(10000, '\n');
}

// Function to get valid integer input
int getValidInput(int min, int max) {
    int choice;
    while (true) {
        cout << "Enter your choice (" << min << "-" << max << "): ";
        if (cin >> choice && choice >= min && choice <= max) {
            clearInputBuffer();
            return choice;
        }
        cout << "Invalid input! Please enter a number between " << min << " and " << max << ".\n";
        clearInputBuffer();
    }
}

// Function to display all stored matrices
void displayStoredMatrices() {
    if (matrixCount == 0) {
        cout << "\nNo matrices stored yet.\n";
        return;
    }
    
    cout << "\n=== Stored Matrices ===\n";
    for (int i = 0; i < matrixCount; i++) {
        cout << "\nMatrix " << (i + 1) << " (" << storedMatrices[i]->r << "x" << storedMatrices[i]->c << "):\n";
        storedMatrices[i]->printMatrix();
    }
}

// Function to select a matrix from stored matrices
int selectMatrix() {
    if (matrixCount == 0) {
        cout << "No matrices available! Please create a matrix first.\n";
        return -1;
    }
    
    displayStoredMatrices();
    cout << "\nSelect matrix (1-" << matrixCount << "): ";
    int choice = getValidInput(1, matrixCount);
    return choice - 1; // Convert to 0-based index
}

// Function to create a new matrix
void createMatrix() {
    if (matrixCount >= 10) {
        cout << "Maximum number of matrices (10) reached!\n";
        return;
    }
    
    int rows, cols;
    cout << "\nEnter number of rows: ";
    cin >> rows;
    cout << "Enter number of columns: ";
    cin >> cols;
    
    if (rows <= 0 || cols <= 0 || rows > 100 || cols > 100) {
        cout << "Invalid dimensions! Rows and columns must be between 1 and 100.\n";
        return;
    }
    
    storedMatrices[matrixCount] = new Matrix(rows, cols);
    cout << "\nMatrix " << (matrixCount + 1) << " created successfully!\n";
    cout << "Matrix contents:\n";
    storedMatrices[matrixCount]->printMatrix();
    matrixCount++;
}

// Function to display matrix properties menu
void matrixPropertiesMenu(int index) {
    Matrix* mat = storedMatrices[index];
    
    while (true) {
        cout << "\n=== Matrix Properties Menu ===\n";
        cout << "1. Check if Square\n";
        cout << "2. Check if Identity\n";
        cout << "3. Calculate Trace\n";
        cout << "4. Get Dimensions\n";
        cout << "5. Check if Idempotent\n";
        cout << "6. Check if Involutory\n";
        cout << "7. Calculate Determinant\n";
        cout << "8. Back to main menu\n";
        
        int choice = getValidInput(1, 8);
        
        switch (choice) {
            case 1:
                cout << "Matrix is " << (mat->isSquare() ? "square" : "not square") << "\n";
                break;
            case 2:
                cout << "Matrix is " << (mat->isIdentity() ? "an identity matrix" : "not an identity matrix") << "\n";
                break;
            case 3:
                if (mat->isSquare()) {
                    cout << "Trace: " << mat->trace() << "\n";
                } else {
                    cout << "Trace can only be calculated for square matrices!\n";
                }
                break;
            case 4: {
                int* dims = mat->dimensions();
                cout << "Dimensions: " << dims[0] << " x " << dims[1] << "\n";
                break;
            }
            case 5:
                if (mat->isSquare()) {
                    cout << "Matrix is " << (mat->isIdempotent() ? "idempotent" : "not idempotent") << "\n";
                } else {
                    cout << "Idempotent check only applies to square matrices!\n";
                }
                break;
            case 6:
                if (mat->isSquare()) {
                    cout << "Matrix is " << (mat->isInvolutory() ? "involutory" : "not involutory") << "\n";
                } else {
                    cout << "Involutory check only applies to square matrices!\n";
                }
                break;
            case 7:
                if (mat->isSquare()) {
                    double det = mat->determinant(mat->r, mat->m);
                    cout << "Determinant: " << det << "\n";
                } else {
                    cout << "Determinant can only be calculated for square matrices!\n";
                }
                break;
            case 8:
                return;
        }
    }
}

// Function to perform matrix operations
void matrixOperationsMenu() {
    while (true) {
        cout << "\n=== Matrix Operations Menu ===\n";
        cout << "1. Matrix Addition (A + B)\n";
        cout << "2. Matrix Subtraction (A - B)\n";
        cout << "3. Matrix Multiplication (A * B)\n";
        cout << "4. Scalar Multiplication\n";
        cout << "5. Scalar Division\n";
        cout << "6. Transpose\n";
        cout << "7. Additive Inverse\n";
        cout << "8. Gaussian Elimination\n";
        cout << "9. Column Space\n";
        cout << "10. LU Decomposition\n";
        cout << "11. Symmetric/Skew-Symmetric Decomposition\n";
        cout << "12. Back to main menu\n";
        
        int choice = getValidInput(1, 12);
        
        switch (choice) {
            case 1: { // Addition
                cout << "\nSelect first matrix:\n";
                int idx1 = selectMatrix();
                if (idx1 == -1) break;
                
                cout << "\nSelect second matrix:\n";
                int idx2 = selectMatrix();
                if (idx2 == -1) break;
                
                Matrix result = (*storedMatrices[idx1]) + (*storedMatrices[idx2]);
                cout << "\nResult of addition:\n";
                result.printMatrix();
                break;
            }
            case 2: { // Subtraction
                cout << "\nSelect first matrix:\n";
                int idx1 = selectMatrix();
                if (idx1 == -1) break;
                
                cout << "\nSelect second matrix:\n";
                int idx2 = selectMatrix();
                if (idx2 == -1) break;
                
                Matrix result = (*storedMatrices[idx1]) - (*storedMatrices[idx2]);
                cout << "\nResult of subtraction:\n";
                result.printMatrix();
                break;
            }
            case 3: { // Multiplication (using existing method)
                cout << "\nSelect matrix for multiplication:\n";
                int idx = selectMatrix();
                if (idx == -1) break;
                
                Matrix result = storedMatrices[idx]->multiplication();
                cout << "\nResult of multiplication:\n";
                result.printMatrix();
                break;
            }
            case 4: { // Scalar multiplication
                cout << "\nSelect matrix:\n";
                int idx = selectMatrix();
                if (idx == -1) break;
                
                int scalar;
                cout << "Enter scalar value: ";
                cin >> scalar;
                
                Matrix result = (*storedMatrices[idx]) * scalar;
                cout << "\nResult of scalar multiplication:\n";
                result.printMatrix();
                break;
            }
            case 5: { // Scalar division
                cout << "\nSelect matrix:\n";
                int idx = selectMatrix();
                if (idx == -1) break;
                
                int scalar;
                cout << "Enter scalar value (non-zero): ";
                cin >> scalar;
                if (scalar == 0) {
                    cout << "Cannot divide by zero!\n";
                    break;
                }
                
                Matrix result = (*storedMatrices[idx]) / scalar;
                cout << "\nResult of scalar division:\n";
                result.printMatrix();
                break;
            }
            case 6: { // Transpose
                cout << "\nSelect matrix:\n";
                int idx = selectMatrix();
                if (idx == -1) break;
                
                Matrix result = storedMatrices[idx]->transpose();
                cout << "\nTranspose:\n";
                result.printMatrix();
                break;
            }
            case 7: { // Additive inverse
                cout << "\nSelect matrix:\n";
                int idx = selectMatrix();
                if (idx == -1) break;
                
                Matrix result = storedMatrices[idx]->additiveInv();
                cout << "\nAdditive Inverse:\n";
                result.printMatrix();
                break;
            }
            case 8: { // Gaussian elimination
                cout << "\nSelect matrix:\n";
                int idx = selectMatrix();
                if (idx == -1) break;
                
                cout << "\nPerforming Gaussian Elimination...\n";
                storedMatrices[idx]->gaussElimination();
                cout << "\nMatrix after Gaussian Elimination:\n";
                storedMatrices[idx]->printMatrix();
                break;
            }
            case 9: { // Column space
                cout << "\nSelect matrix:\n";
                int idx = selectMatrix();
                if (idx == -1) break;
                
                Matrix result = storedMatrices[idx]->columnSpace();
                cout << "\nColumn Space:\n";
                result.printMatrix();
                break;
            }
            case 10: { // LU decomposition
                cout << "\nSelect square matrix:\n";
                int idx = selectMatrix();
                if (idx == -1) break;
                
                if (!storedMatrices[idx]->isSquare()) {
                    cout << "LU decomposition requires a square matrix!\n";
                    break;
                }
                
                cout << "\nPerforming LU Decomposition...\n";
                storedMatrices[idx]->LUdecomposition();
                break;
            }
            case 11: { // Symmetric/Skew-symmetric decomposition
                cout << "\nSelect square matrix:\n";
                int idx = selectMatrix();
                if (idx == -1) break;
                
                if (!storedMatrices[idx]->isSquare()) {
                    cout << "Symmetric/Skew-symmetric decomposition requires a square matrix!\n";
                    break;
                }
                
                cout << "\nPerforming Symmetric/Skew-symmetric Decomposition...\n";
                Matrix* result = storedMatrices[idx]->symmskew();
                if (result != NULL) {
                    cout << "\nSymmetric part:\n";
                    result[0].printMatrix();
                    cout << "\nSkew-symmetric part:\n";
                    result[1].printMatrix();
                }
                break;
            }
            case 12:
                return;
        }
    }
}

// Main menu function
void displayMainMenu() {
    cout << "\n===============================================\n";
    cout << "         LINEAR ALGEBRA LIBRARY MENU          \n";
    cout << "===============================================\n";
    cout << "1. Create New Matrix\n";
    cout << "2. Display All Matrices\n";
    cout << "3. Matrix Properties\n";
    cout << "4. Matrix Operations\n";
    cout << "5. Delete Matrix\n";
    cout << "6. Exit\n";
    cout << "===============================================\n";
}

// driver code
int main()
{
    cout << "Welcome to the Linear Algebra Library!\n";
    cout << "This program allows you to create matrices and perform various operations.\n";
    
    while (true) {
        displayMainMenu();
        int choice = getValidInput(1, 6);
        
        switch (choice) {
            case 1:
                createMatrix();
                break;
                
            case 2:
                displayStoredMatrices();
                break;
                
            case 3: {
                if (matrixCount == 0) {
                    cout << "No matrices available! Please create a matrix first.\n";
                    break;
                }
                cout << "\nSelect matrix to view properties:\n";
                int idx = selectMatrix();
                if (idx != -1) {
                    matrixPropertiesMenu(idx);
                }
                break;
            }
            
            case 4:
                if (matrixCount == 0) {
                    cout << "No matrices available! Please create a matrix first.\n";
                } else {
                    matrixOperationsMenu();
                }
                break;
                
            case 5: {
                if (matrixCount == 0) {
                    cout << "No matrices to delete!\n";
                    break;
                }
                cout << "\nSelect matrix to delete:\n";
                int idx = selectMatrix();
                if (idx != -1) {
                    delete storedMatrices[idx];
                    // Shift remaining matrices
                    for (int i = idx; i < matrixCount - 1; i++) {
                        storedMatrices[i] = storedMatrices[i + 1];
                    }
                    matrixCount--;
                    cout << "Matrix deleted successfully!\n";
                }
                break;
            }
            
            case 6:
                cout << "\nCleaning up memory...\n";
                // Clean up allocated memory
                for (int i = 0; i < matrixCount; i++) {
                    delete storedMatrices[i];
                }
                cout << "Thank you for using the Linear Algebra Library!\n";
                return 0;
        }
    }
    
    return 0;
}

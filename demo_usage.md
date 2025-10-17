# Linear Algebra Library - Menu Interface Demo

## Overview
The Linear Algebra Library now features a comprehensive menu-driven interface that allows users to interactively create matrices and perform various operations without any hardcoded limitations.

## How to Run
```bash
# Compile the library
g++ -o library_menu library.cpp

# Run the interactive menu
./library_menu
```

## Menu Structure

### Main Menu Options:
1. **Create New Matrix** - Create up to 10 matrices with custom dimensions
2. **Display All Matrices** - View all stored matrices
3. **Matrix Properties** - Check various properties of matrices
4. **Matrix Operations** - Perform operations between matrices
5. **Delete Matrix** - Remove a matrix from storage
6. **Exit** - Clean up and exit the program

### Matrix Properties Submenu:
- Check if Square
- Check if Identity
- Calculate Trace
- Get Dimensions
- Check if Idempotent
- Check if Involutory
- Calculate Determinant

### Matrix Operations Submenu:
- Matrix Addition (A + B)
- Matrix Subtraction (A - B)
- Matrix Multiplication (A * B)
- Scalar Multiplication
- Scalar Division
- Transpose
- Additive Inverse
- Gaussian Elimination
- Column Space
- LU Decomposition
- Symmetric/Skew-Symmetric Decomposition

## Features Added:

### ✅ User Interface Improvements:
- **Interactive Menu System**: No more hardcoded operations
- **Matrix Storage**: Store up to 10 matrices simultaneously
- **Input Validation**: Robust error handling for invalid inputs
- **Clear Navigation**: Hierarchical menu structure with easy navigation

### ✅ Enhanced Functionality:
- **Dynamic Matrix Creation**: User-specified dimensions (up to 100x100)
- **Matrix Selection**: Choose which matrices to operate on
- **Memory Management**: Proper cleanup of allocated matrices
- **Error Prevention**: Checks for compatible matrix dimensions

### ✅ Safety Features:
- **Input Buffer Clearing**: Prevents input stream corruption
- **Dimension Validation**: Ensures valid matrix sizes
- **Operation Validation**: Checks matrix compatibility before operations
- **Memory Cleanup**: Automatic memory deallocation on exit

## Example Usage Session:

```
1. Create a 3x3 matrix
   Menu: 1 → Enter 3 → Enter 3 → Input 9 values

2. Create a 2x2 matrix  
   Menu: 1 → Enter 2 → Enter 2 → Input 4 values

3. View all matrices
   Menu: 2

4. Check if first matrix is square
   Menu: 3 → Select Matrix 1 → Option 1

5. Calculate transpose of first matrix
   Menu: 4 → Option 6 → Select Matrix 1

6. Add two compatible matrices
   Menu: 4 → Option 1 → Select Matrix A → Select Matrix B

7. Exit
   Menu: 6
```

## Benefits:
- **No More Hardcoding**: Users can create any matrices they want
- **Flexible Operations**: Choose which matrices to operate on
- **Educational Value**: Step-by-step exploration of linear algebra concepts
- **Extensible**: Easy to add new operations to the menu system
- **User-Friendly**: Clear prompts and error messages

## Technical Notes:
- The interface maintains backward compatibility with all existing Matrix class methods
- Global matrix storage allows for complex multi-matrix operations
- Input validation prevents crashes from invalid user input
- Memory management ensures no memory leaks
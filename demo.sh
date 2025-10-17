#!/bin/bash

# Demo script for Linear Algebra Library Menu Interface
echo "=== Linear Algebra Library Demo ==="
echo "This script demonstrates the menu interface functionality."
echo ""
echo "To run the program manually:"
echo "1. Compile: g++ -o library library.cpp"
echo "2. Run: ./library"
echo ""
echo "Sample interaction:"
echo "Choice 1: Create matrix"
echo "  Rows: 2"
echo "  Columns: 2" 
echo "  Elements: 1 2 3 4"
echo ""
echo "Choice 2: Matrix addition"
echo "  First matrix: 1 2 3 4"
echo "  Second matrix: 5 6 7 8"
echo "  Result: 6 8 10 12"
echo ""
echo "Choice 6: Transpose (after creating matrix)"
echo "Choice 0: Exit"
echo ""
echo "The program now has a complete menu interface with:"
echo "- 17 different linear algebra operations"
echo "- Input validation"
echo "- Error handling"
echo "- User-friendly prompts"
echo "- Dynamic matrix creation"

# Make this file executable
chmod +x demo.sh
classdef MA1508E
    properties (Constant)
        ZERO_TOLERANCE = 10^-6;
    end

    methods
        function obj = MA1508E()
            % Constructor
            fprintf('MA1508E Toolkit Initialized.\n');
        end
        
        % Chapter 1: Linear Systems

        function res = isValidERO(~, str)
            % Standardize spaces and case
            str = upper(str); 
            str = regexprep(str, '([+-])\s*(R\d)', '$1 $2'); 
            str = regexprep(str, '\s*([+-])\s*', ' $1 '); 
            str = regexprep(str, '\s*S\s*', ' S ');       
            str = strtrim(str);                           
            
            % Patterns that allow optional spaces (\s*) after the operators
            addPattern = "^R\d\s*[+-]\s*[\d./]*\*?R\d$";
            swapPattern = "^R\d\s*S\s*R\d$";
            multiplePattern = "^[+-]?\s*[\d./]*\*?R\d$"; % Now allows "- 1/2R2"
            
            isAdd = ~isempty(regexp(str, addPattern, 'once'));
            isSwap = ~isempty(regexp(str, swapPattern, 'once'));
            isMultiple = ~isempty(regexp(str, multiplePattern, 'once'));
            
            res = isAdd || isSwap || isMultiple;
        end

        function M = generateElemAdd(~, data)
            M = sym(eye(data.size));
            leftVal = data.left(2);
            rightCOE = data.right(1);
            rightVal = data.right(2);
            M(leftVal, rightVal) = rightCOE;
        end

        function M = generateElemSwap(~, data)
            M = sym(eye(data.size));
            leftVal = data.left(2);
            rightVal = data.right(2);
    
            M(leftVal, leftVal) = 0;
            M(rightVal, rightVal) = 0;
            M(leftVal, rightVal) = 1;
            M(rightVal, leftVal) = 1;
        end

        function res = generateElemMatrix(obj, str, n)
            sp = split(str);
            [numParts, ~] = size(sp);
            
            % --- UPDATED SCALAR CASE ---
            % Handles "-1/2R2" (1 part) or "- 1/2R2" (2 parts)
            if numParts == 1 || (numParts == 2 && (sp(1) == "-" || sp(1) == "+"))
                if numParts == 2
                    % If it's "- 1/2R2", combine them back or handle sign
                    opSign = sp(1);
                    target = sp(2);
                else
                    target = sp(1);
                    opSign = "+";
                end
        
                parts = split(target, "R");
                coeStr = strrep(parts(1), '*', '');
                
                if coeStr == ""
                    coe = sym(1);
                else
                    coe = str2sym(coeStr);
                end
                
                if opSign == "-", coe = -coe; end
                
                val = str2double(parts(2));
                res = sym(eye(n));
                res(val, val) = coe;
                return;
            end
            
            % CASE 2: Addition, Subtraction, or Swap (3 parts like "R2 - R1")
            ero = struct('left', zeros(2, 1), 'op', '', 'right', zeros(2, 1), 'size', n);
            
            ero.left = str2double(split(string(sp(1)), "R"));
            ero.op = char(sp(2));
            
            % Process the right side (e.g., "1/2R1" or "R1")
            rightParts = split(string(sp(3)), "R");
            rightCoeStr = strrep(rightParts(1), '*', '');
            
            % Use == "" to correctly identify the "implied 1"
            if rightCoeStr == ""
                ero.right(1) = sym(1);
            else
                ero.right(1) = str2sym(rightCoeStr);
            end
            ero.right(2) = str2double(rightParts(2));
            
            switch ero.op
                case '+'
                    res = obj.generateElemAdd(ero);
                case '-'
                    % Multiply by -1 and treat as addition
                    ero.right(1) = -ero.right(1);
                    res = obj.generateElemAdd(ero);
                case 'S'
                    res = obj.generateElemSwap(ero);
                otherwise
                    fprintf("Something went wrong.\n");
            end
        end

        function res = performERO(obj, A)
            [rows, ~] = size(A);
            M = A;
            command = "";
            while true
                command = input("Enter the elementary row operation: ", 's'); 
                
                if command == "quit" || command == "exit"
                    res = M;
                    return;
                end
                
                % --- SANITIZATION START ---
                command = upper(command); 
                % Fix R2-R1 to R2 - R1
                command = regexprep(command, '([+-])\s*(R\d)', '$1 $2'); 
                % Force single spaces
                command = regexprep(command, '\s*([+-])\s*', ' $1 '); 
                command = regexprep(command, '\s*S\s*', ' S ');      
                command = strtrim(command);                 
                % --- SANITIZATION END ---
                
                if ~obj.isValidERO(command)
                    fprintf("Invalid ERO.\n");
                    continue;
                end
                
                M = obj.generateElemMatrix(command, rows) * M;
                fprintf("--> %s:\n\n", command);
                disp(M);
            end
        end

        function E = getRowTransformation(obj, A, B)
            % Finds the transformation matrix E that represents the sequence of 
            % elementary row operations to turn A into B, such that E*A = B.
            arguments
                obj;
                A;
                B;
            end
            
            [m_A, n_A] = size(A);
            [m_B, n_B] = size(B);
            
            fprintf('\n==================================================\n');
            fprintf('          ROW EQUIVALENCE & TRANSFORMATION        \n');
            fprintf('==================================================\n');
            
            if m_A ~= m_B || n_A ~= n_B
                fprintf('  [!] Error: Matrices must be the exact same dimension.\n');
                fprintf('      A is %dx%d, B is %dx%d.\n', m_A, n_A, m_B, n_B);
                fprintf('==================================================\n\n');
                E = [];
                return;
            end
            
            A_sym = sym(A);
            B_sym = sym(B);
            I = sym(eye(m_A));
            
            % Augment with Identity and find RREF
            rref_A = rref([A_sym, I]);
            rref_B = rref([B_sym, I]);
            
            % Split into the reduced matrix R and the operation matrix E
            R_A = rref_A(:, 1:n_A);
            E_A = rref_A(:, n_A+1:end);
            
            R_B = rref_B(:, 1:n_B);
            E_B = rref_B(:, n_B+1:end);
            
            % Check if they are actually row equivalent
            if ~isequal(simplify(R_A), simplify(R_B))
                fprintf('  [!] Matrices are NOT row equivalent.\n');
                fprintf('      They do not share the same Reduced Row Echelon Form (RREF).\n');
                fprintf('      No sequence of elementary row operations can turn A into B.\n');
                E = [];
            else
                fprintf('  [v] Matrices ARE row equivalent (They share the same RREF).\n');
                
                % Calculate the exact transformation matrix
                % B = (E_B^-1 * E_A) * A
                E = simplify(inv(E_B) * E_A);
                
                fprintf('  ➜ The Net Transformation Matrix (E) where E * A = B is:\n');
                disp(E);
                
                % --- NEW: ELEMENTARY DECOMPOSITION (LU FACTORIZATION) ---
                fprintf('  ➜ Elementary Matrix Decomposition (E = P^T * L * U):\n');
                [L, U, P] = lu(sym(E));
                PT = transpose(P);
                
                fprintf('    1. Permutations / Row Swaps (P^T):\n'); 
                disp(PT);
                fprintf('    2. Forward Row Additions (L):\n'); 
                disp(L);
                fprintf('    3. Row Scaling & Back-Substitutions (U):\n'); 
                disp(U);
                
                fprintf('\n  ➜ DISCRETE ELEMENTARY ROW OPERATIONS:\n');
                fprintf('    (Note: Because E = P^T * L * U, the chronological application on A \n');
                fprintf('     is right-to-left: first U operations, then L, then P^T)\n\n');
                
                stepCounter = 1;
                
                % 1. Extract operations from U (Scaling and Backward Adds)
                for j = n_A:-1:1
                    for i = j:-1:1
                        val = simplify(U(i, j));
                        if i == j && val ~= 1 && val ~= 0
                            fprintf('    Step %d: Scale Row %d ( R%d = %s * R%d )\n', stepCounter, i, i, char(val), i);
                            stepCounter = stepCounter + 1;
                        elseif i ~= j && val ~= 0
                            % Normalize the addition by the pivot to get the pure operation
                            pivot = simplify(U(i,i));
                            if pivot ~= 0
                                addVal = simplify(val / pivot);
                                fprintf('    Step %d: Add to Row %d ( R%d = R%d + (%s)*R%d )\n', stepCounter, i, i, i, char(addVal), j);
                                stepCounter = stepCounter + 1;
                            end
                        end
                    end
                end
                
                % 2. Extract operations from L (Forward Adds)
                for j = 1:n_A
                    for i = j+1:m_A
                        val = simplify(L(i, j));
                        if val ~= 0
                            fprintf('    Step %d: Add to Row %d ( R%d = R%d + (%s)*R%d )\n', stepCounter, i, i, i, char(val), j);
                            stepCounter = stepCounter + 1;
                        end
                    end
                end
                
                % 3. Extract operations from P^T (Row Swaps)
                for i = 1:m_A
                    % Find where the 1 is in this row of P^T
                    swapTarget = find(PT(i, :) == 1);
                    % Only print if it actually moved (and avoid double printing swaps)
                    if swapTarget ~= i && swapTarget > i 
                        fprintf('    Step %d: Swap Rows ( R%d <-> R%d )\n', stepCounter, i, swapTarget);
                        stepCounter = stepCounter + 1;
                    end
                end
                
                if stepCounter == 1
                    fprintf('    (No operations required. Matrices are identical.)\n');
                end
                % --------------------------------------------------------
                
                % Final verification check
                fprintf('  [v] Verification Check: isequal(E * A, B) -> ');
                if isequal(simplify(E * A_sym), B_sym)
                    fprintf('TRUE\n');
                else
                    fprintf('FALSE\n');
                end
            end
            fprintf('==================================================\n\n');
        end

        function [x_gen, basis, x_p] = solveSystem(obj, A, b)
            % Solves Ax = b and returns the parametric general solution
            A_sym = sym(A);
            b_sym = sym(b);
            [rows, cols] = size(A_sym);
            
            % 1. Form Augmented Matrix and find RREF
            Aug = [A_sym, b_sym];
            R = rref(Aug);
            
            % 2. Check for Consistency
            consistent = true;
            for i = 1:size(R, 1)
                if all(R(i, 1:end-1) == 0) && R(i, end) ~= 0
                    consistent = false;
                    break;
                end
            end
            
            if ~consistent
                fprintf('\n[!] RESULT: System is INCONSISTENT (No Solution).\n');
                x_gen = []; basis = []; x_p = [];
                return;
            end
            
            % 3. Extract Particular Solution (setting free variables to 0)
            x_p = sym(zeros(cols, 1));
            [~, pivotCols] = rref(double(A_sym)); 
            for i = 1:length(pivotCols)
                x_p(pivotCols(i)) = R(i, end);
            end
            
            % 4. Extract Null Space Basis (Homogeneous Solution)
            % sym/null automatically returns the rational basis
            basis = null(A_sym); 
            [~, k] = size(basis);
            
            % 5. Format outputs
            if k > 0
                s_vars = sym('s', [k, 1], 'real');
                x_gen = simplify(x_p + (basis * s_vars));
            else
                x_gen = x_p;
            end
            
            % --- BEAUTIFIED OUTPUT ---
            fprintf('\n==================================================\n');
            fprintf('          PARAMETRIC SYSTEM SOLUTION              \n');
            fprintf('==================================================\n');
            if k == 0
                fprintf('RESULT: Unique Solution Found\n');
                disp(x_p);
            else
                fprintf('RESULT: Infinite Solutions (Parametric Form)\n');
                fprintf('Particular Solution (x_p):\n');
                disp(x_p);
                fprintf('Null Space Basis (v_i):\n');
                disp(basis);
                fprintf('--------------------------------------------------\n');
                fprintf('General Equation: x = x_p');
                for i = 1:k
                    fprintf(' + s%d*v%d', i, i);
                end
                fprintf('\n');
            end
            fprintf('==================================================\n\n');
        end

        % Chapter 2: Matrices

        function check = isSquare(~, A)
            % Evaluates if a matrix is square based purely on its dimensions
            arguments
                ~; % Ignore obj since it doesn't use class properties
                A;
            end
            
            [rows, cols] = size(A);
            if rows == cols
                check = true;
                fprintf('  [v] Matrix is square (%d x %d).\n', rows, cols);
            else
                check = false;
                fprintf('  [!] Matrix is NOT square (%d x %d).\n', rows, cols);
            end
        end

        function evaluateIMT(obj, A)
            % Evaluates the Invertible Matrix Theorem based on size and determinant
            arguments
                obj;
                A;
            end
            
            fprintf('\n==================================================\n');
            fprintf('        INVERTIBLE MATRIX THEOREM (IMT) ANALYSIS  \n');
            fprintf('==================================================\n');
            
            % 1. Size Check
            [n, cols] = size(A);
            if n ~= cols
                fprintf('[!] ERROR: Matrix is %dx%d.\n', n, cols);
                fprintf('The Invertible Matrix Theorem ONLY applies to square (n x n) matrices.\n');
                fprintf('==================================================\n\n');
                return; % Exit the function early
            end
            
            % 2. Determinant Check
            matrix_det = det(sym(A));
            fprintf('➜ Dimension: %dx%d\n', n, n);
            fprintf('➜ Determinant: %s\n\n', char(matrix_det));
            
            if matrix_det ~= 0
                fprintf('Because det(A) ≠ 0, ALL of the following are PROVEN TRUE:\n');
                fprintf('  1. Matrix A is invertible (A^-1 exists).\n');
                fprintf('  2. The columns of A form a basis for R^%d.\n', n);
                fprintf('  3. The equation Ax = 0 has ONLY the trivial solution (x = 0).\n');
                fprintf('  4. The equation Ax = b has exactly one solution for every b in R^%d.\n', n);
                fprintf('  5. The linear transformation x ↦ Ax is one-to-one and onto.\n');
                fprintf('  6. Zero is NOT an eigenvalue of A.\n');
            else
                fprintf('Because det(A) = 0, ALL of the following are PROVEN TRUE:\n');
                fprintf('  1. Matrix A is singular (NOT invertible).\n');
                fprintf('  2. The columns of A are linearly dependent.\n');
                fprintf('  3. The equation Ax = 0 has non-trivial (infinite) solutions.\n');
                fprintf('  4. The equation Ax = b may have zero or infinite solutions (never unique).\n');
                fprintf('  5. Zero IS an eigenvalue of A.\n');
            end
            fprintf('==================================================\n\n');
        end

        function checkSymmetry(obj, A)
            % Evaluates if a matrix is symmetric (A = A^T) or skew-symmetric (A = -A^T)
            arguments
                obj;
                A;
            end
            
            [rows, cols] = size(A);
            if rows ~= cols
                fprintf('  [!] Dimension Error: Matrix is %dx%d. Symmetry requires a square matrix.\n', rows, cols);
                return;
            end
            
            A_sym = sym(A);
            A_trans = transpose(A_sym);
            
            fprintf('\n==================================================\n');
            fprintf('               SYMMETRY ANALYSIS                  \n');
            fprintf('==================================================\n');
            
            if isequal(simplify(A_sym), simplify(A_trans))
                fprintf('  [v] Matrix is SYMMETRIC (A = A^T).\n\n');
                fprintf('  PROVEN THEORETICAL RELATIONSHIPS:\n');
                fprintf('  1. Spectral Theorem: A is orthogonally diagonalizable.\n');
                fprintf('  2. All eigenvalues of A are real numbers.\n');
                fprintf('  3. Eigenvectors from different eigenspaces are orthogonal.\n');
            elseif isequal(simplify(A_sym), simplify(-A_trans))
                fprintf('  [v] Matrix is SKEW-SYMMETRIC (A = -A^T).\n\n');
                fprintf('  PROVEN THEORETICAL RELATIONSHIPS:\n');
                fprintf('  1. All main diagonal entries are exactly zero.\n');
                fprintf('  2. If the dimension is odd, the determinant is zero.\n');
            else
                fprintf('  [!] Matrix is NEITHER symmetric nor skew-symmetric.\n');
            end
            fprintf('==================================================\n\n');
        end

        % Chapter 3: Vector Spaces & Change of Basis

        function check = isInSpan(~, S, v)
            % Checks if vector(s) v exist within the span of set S
            arguments
                ~;
                S; % Matrix where columns are the spanning set
                v; % Vector(s) to check
            end
            
            fprintf('\n==================================================\n');
            fprintf('                 SPAN ANALYSIS                    \n');
            fprintf('==================================================\n');
            
            rankS = rank(sym(S));
            rankSv = rank(sym([S, v]));
            
            if rankS == rankSv
                check = true;
                fprintf('  [v] The vector IS in the span.\n');
                fprintf('      Rank(S) = %d matches Rank([S, v]) = %d.\n', rankS, rankSv);
            else
                check = false;
                fprintf('  [!] The vector is NOT in the span.\n');
                fprintf('      Rank(S) = %d does not match Rank([S, v]) = %d.\n', rankS, rankSv);
            end
            fprintf('==================================================\n\n');
        end

        function compareSubspaces(~, U, W)
            % Compares two subspaces to determine containment and equality
            arguments
                ~;
                U; % Matrix where columns span Subspace U
                W; % Matrix where columns span Subspace W
            end
            
            fprintf('\n==================================================\n');
            fprintf('              SUBSPACE CONTAINMENT                \n');
            fprintf('==================================================\n');
            
            rU = rank(sym(U));
            rW = rank(sym(W));
            rUW = rank(sym([U, W]));
            
            fprintf('  ➜ Dimension of U: %d\n', rU);
            fprintf('  ➜ Dimension of W: %d\n', rW);
            fprintf('  ➜ Dimension of Union Span [U, W]: %d\n\n', rUW);
            
            if rU == rW && rUW == rU
                fprintf('  [v] U and W are the EXACT SAME subspace (U = W).\n');
            elseif rUW == rW
                fprintf('  [v] U is a PROPER SUBSPACE of W (U ⊂ W).\n');
                fprintf('      Every vector in U can be built using vectors from W.\n');
            elseif rUW == rU
                fprintf('  [v] W is a PROPER SUBSPACE of U (W ⊂ U).\n');
                fprintf('      Every vector in W can be built using vectors from U.\n');
            else
                fprintf('  [!] Neither is a subspace of the other.\n');
                fprintf('      They may intersect, but one does not fully contain the other.\n');
            end
            fprintf('==================================================\n\n');
        end

        function intersectionBasis = intersectSubspaces(~, U, W)
            % Finds the basis for the intersection of two subspaces
            arguments
                ~;
                U; % Matrix where columns span Subspace U
                W; % Matrix where columns span Subspace W
            end
            
            fprintf('\n==================================================\n');
            fprintf('             SUBSPACE INTERSECTION                \n');
            fprintf('==================================================\n');
            
            % Setup block matrix [U, -W] to solve Ux - Wy = 0
            blockMatrix = [sym(U), -sym(W)];
            nullBasis = null(blockMatrix);
            
            if isempty(nullBasis)
                fprintf('  [i] The intersection is only the ZERO VECTOR {0}.\n');
                fprintf('      The subspaces are disjoint.\n');
                intersectionBasis = [];
            else
                % Extract the 'x' portion of the null space vectors
                colsU = size(U, 2);
                x_components = nullBasis(1:colsU, :);
                
                % Multiply U by the x_components to get the physical basis vectors
                intersectionBasis = simplify(sym(U) * x_components);
                
                dimIntersect = size(intersectionBasis, 2);
                fprintf('  [v] The intersection is a subspace of dimension %d.\n', dimIntersect);
                fprintf('  ➜ Basis for the intersection:\n');
                disp(intersectionBasis);
            end
            fprintf('==================================================\n\n');
        end

        function check = isBasis(obj, V)
            % Evaluates if a set of column vectors forms a valid basis in Rn
            arguments
                obj;
                V;
            end
            
            [rows, cols] = size(V);
            if rows ~= cols
                fprintf('  [!] Check Failed: Not a square matrix. You need %d vectors for R^%d.\n', rows, rows);
                check = false;
            elseif det(sym(V)) == 0
                fprintf('  [!] Check Failed: Vectors are linearly dependent (Determinant is 0).\n');
                check = false;
            else
                fprintf('  [v] Check Passed: Set is a valid basis.\n');
                check = true;
            end
        end

        function P = getTransitionMatrix(obj, B, C)
            % Calculates the transition matrix from basis B to basis C
            arguments
                obj;
                B;
                C;
            end
            
            [rB, cB] = size(B);
            [rC, cC] = size(C);
            
            if rB ~= rC || cB ~= cC
                error('Dimension Mismatch: Bases B and C must have the same dimensions.');
            end
            if det(sym(C)) == 0
                error('Math Error: Target set C is linearly dependent and cannot form a basis.');
            end
            
            P = inv(sym(C)) * sym(B);
            
            fprintf('\n==================================================\n');
            fprintf('        TRANSITION MATRIX: [I]_{C <- B}           \n');
            fprintf('==================================================\n');
            disp(P);
            fprintf('==================================================\n\n');
        end

        function v_C = changeVectorBasis(obj, v_B, B, C)
            % Converts coordinates of a vector from basis B to basis C
            arguments
                obj;
                v_B;
                B;
                C;
            end
            
            % Force v_B to be a column vector just in case the student passed a row vector
            if isrow(v_B)
                v_B = v_B.';
            end
            
            P = obj.getTransitionMatrix(B, C);
            v_C = P * sym(v_B);
            
            fprintf('\n➜ Coordinates of vector in Target Basis C:\n');
            disp(v_C);
        end

        function A_C = similarityTransform(obj, A_B, B, C)
            % Converts a linear transformation matrix A from basis B to basis C
            arguments
                obj;
                A_B;
                B;
                C;
            end
            
            [rA, cA] = size(A_B);
            if rA ~= cA
                error('Dimension Error: Transformation matrix A must be square.');
            end
            
            % For similarity, P is the transition from C to B
            P = obj.getTransitionMatrix(C, B); 
            A_C = inv(P) * sym(A_B) * P;
            
            fprintf('\n==================================================\n');
            fprintf('      SIMILARITY TRANSFORM: A_C = P^-1 * A_B * P  \n');
            fprintf('==================================================\n');
            disp(A_C);
            fprintf('==================================================\n\n');
        end

        % Chapter 4: Subspaces

        function LI = leftInverse(~, A)
            [rows, cols] = size(A);
            if rows < cols 
                fprintf("The matrix does not have a left inverse!\n");
                return;
            end
            if rows == cols
                fprintf("The matrix is square, its inverse is:\n");
                LI = inv(A);
                return
            end
            
            fprintf("The left inverse of the matrix exists.\n")
            LI = inv(A' * A) * A';
        end

        function RI = rightInverse(~, A)
            [rows, cols] = size(A);
            if rows > cols 
                fprintf("The matrix does not have a right inverse!\n");
                return;
            end
            if rows == cols
                fprintf("The matrix is square, its inverse is:\n");
                disp(inv(A));
                return
            end
            
            fprintf("The right inverse of the matrix exists.\n")
            RI = A' * inv(A * A');
        end
        
        function getMatrixSpaces(~, A)
            % Finds the basis for Col(A), Row(A), and Null(A)
            A_sym = sym(A);
            [R, pivots] = rref(A);
            
            fprintf('\n--- COLUMN SPACE BASIS (Pivot Columns of A) ---\n');
            disp(A_sym(:, pivots));
            
            fprintf('--- ROW SPACE BASIS (Non-zero rows of RREF) ---\n');
            disp(R(1:length(pivots), :));
            
            fprintf('--- NULL SPACE BASIS (Solution to Ax = 0) ---\n');
            disp(null(A_sym));
        end

        % Chapter 5: Orthogonal Projection

        function checkOrthogonal(obj, A)
            % Evaluates if a matrix is orthogonal (A^T = A^-1)
            arguments
                obj;
                A;
            end
            
            [rows, cols] = size(A);
            if rows ~= cols
                fprintf('  [!] Dimension Error: Matrix is %dx%d. Orthogonality requires a square matrix.\n', rows, cols);
                return;
            end
            
            A_sym = sym(A);
            I = sym(eye(rows));
            
            fprintf('\n==================================================\n');
            fprintf('             ORTHOGONALITY ANALYSIS               \n');
            fprintf('==================================================\n');
            
            % A matrix is orthogonal if A^T * A = I. 
            % We use this instead of A^T == A^-1 because calculating inverses is slow and prone to symbolic errors.
            if isequal(simplify(transpose(A_sym) * A_sym), I)
                fprintf('  [v] Matrix is ORTHOGONAL (A^T = A^-1).\n\n');
                fprintf('  PROVEN THEORETICAL RELATIONSHIPS:\n');
                fprintf('  1. The columns form an orthonormal basis for R^%d.\n', rows);
                fprintf('  2. The rows form an orthonormal basis for R^%d.\n', rows);
                fprintf('  3. The determinant is exactly %s.\n', char(det(A_sym)));
                fprintf('  4. The transformation preserves lengths (||Ax|| = ||x||).\n');
                fprintf('  5. The transformation preserves angles (Ax · Ay = x · y).\n');
            else
                fprintf('  [!] Matrix is NOT orthogonal (A^T ≠ A^-1).\n');
            end
            fprintf('==================================================\n\n');
        end

        function res = isOrthogonalSet(obj, S)    
            [~, w] = size(S);
            for i = 1:w
                for j = i:w
                    if (i == j) continue; end
                    if dot(S(:, i), S(:, j)) > obj.ZERO_TOLERANCE
                        fprintf("The set is not orthogonal\n");
                        res = false;
                        return
                    end
                end
            end
                fprintf("The set is orthogonal\n");
                res = true;
        end
        
        function isOrthonormalSet(obj, S)  
            if (~obj.isOrthogonalSet(S)) return; end
            
            [~, w] = size(S);
            for i = 1:w
                if (abs(norm(S(:, i)) - 1) > obj.ZERO_TOLERANCE)
                    fprintf("The set is not orthonormal\n");
                    return
                end
            end
            fprintf("The set is orthonormal\n");
        end
        
        function isOrthogonalTo(~, S, target)
            %M: The matrix containing the set of vectors to check against
            %target: The vector to check if orthogonal
            [~, w] = size(S);
            for i = 1:w
                if dot(S(:, i), target) ~= 0
                    fprintf("The vector is not orthogonal to the matrix\n");
                    return
                end
            end
            fprintf("The vector is orthogonal to the matrix\n");
        end  
        
        function ON = toOrthonormalSet(~, OG)
            [rows, cols] = size(OG);
            ON = zeros(rows, cols);
            for i = 1:cols
                ON(:, i) = OG(:, i) / norm(OG(:, i));
            end
        end

        function T = dotWithSet(~, S, v)
            [rowsS, colsS] = size(S);
            [rowsV, ~] = size(v);
            
            if rowsS ~= rowsV
               fprintf("Error: v does not have the same number of rows as S!\n");
               return;
            end
            
            T = zeros(1, colsS);
            for i = 1:colsS
                T(:, i) = dot(S(:, i), v);
            end
        end

        function b = orthogonalProj(~, A, u)
            %A: The matrix to be projected onto
            %u: The vector to project
            %Returns: The projection of u onto A
            [~, w] = size(u);
            if (w ~= 1) 
                fprintf("Warning: The input vector u has more than one column\n"); 
            end
            b = A * inv(A'*A) * A' * u;
        end
        
        function [e, r] = gramSchmidt(obj, v, showSteps)
            % Set default for showSteps if not provided
            if nargin < 3
                showSteps = true;
            end
            
            [M, N] = size(v);
            v_sym = sym(v);        % Force symbolic for exact surds
            e = sym(zeros(M, N));  % Orthonormal vectors (Q matrix)
            r = sym(zeros(N, N));  % Upper triangular coefficients (R matrix)
            
            if showSteps
                fprintf("\n==================================================\n");
                fprintf("       GRAM-SCHMIDT ORTHONORMALIZATION        \n");
                fprintf("==================================================\n");
            end
        
            for j = 1:N
                % Start with the original vector v_j
                temp_u = v_sym(:, j);
                
                if showSteps
                    fprintf("\nStep %d: Processing Input Vector v_%d\n", j, j);
                    disp(v_sym(:, j));
                end
        
                % 1. Orthogonalization (Subtract projections onto previous e vectors)
                for i = 1:j-1
                    % Calculate R_ij (projection of original v_j onto current e_i)
                    % This must be done using the orthonormal vectors for V = QR
                    r(i, j) = simplify(dot(e(:, i), v_sym(:, j)));
                    
                    % Subtract the projection component
                    temp_u = temp_u - (r(i, j) * e(:, i));
                    
                    if showSteps
                        fprintf("  > Projection onto orthonormal e_%d:\n", i);
                        fprintf("    Coefficient R(%d,%d) = <e_%d, v_%d> = %s\n", i, j, i, j, string(r(i,j)));
                    end
                end
                
                % 2. Calculate the norm of the remaining orthogonal vector
                u_orth = simplify(temp_u);
                mag = simplify(sqrt(dot(u_orth, u_orth)));
                r(j, j) = mag; % This is the diagonal of the R matrix
                
                % 3. Normalization (Generating the Orthonormal basis)
                if mag ~= 0
                    e(:, j) = simplify(u_orth / mag);
                    
                    if showSteps
                        fprintf("  > Resulting Orthogonal Vector u_%d:\n", j);
                        disp(u_orth);
                        fprintf("  > Normalizing... Magnitude ||u_%d|| = %s\n", j, string(mag));
                        fprintf("  > Final Orthonormal Vector e_%d (Surd Form):\n", j);
                        disp(e(:, j));
                    end
                else
                    % Handle cases where vectors are linearly dependent
                    fprintf("  > Warning: Vector v_%d is linearly dependent. Result is 0.\n", j);
                    e(:, j) = sym(zeros(M, 1));
                end
            end
            
            if showSteps
                fprintf("==================================================\n\n");
            end
        end

        function orthoBasis = getOrthogonalBasis(obj, A)
            % Finds an orthogonal basis for the subspace spanned by the columns of A
            
            % 1. Find the linearly independent basis vectors (Pivot columns)
            [~, pivotCols] = rref(double(A)); % Use double() to safely find pivots
            basis = sym(A(:, pivotCols));     % Extract the exact symbolic columns
            
            fprintf("\nThe original subspace is defined by pivot columns: ");
            disp(pivotCols);
            
            fprintf("The linearly independent basis for this subspace is:\n");
            disp(basis);
            
            % 2. Apply Gram-Schmidt to orthogonalize the basis
            [orthoBasis, ~] = obj.gramSchmidt(basis);
            
            fprintf("The Orthogonal Basis for the subspace is:\n");
            disp(orthoBasis);
        end

        function P = getProjectionMatrix(~, A)
            % Returns the projection matrix P = A(A^T A)^-1 A^T
            A_sym = sym(A);
            % Use transpose() for symbolic compatibility
            P = A_sym * inv(transpose(A_sym) * A_sym) * transpose(A_sym);
            P = simplify(P);
            
            fprintf('\n--- Projection Matrix (P) ---\n');
            disp(P);
            fprintf('Properties: P^2 = P and P^T = P\n');
        end
        
        function [x_gen, basis, x_p] = calcLSS(obj, A, b)
            % 1. Symbolic Conversion
            A_sym = sym(A);
            b_sym = sym(b);
            [rowsA, colsA] = size(A_sym);
        
            % 2. Form Normal Equations 
            % Using transpose() instead of ' to avoid "Invalid use of operator"
            left = transpose(A_sym) * A_sym;
            right = transpose(A_sym) * b_sym;
            
            % 3. Solve via RREF
            Aug = [left, right];
            R = rref(Aug);
            
            % 4. Extract Particular Solution (x_p)
            x_p = sym(zeros(colsA, 1));
            [numRows, ~] = size(R);
            pivotCount = 1;
            for j = 1:colsA
                if pivotCount <= numRows && R(pivotCount, j) == 1
                    x_p(j) = R(pivotCount, end);
                    pivotCount = pivotCount + 1;
                end
            end
        
            % 5. Homogeneous Solution (Null Space)
            basis = null(left);
            [~, k] = size(basis);
        
            % 6. Simplify Building Blocks
            x_p = simplify(x_p);
            basis = simplify(basis);
            
            % 7. Generate x_gen for return (not for display)
            if k > 0
                % Create s_vars as a column vector to avoid transpose operator
                s_vars = sym('s', [k, 1], 'real');
                x_gen = simplify(x_p + (basis * s_vars));
            else
                x_gen = x_p;
            end
        
            % --- FINAL BEAUTIFIED OUTPUT ---
            fprintf('\n==================================================\n');
            fprintf('             LEAST SQUARES ANALYSIS               \n');
            fprintf('==================================================\n');
            fprintf('System: (A^T * A)x = A^T * b\n');
        
            if k == 0
                fprintf('\nRESULT: Unique Solution Found\n');
                fprintf('--------------------------------------------------\n');
                fprintf('x_lss =\n');
                disp(x_p);
            else
                fprintf('\nRESULT: Infinite Solutions (Parametric Form)\n');
                fprintf('--------------------------------------------------\n');
                fprintf('Particular Solution (x_p):\n');
                disp(x_p);
                
                fprintf('Null Space Basis (v_i):\n');
                for i = 1:k
                    fprintf('  Vector v%d:\n', i);
                    disp(basis(:, i));
                end
                
                fprintf('--------------------------------------------------\n');
                fprintf('Parametric Equation:\n');
                fprintf('  x = x_p');
                % Loop-based printing avoids all string concatenation errors
                for i = 1:k
                    fprintf(' + s%d*v%d', i, i);
                end
                fprintf('\n');
            end
            fprintf('==================================================\n\n');
        end
        
        % Chapter 6: Diagonalisation

        function ev = getEigenvalues(~, A)
            syms x;
            px = charpoly(A, x);
            f = factor(px, x);
            [~, nFactors] = size(f);
            
            j = 0;
            for i = 1:nFactors
                roots = solve(f(i) == 0, x);
                nRoots = polynomialDegree(f(i));
                
                for k = 1:nRoots
                    ev(i + j) = roots(k);
                    % Complex roots
                    if nRoots > 1 && k < nRoots
                        j = j + 1;
                    end
                end
            end
        end
        
        function res = isAssociatedEigenvalue(obj, A, lambda)
            ev = obj.getEigenvalues(A);
            valid = false;
            for i = 1:length(ev) % FIXED: Changed size(ev') to length(ev)
                if lambda == ev(i)
                    valid = true;
                    break;
                end
            end
            if ~valid
                fprintf("Eigenvalue %i is not associated with the matrix!\n", lambda);
                fprintf("The following eigenvalues are associated with the matrix: ");
                disp(ev);
            end
            res = valid;
        end
        
        function [V, genSol] = getEigenvector(obj, A, lambda, output)
            arguments
                obj;
                A; 
                lambda; 
                output logical = true;
            end
            
            [rows, cols] = size(A);
            if rows ~= cols
                fprintf("Input matrix is not square!\n");
                return;
            end
            if ~obj.isAssociatedEigenvalue(A, lambda)
                return;
            end
            
            % 1. Construct the characteristic matrix (lambda*I - A)
            M = lambda * sym(eye(cols)) - sym(A); 
            
            % 2. Get the 'Rational' Basis (matches RREF manual working)
            V = null(M); 
            [~, k] = size(V); % k is the geometric multiplicity
            
            % 3. Generate the general solution: x = s1*v1 + s2*v2...
            s = sym('s', [1 k], 'real'); % Creates s1, s2, ... sk
            genSol = V * s';      % Matrix-vector multiplication for the general form
            
            if output
                fprintf("\n--- Eigenspace Analysis for lambda = %s ---\n", string(lambda));
                fprintf("Characteristic Matrix (lambda*I - A) reduced to RREF:\n");
                disp(rref(M));
                
                fprintf("Basis for the Eigenspace:\n");
                disp(V);
                
                fprintf("Parameterized General Solution (x = s1*v1 + ...):\n");
                disp(genSol);
            end
        end

        function checkDiagonalizable(obj, A)
            fprintf('\n==================================================\n');
            fprintf('          DIAGONALIZABILITY DIAGNOSTIC            \n');
            fprintf('==================================================\n');
            
            A_sym = sym(A);
            % Get ALL eigenvalues including repetitions
            full_evals = eig(A_sym);
            % Get unique eigenvalues to loop through
            unique_evals = unique(full_evals);
            
            is_diag = true;
            for i = 1:length(unique_evals)
                curr_lambda = unique_evals(i);
                
                % --- FIX: Use isAlways for symbolic comparison ---
                % This correctly counts multiplicities
                curr_am = sum(isAlways(full_evals == curr_lambda));
                
                % Geometric Multiplicity (Dimension of Null Space)
                [V, ~] = obj.getEigenvector(A_sym, curr_lambda, false);
                [~, curr_gm] = size(V);
                
                fprintf('λ = %s:\n', char(curr_lambda));
                fprintf('  ➜ Algebraic Multiplicity (AM): %d\n', curr_am);
                fprintf('  ➜ Geometric Multiplicity (GM): %d\n', curr_gm);
                
                if curr_am ~= curr_gm
                    fprintf('  [!] DEFECTIVE: AM > GM. This eigenvalue prevents diagonalization.\n');
                    is_diag = false;
                else
                    fprintf('  [v] Healthy: AM = GM.\n');
                end
            end
            
            if is_diag
                fprintf('\nCONCLUSION: Matrix is DIAGONALIZABLE.\n');
            else
                fprintf('\nCONCLUSION: Matrix is NOT diagonalizable (Defective).\n');
            end
            fprintf('==================================================\n\n');
        end

        function q = solveMarkovSteadyState(~, P)
            % Finds the steady-state vector q for a transition matrix P
            % Solves (P - I)q = 0 and normalizes so sum(q) = 1
            [n, ~] = size(P);
            P_sym = sym(P);
            
            % Solve the homogeneous system (P - I)x = 0
            M = P_sym - eye(n);
            v = null(M); 
            
            if isempty(v)
                fprintf('[!] Error: No steady state found. Is this a valid Stochastic Matrix?\n');
                return;
            end
            
            % Normalize so the sum of components is 1
            q = simplify(v / sum(v));
            
            fprintf('\n--- Markov Chain Steady-State Vector (q) ---\n');
            disp(q);
            fprintf('Sum of components: %s\n', char(sum(q)));
        end

        function computeSVD(obj, A)
            % Analyzes the Singular Value Decomposition components for any m x n matrix
            arguments
                obj;
                A;
            end
            
            A_sym = sym(A);
            [m, n] = size(A_sym);
            
            fprintf('\n==================================================\n');
            fprintf('      SINGULAR VALUE DECOMPOSITION (SVD)          \n');
            fprintf('==================================================\n');
            
            if m == n
                fprintf('  [i] Dimension: %dx%d (Square).\n', m, n);
            else
                fprintf('  [!] Dimension: %dx%d (Rectangular).\n', m, n);
                fprintf('      Standard eigenvalues are strictly UNDEFINED.\n');
                fprintf('      Analyzing A^T*A and A*A^T instead...\n');
            end
            
            % 1. Analyze A^T * A (Right Singular Vectors)
            AtA = transpose(A_sym) * A_sym;
            fprintf('\n--- 1. Right Singular Matrix (V) from A^T * A ---\n');
            [V, D_AtA] = eig(AtA);
            disp(simplify(V));
            
            % 2. Analyze A * A^T (Left Singular Vectors)
            AAt = A_sym * transpose(A_sym);
            fprintf('\n--- 2. Left Singular Matrix (U) from A * A^T ---\n');
            [U, ~] = eig(AAt);
            disp(simplify(U));
            
            % 3. Extract Singular Values (Sigma)
            fprintf('\n--- 3. Non-Zero Singular Values (σ) ---\n');
            fprintf('       (Square roots of the non-zero eigenvalues)\n');
            
            evals = diag(D_AtA);
            % Filter for non-zero eigenvalues (handling symbolic floating point quirks)
            non_zero_evals = evals(evals > 1e-10); 
            sigma = simplify(sqrt(non_zero_evals));
            
            disp(sigma);
            fprintf('==================================================\n\n');
        end

        % Chapter 7: System of Linear Differential Equations

        function v2 = getGeneralisedEigenvector(obj, A, lambda)
            [rows, cols] = size(A);
            if rows ~= cols
                fprintf("Input matrix is not square!\n");
                return;
            end
            
            if ~obj.isAssociatedEigenvalue(A, lambda)
                return;
            end
            M = sym(A) - lambda * sym(eye(cols));
            v1 = obj.getEigenvector(A, lambda, false);
            
            % Calculate RREF of augmented matrix
            R = rref([M v1]);
            fprintf("\nThe augmented matrix [M | v1] is reduced to:\n")
            disp(R);
            
            % Automate the extraction (set free variables to 0)
            [~, pivotCols] = rref(double(M)); % Cast to double just to find the pivots
            v2 = sym(zeros(cols, 1)); % Start with a blank zero vector
            
            % Map the solutions to the pivot positions
            for i = 1:length(pivotCols)
                v2(pivotCols(i)) = R(i, end); 
            end
            
            fprintf("By setting free variables to 0, a valid generalised eigenvector is:\n");
            disp(v2);
        end

        function y0 = generateInitialConditions(obj, n)
            % Prompts user for initial values at t=0 and returns a column vector
            fprintf('\n--- Entering Initial Conditions for y(0) ---\n');
            y0 = sym(zeros(n, 1)); 
            for i = 1:n
                val = input(sprintf('  Enter value for y%i(0): ', i));
                y0(i) = sym(val);
            end
        end
        
        function res = solveDifferentialSystem(obj, A, isInitial)
            % Solves systems of linear ODEs: y' = Ay
            arguments
                obj;
                A;
                isInitial logical = false;
            end

            [n, ~] = size(A);
            syms t;
            
            % Step 1: Characteristic Analysis (Eigenvalues)
            evals = eig(sym(A));
            
            % Step 2: Create array of symbolic functions
            y = sym(zeros(n, 1)); 
            for i = 1:n
                % Create symfun in workspace and store in array
                eval(sprintf('syms y%d(t)', i));
                eval(sprintf('y(i) = y%d;', i));
            end
            
            % Step 3: Define the ODE system (y' = Ay)
            ode = diff(y, t) == A * y;

            if ~isInitial
                % Solve for the General Solution
                res = dsolve(ode);
            else
                % Interactive step for y0 vector
                y0_vec = obj.generateInitialConditions(n);
                
                % FIX: Call y1(0), y2(0) directly from workspace using eval
                % This bypasses MATLAB's habit of breaking functions inside arrays
                conds = [];
                for i = 1:n
                    y_at_zero = eval(sprintf('y%d(0)', i)); 
                    conds = [conds, y_at_zero == y0_vec(i)];
                end
                
                % Solve the Initial Value Problem
                res = dsolve(ode, conds);
            end

            % --- ORGANIZED WORKFLOW OUTPUT ---
            fprintf('\n==================================================\n');
            fprintf('        SYSTEM ANALYSIS: y'' = Ay\n');
            fprintf('==================================================\n');
            
            fprintf('1. EIGENVALUES (λ):\n');
            disp(evals);
            fprintf('--------------------------------------------------\n');

            fprintf('2. ANALYTICAL SOLUTION:\n');
            
            % Handle structure return from dsolve
            if isstruct(res)
                sol_fields = fieldnames(res);
                for i = 1:length(sol_fields)
                    fprintf('  ➜ %s(t) = \n', sol_fields{i});
                    disp(res.(sol_fields{i}));
                    if i < length(sol_fields)
                        fprintf('  \n');
                    end
                end
            else
                % Fallback for 1x1 systems or single variables
                disp(res);
            end
            fprintf('==================================================\n\n');
        end

        % --- Power User Functions ---
        
        function [P, D] = getPD(~, A)
            % Instantly builds the P and D matrices for A = PDP^-1
            [rows, cols] = size(A);
            if rows ~= cols
                fprintf("Error: Matrix must be square!\n");
                return;
            end
            
            % The 'jordan' function safely handles both normal and defective matrices
            [P_sym, D_sym] = jordan(sym(A)); 
            
            fprintf("\nThe Diagonal (or Jordan) matrix 'D' containing eigenvalues:\n");
            disp(D_sym);
            fprintf("The Change of Basis matrix 'P' containing exact eigenvectors:\n");
            disp(P_sym);
            
            P = P_sym;
            D = D_sym;
        end
        
        function plotPhasePortrait(~, A)
            % Plots the 2D vector field and eigenvectors for x' = Ax
            if size(A,1) ~= 2 || size(A,2) ~= 2
                fprintf("Error: Phase portraits can only be plotted for 2x2 matrices.\n");
                return;
            end
            
            A_num = double(A); % Graphics require standard floating-point numbers
            [x, y] = meshgrid(-10:1.5:10, -10:1.5:10);
            
            % Calculate the directional vectors for the grid
            u = A_num(1,1)*x + A_num(1,2)*y;
            v = A_num(2,1)*x + A_num(2,2)*y;
            
            figure('Name', 'Chapter 7: Phase Portrait');
            quiver(x, y, u, v, 'b', 'AutoScaleFactor', 1.5);
            hold on;
            
            % Calculate and plot real eigenvectors as red lines
            [V, D] = eig(A_num);
            if isreal(diag(D))
                for i = 1:2
                    % Extend the eigenvector to draw a full line across the graph
                    p_x = [-15*V(1,i), 15*V(1,i)];
                    p_y = [-15*V(2,i), 15*V(2,i)];
                    plot(p_x, p_y, 'r', 'LineWidth', 2);
                end
                legend('Vector Field (Movement)', 'Eigenvectors (Straight-Line Solutions)', 'Location', 'best');
            else
                legend('Vector Field (Spirals - Complex Eigenvalues)', 'Location', 'best');
            end
            
            title("Phase Portrait of x' = Ax");
            xlabel('x_1'); ylabel('x_2');
            axis([-10 10 -10 10]);
            grid on; hold off;
        end

    end
end

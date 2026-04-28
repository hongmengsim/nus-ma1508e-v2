classdef MA1508E
    properties (Constant)
        ZERO_TOLERANCE = 10^-6;
    end
    methods
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
        % Chapter 2: Matrices
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
        
        % Chapter 5: Orthogonal Projection
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

        % Chapter 7: System of Linear Differential Equations
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

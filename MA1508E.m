classdef MA1508E
    properties (Constant)
        ZERO_TOLERANCE = 10^-6;
    end
    methods
        % Chapter 1: Linear Systems
        function res = isValidERO(~, str)
            % Match regular expression, includes floating point coefficient            
            % ERO Examples:
            % Add - "R1 + 3R2"; "R4 - 1R2"
            % Swap - "R2 S R4"
            % Multiple - "3R1"; "0.25R3"
            [addStart, addEnd] = regexp(str, "R\d [+-] \d+.?\d*R\d");
            [swapStart, swapEnd] = regexp(str, "R\d S R\d");
            [multipleStart, multipleEnd] = regexp(str, "-?\d+.?\d*R\d");
            
            % Check if input string is exactly the length of the match
            len = strlength(str);
            isAdd = any(addStart == 1) && ((addEnd - addStart + 1) == len);
            isSwap = any(swapStart == 1) && ((swapEnd - swapStart + 1) == len);
            isMultiple = any(multipleStart == 1) && ((multipleEnd - multipleStart + 1) == len);
            if ~isAdd && ~isSwap && ~isMultiple
                res = false;
                return;
            end
            res = true;
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
            % Process input string
            sp = split(str);
            [rows, cols] = size(sp);
            if rows == 1 && cols == 1
                data = str2double(split(sp(1), "R"));
                coe = data(1);
                val = data(2);
                res = sym(eye(n));
                res(val, val) = coe;
                return;
            end
            ero = struct('left', zeros(2, 1), 'op', '', 'right', zeros(2, 1), 'size', n);
            
            ero.left = str2double(split(string(sp(1)), "R"));
            ero.op = char(sp(2)); % This converts the "-" from a string to a character
            ero.right = str2double(split(string(sp(3)), "R"));
            
            switch ero.op
                case '+'
                    res = obj.generateElemAdd(ero);
                case '-'
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
                command = input("Enter the elementary row operation: ", 's'); % 's' forces string input
                if command == "quit"
                    res = M;
                    return;
                end
                
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
        
        function [u,r] = gramSchmidt(obj, v, showSteps)
            arguments
                obj;
                v;
                showSteps logical = true; 
            end
            
            [M, N] = size(v);
            
            % FIX 1: Force the container to be symbolic so it can't convert to decimals
            u = sym(zeros(M, N)); 
            r = sym(eye(N)); 
            
            u(:,1) = v(:,1);
            
            if showSteps
                fprintf("\n--- Gram-Schmidt Step-by-Step ---\n");
                fprintf("Step 1: u_1 is simply the first vector v_1\n");
                disp(u(:,1));
            end
            
            for n = 2:N
                acu = sym(zeros(M, 1)); 
                
                if showSteps
                    fprintf("\nStep %d: Calculating u_%d\n", n, n);
                    fprintf("Formula: u_%d = v_%d", n, n);
                    for m = 1:n-1
                        fprintf(" - proj_u%d(v_%d)", m, n);
                    end
                    fprintf("\n");
                end
                
                for m = 1:n-1
                    num = dot(v(:,n),u(:,m)); 
                    den = dot(u(:,m),u(:,m)); 
                    r(n,m) = num/den;
                    acu = acu + r(n,m)*u(:,m);
                    
                    if showSteps
                        fprintf("  > Projection of v_%d onto u_%d:\n", n, m);
                        % FIX 2: Use string() instead of char() to stop ASCII conversion
                        fprintf("    Inner product <v_%d, u_%d> = %s\n", n, m, string(num));
                        fprintf("    Inner product <u_%d, u_%d> = %s\n", m, m, string(den));
                        fprintf("    Projection vector:\n");
                        disp(r(n,m)*u(:,m));
                    end
                end
                
             u(:,n) = v(:,n) - acu;
             
             if showSteps
                 fprintf("  > Resulting orthogonal vector u_%d:\n", n);
                 disp(u(:,n));
             end
            end
            
            if showSteps
                fprintf("--- End of Process ---\n\n");
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
        
        function res = calcLSS(~, A, b)
            [rowsA, ~] = size(A);
            [rowsB, ~] = size(b);
            if rowsA ~= rowsB
                fprintf("Error: Inputs need to have the same number of rows!\n");
                return;
            end
            % Calculate least squares solution of Ax = b
            left = A' * A;
            right = A' * b;
            R = rref([left right]);
            [~, cols] = size(R);
            v = R(:, cols);
            p = null(R(:, 1:cols - 1), "r");
            
            param = 's';
            fprintf("The general solution has parameters ");
            [~, numP] = size(p);
            for i = 1:numP
                fprintf("%c, ", param);
                param = char(param + 1);
            end
            fprintf("from columns %i to %i", 2, 1 + numP);
           
            res = [v p];
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
        
        function v1 = getEigenvector(obj, A, lambda, output)
            arguments
                obj;
                A; % FIXED: Removed 'double' to allow symbolic matrices
                lambda; % FIXED: Removed 'double' to allow symbolic eigenvalues
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
            
            M = lambda * sym(eye(cols)) - sym(A); % FIXED: Added sym() for exact tracking
            v1 = null(M);
            if output
                fprintf("The matrix\n");
                disp(M);
                fprintf("is reduced to\n")
                disp(rref(M));
                fprintf("Hence the basis for the eigenspace associated with %i is", lambda);
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
        function s = generateInitialConditions(~, n)
            conds = zeros(n, 2);
            for i = 1:n
                fprintf("Enter the t value for y%i: ", i);
                conds(i, 1) = input("");
                fprintf("Enter the result of y%i(%i): ", i, conds(i, 1));
                conds(i, 2) = input("");
            end
            
            s = "[";
            for i = 1:n
                s = s + "y" + i + "(" + conds(i, 1) + ")==" + conds(i, 2);
                if i ~= n
                    s = s + ", ";
                end
            end
            s = s + "]";
        end
        
        function res = solveDifferentialSystem(obj, A, isInitial)
            [rows, ~] = size(A);
            switch rows
                case 2
                    syms y1(t) y2(t);
                    y = [y1; y2];
                case 3
                    syms y1(t) y2(t) y3(t);
                    y = [y1; y2; y3];
                case 4
                    syms y1(t) y2(t) y3(t) y4(t);
                    y = [y1; y2; y3; y4];
                case 5
                    syms y1(t) y2(t) y3(t) y4(t) y5(t);
                    y = [y1; y2; y3; y4; y5];
                otherwise
                    fprintf("Unsupported size.\n");
                    return;
            end
            if ~isInitial
                res = dsolve(diff(y,t) == A * y);
            else
                conds = eval(obj.generateInitialConditions(rows));
                res = dsolve(diff(y,t) == A * y, conds);
            end
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
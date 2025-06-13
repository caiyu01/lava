classdef (InferiorClasses = {?hpf}) lava
% LAserre VAriables objects:
% a lava object represents a matrix of polynomials, M(i,j)
% where the powers of the monomials are stored in opVar(i,j,:,:), [1 1] for x1^2 etc
% then the monomials are summed vertically in opVar, with coefficients
% specified in coeff(i,j,:)
% dim 1: matrix dim 1
% dim 2: matrix dim 2
% dim 3: addition
% dim 4: multiplication
% 
% Some ways to declare a lava object:
% x = lava; % empty
% x = lava([1 2;3 4]); % specifies the variables, but coeff=1
% x = lava({[1 2] 2; 3 4},[1 1;-1 1i]); % specifies the variables and coefficients
% x = lava.num2lava([1 2;3 4]); % a numeric matrix 
% written by Cai Yu, 2020-02-25

% last updated 2020-3-10, new version
% get rid of cells, use 4-D (3-D) array representation for opVar (coeff)
% updated the function accordingly

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           TO DO                              %
%                                                              %
% - dependence: severalBases2dec, uniquecell
% - from double to lava: a = rand(2); a = lava(zeros(size(a)),a);
% - speed up times, mtimes and simplify
% - check speed for assignSdpVar, assign hermitian matrices
% - ambiguity in constructing lava from two matrices
%   when some of the dimension is 1...
% - suport multidimensional arrays
%                                                              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


    properties
% op contains an m by n cell of row vectors (of integers?), opVar
% and a cell of the same size denoting the corresponding coefficients
% for commutative polynomial optimization
        opVar % a four dimensional array of integers
        coeff % a three dimensional array of doubles
    end
    
    methods
        % constructor
        % formally this is a static method
        function opOut = lava(varargin)
        % constructs a lava object
        %
        % To construct a constant lava (i.e. with no variables), use the
        % method lava.num2lava .
        %
        % Examples:
        %     a = lava             % empty
        %     a = lava([1 2;3 4])  % specifies the variables, but coeff=1
        %     a = lava([1 2;3 4], [1 1;-1 1i])  % specifies the variables and coefficients
        %     a = lava({[1 5] 2;3 4}, {1 1;-1 1i})  % including a monomial of degree 2
        %     a = lava({[1; 5] 2;3 4}, {[1; 2] 1;-1 1i})  % including a polynomial term
        %     a = lava(1)                    % x1
        %     a = lava.num2lava(1)           % constant 1
        %     a = lava.num2lava([1 2; 3 4])  % constant matrix [1 2; 3 4]
            switch nargin
                case 0
                    % if coeff1 is not specified
                    opOut = lava.num2lava([]); % size is zero
                case 1
                    if ~iscell(varargin{1})
                        opVar = num2cell(varargin{1});
                    else
                        opVar = varargin{1};
                    end
                    % if coeff1 not specified, assume coeff1 = ones
                    [m1,n1] = size(varargin{1});
                    coeff = cell(m1,n1);
                    for ii=1:numel(coeff)
                        coeff{ii} = ones(size(opVar{ii},1),1);
                    end
                    opOut = lava(opVar, coeff);
                case 2
                    % the main constructor
                    % constructs a m by n lava object
                    % if it is a cell input
                    opVar = varargin{1}; coeff = varargin{2};
                    if iscell(opVar) && iscell(coeff)
                        [m1,n1] = size(varargin{1});
                        [m2,n2] = size(varargin{2});
                        if m1~=m2 || n1~=n2
                            error('Wrong size. opVar and coeff must have same size.')
                        end
                        % wrap it into a 4-matrix
                        % homogenize every cell element to have the maximum width and depth
                        for ii=1:numel(opVar)
                            if (size(opVar{ii},1)~=size(coeff{ii},1) || size(coeff{ii},2)~=1)... % if size does not match
                                    && ~isempty(opVar{ii}) && ~isempty(coeff{ii}) % and it is not empty
                                error('Wrong size. coeff must be an n-by-1 vector, where n is the height of opVar')
                            end
                        end
                        maxDepth = 0; % addition: number of terms in the polynomial
                        maxWidth = 0; % multiplication: maximum degree of monomials
                        if ~isempty(opVar)
                            maxDepth = max(cellfun(@(x) size(x,1), opVar),[],'all');
                            maxWidth = max(cellfun(@(x) size(x,2), opVar),[],'all');
                        end
                        for jj=1:size(opVar,2)
                            for ii=1:size(opVar,1)
                                currDepth = size(opVar{ii,jj},1);
                                currWidth = size(opVar{ii,jj},2);
                                if currDepth<maxDepth || currWidth<maxWidth
                                    opVar{ii,jj} = sort([zeros(maxDepth-currDepth,maxWidth);...
                                                          zeros(currDepth,maxWidth-currWidth), opVar{ii,jj}],2); % commute
                                    coeff{ii,jj} = [zeros(maxDepth-currDepth,1);coeff{ii,jj}];
                                end
                            end
                        end
                        opVar = permute(reshape(cell2mat(opVar),maxDepth,m1,maxWidth,n1),[2,4,1,3]);
                        % sort, since multiplication commutes:
                        opVar = reshape(sort(reshape(opVar,maxDepth*m1*n1,maxWidth),2),m1,n1,maxDepth,maxWidth);
                        % sort, since addition commutes:

                        coeff = permute(reshape(cell2mat(coeff),maxDepth,m1,n1),[2,3,1]);
                        opOut.opVar = opVar;
                        opOut.coeff = coeff;
                    elseif isnumeric(opVar) && isnumeric(coeff)
                        % should allow construction from 4-D opVar with 3-D coeff
                        if ismatrix(opVar) && ismatrix(coeff)
                            % both are 2-D matrix, 
                            % interpreted depth=1 and width = 1
                            if size(opVar,1)~=size(coeff,1) || size(opVar,2)~=size(coeff,2)
                                error('Wrong size. opVar and coeff must have same size.')
                            end
                            opOut.opVar = opVar;
                            opOut.coeff = coeff;
                        elseif ~ismatrix(opVar) % && ~ismatrix(coeff)
                            % in case some dimensions has length 1
                            if (size(opVar,1)==1 || size(opVar,2)==1) && (numel(coeff) > 0)
                                coeff = reshape(coeff,size(opVar,1),size(opVar,2),numel(coeff)/(size(opVar,1)*size(opVar,2)));
                            end
                            if size(opVar,1)~=size(coeff,1) || size(opVar,2)~=size(coeff,2)
                                error('Wrong size. opVar and coeff must have same size.')
                            end
                            % Failed example: lava(ceil(rand(1,3,1,2)), [1 2 3])
                            coeff = reshape(coeff,size(opVar,1),size(opVar,2),size(opVar,3));
                            if size(opVar,3)~=size(coeff,3)
                                error('Wrong size. Depth of opVar and coeff must be the same')
                            end
                            opOut.opVar = opVar;
                            opOut.coeff = coeff;
                        else
                            error(['Ambigous input. Use either 4-D opVar with 3-D coeff'...
                                newline 'or cell input.']);
                        end
                    else
                        error('Invalid input. Input should have the same type.')
                    end                     
                otherwise
                    error('Invalid input.')
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % size and numel
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function varargout = size(op,varargin)
            % return all four dimensions
            % size(op,3) = maxDepth; size(op,4) = maxWidth;
            if nargin==1 && nargout<=1
                % output a vector, if only the matrix dimension
                % not depth and width of opVar
                opVar1 = op.opVar;
                varargout{1} = [size(opVar1,1),size(opVar1,2)];
            else
                [varargout{1:nargout}] = size(op.opVar,varargin{:});
            end
        end
        
        function result = isempty(op)
            result = any(size(op) == 0);
        end
        
        function n = numel(op)
            % only the first two matrix dimensions
            n = size(op,1)*size(op,2);
        end
        
        function l = length(op)
            % need dummy third output to get the correct size
            [m,n,~,~] = size(op);
            if m>=n
                l = m;
            else
                l = n;
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    subscripted referencing and assignment
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Squeezed operators and coefficients access
        function result = sqOpVar(op1)
            result = squeeze(op1.opVar);
        end
        
        function result = sqCoeff(op1)
            result = squeeze(op1.coeff);
        end
        
        % subscripted referencing
        function varargout = subsref(op1,s)
            % Supported indexing
            % a.opVar, a.coeff
            % a.opVar(1), a.opVar(1,1) (not squeezed, use a(1,1).sqOpVar for squeezed opVars!)
            % a(1), a(1,2), a(1,:)
            % a(1).opVar, a(1,2).coeff (not squeezed, 4-D or 3-D)
            
            if length(s) > 1
                % We only need to treat the case a.something(b) in one
                % shot, all other cases, like a.b.c or a(b).c are treated
                % sequentially.
                if ~(isequal(s(1).type, '.') && isequal(s(2).type, '()'))
                    % We apply one step
                    op2 = subsref(op1, s(1));
                    % and then call the function again with the rest
                    if length(s) > 2
                        % not the final call, we only keep one output
                        varargout{1} = subsref(op2, s(2:end));
                    else
                        % final call, we keep all required outputs
                        [varargout{1:nargout}] = subsref(op2, s(2:end));
                    end
                    return;
                end
            end
            
            % Now we can have only one subsref level, or the pattern is
            % something like a.b(c)
            switch s(1).type
              case '.'
                  if length(s) == 1
                      % a.opVar, a.coeff or a.method
                      prop = s(1).subs;
                      [varargout{1:nargout}] = op1.(prop);
                  elseif (length(s) == 2) && isequal(s(2).type, '()')
                      % a.something(b)
                      prop = s(1).subs;
                      arguments = s(2).subs{:};
                      [varargout{1:nargout}] = op1.(prop)(arguments);
                  else
                      error('invalid indexing expression');
                  end
              case '()'
                 % example a(1)
                 [m1, n1, d1, w1] = size(op1);
                 % extract only the matrix dimension
                 t = s;
                 t.subs = [t.subs ':' ':'];
                 if length(s.subs) == 2
                    m2 = length(s.subs{1});
                    n2 = length(s.subs{2});
                    if strcmp(s.subs{1},':')
                        m2 = m1;
                    end
                    if strcmp(s.subs{2},':')
                        n2 = n1;
                    end
                    opVar1 = subsref(op1.opVar,t);
                    opVar1 = reshape(opVar1,[m2,n2,d1,w1]);
                    coeff1 = subsref(op1.coeff,t);
                    coeff1 = reshape(coeff1,[m2,n2,d1]);
                 elseif length(s.subs)==1
                     m2 = length(s.subs{1});
                     if strcmp(s.subs{1},':')
                         m2 = numel(op1);
                     end
                     opVar1 = subsref(reshape(op1.opVar,m1*n1,d1,w1),t);
                     opVar1 = reshape(opVar1,[m2,1,d1,w1]);
                     coeff1 = subsref(reshape(op1.coeff,m1*n1,d1),t);
                     coeff1 = reshape(coeff1,[m2,1,d1]);
                     
                     % Permute if needed
                     if ~strcmp(s.subs{1},':') && (m1 == 1) && (n1 > 1)
                         opVar1 = permute(opVar1, [2 1 3 4]);
                         coeff1 = permute(coeff1, [2 1 3]);
                     end
                 else
                     error('Indexing not supported. a(1) or a(1,2)')
                 end                     
                 varargout{1} = lava(opVar1,coeff1);
               otherwise
                 error('Invalid indexing expression.')
           end
        end
        
        function obj = subsasgn(obj,s,varargin)
            error('subsasgn not supported yet')
%            % Allow subscripted assignment to uninitialized variable
%            if isequal(obj,[])
%               % obj = ClassName.empty;
%            end
% 
%            switch s(1).type
%               case '.'
%                  if length(s) == 1
%                     % Implement obj.PropertyName = varargin{:};
%                     % example: a.coeff = [1 2 3;4 5 6];
%                     prop = s.subs;
%                     if strcmp(prop,'opVar')
%                         obj = lava(varargin{1},obj.coeff);
%                     elseif strcmp(prop,'coeff')
%                         obj = lava(obj.opVar,varargin{1});
%                     else
%                         error('Invalid assignment')
%                     end
%                  elseif length(s) == 2 && strcmp(s(2).type,'()')
%                     % Implement obj.PropertyName(indices) = varargin{:};
%                     % example: a.coeff(2) = 1
%                     prop = s(1).subs;
%                     if strcmp(prop,'opVar')
%                         obj.opVar = subsasgn(obj.opVar,s(2),varargin{1});
%                         obj = lava(obj.opVar,obj.coeff);
%                     elseif strcmp(prop,'coeff')
%                         obj.coeff = subsasgn(obj.coeff,s(2),varargin{1});
%                         obj = lava(obj.opVar,obj.coeff);
%                     else
%                         error('Invalid assignment')
%                     end
%                  else
%                     % Call built-in for any other case
%                     obj = builtin('subsasgn',obj,s,varargin{:});
%                  end
%               case '()'
%                  if length(s) == 1
%                      if ~isa(varargin{1},'lava')
%                          arginClass = class(varargin{1});
%                          error(['cannot assign ' arginClass ' to lava array'])
%                      end
%                     % example a(1,1:2) = lava([1 2]);
%                     obj.opVar = subsasgn(obj.opVar,s,varargin{1}.opVar);
%                     obj.coeff = subsasgn(obj.coeff,s,varargin{1}.coeff);
%                  elseif length(s) == 2 && strcmp(s(2).type,'.')
%                     % example a(1).coeff = {2};
%                     prop = s(2).subs;
%                     if strcmp(prop,'opVar')
%                         obj.opVar = subsasgn(obj.opVar,s(1),varargin{1});
%                         obj = lava(obj.opVar,obj.coeff);
%                     elseif strcmp(prop,'coeff')
%                         obj.coeff = subsasgn(obj.coeff,s(1),varargin{1});
%                         obj = lava(obj.opVar,obj.coeff);
%                     else
%                         error('Invalid assignment')
%                     end
%                  else
%                     % Use built-in for any other expression
%                     obj = builtin('subsasgn',obj,s,varargin{:});
%                  end       
%               otherwise
%                  error('Not a valid indexing expression')
%            end
        end
        
        % changes the default number of output argument?
        % see: https://ch.mathworks.com/matlabcentral/answers/303067-getting-error-output-argument-varargout-2-and-maybe-others-not-assigned-during-call-to-myclas
        % and also: https://ch.mathworks.com/help/matlab/ref/numargumentsfromsubscript.html
        function n = numArgumentsFromSubscript(obj,s,indexingContext)
            n = 1;
        end
        
        function result = end(this, K, N)
            % end - returns the last index
            
            if N == 1
                % Then the indexing is done with only one index, as in a(end)
                result = prod(size(this));
            else
                % Then the indexing is done with two indices, as in a(end,1:2)
                if (K < 1) || (K > 2)
                    error('lava objects only have two dimensions');
                end
                s = size(this);
                result = s(K);
            end
        end
        
        % real part
        function opOut = real(op1)
            opOut = lava(op1.opVar, real(op1.coeff));
            opOut = simplify(opOut, true);
        end
        
        % imaginary part
        function opOut = imag(op1)
            opOut = lava(op1.opVar, imag(op1.coeff));
            opOut = simplify(opOut, true);
        end
        
%         % constant part
%         function opOut = constPart(op1)
%             op1 = simplify(op1); % put all constants together
%             
%             % TODO
%         end
    
        % diagonal part
        function opOut = diag(op1)
            
            [m, n, d, w] = size(op1);
            if m ~= n
                error('Matrix is not square');
            end
            
            opVar = reshape(op1.opVar, [m*n, d*w]);
            coeff = reshape(op1.coeff, [m*n, d]);

            indices = 1:m+1:m^2;
            opVar = reshape(opVar(indices,:), [m 1 d w]);
            coeff = reshape(coeff(indices,:), [m 1 d]);
            
            opOut = simplify(lava(opVar, coeff), true);
        end


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    basic arithmetic operations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function opOut = plus(op1,op2)
            % addition
            
            % Make sure both object are lava objects
            if ~isa(op1, 'lava')
                op1 = lava.num2lava(op1);
            end
            if ~isa(op2, 'lava')
                op2 = lava.num2lava(op2);
            end
        
            [m1, n1, d1, w1] = size(op1);
            [m2, n2, d2, w2] = size(op2);
            % in case one of them is empty
            if isempty(op1)
                opOut = op2;
                return;
            elseif isempty(op2)
                opOut = op1;
                return;
            elseif (max(m1,n1)==1) && ((m1 ~= m2) || (n1 ~= n2))
                % scalar with matrix
                op1 = kron(op1,ones(size(op2)));
                [m1, n1, d1, w1] = size(op1);
            elseif (max(m2,n2)==1) && ((m1 ~= m2) || (n1 ~= n2))
                % matrix with scalar
                op2 = kron(op2,ones(size(op1)));
                [m2, n2, d2, w2] = size(op2);
            end
            if m1~=m2 || n1~=n2
                error('wrong size');
            end
            
            % check maximum width; pad with zeros if not match
            if w1<w2
                wDiff = w2-w1;
                op1.opVar = reshape([zeros(m1*n1*d1,wDiff) reshape(op1.opVar,m1*n1*d1,w1)],m1,n1,d1,w2);
            elseif w2<w1
                wDiff = w1-w2;
                op2.opVar = reshape([zeros(m2*n2*d2,wDiff) reshape(op2.opVar,m2*n2*d2,w2)],m2,n2,d2,w1);
            end
            
            % perform addition
            opVar3 = cat(3,op1.opVar,op2.opVar);
            coeff3 = cat(3,op1.coeff,op2.coeff);
            
            opOut = lava(opVar3,coeff3);
        end
        
        % times
        function opOut = times(op1,op2)   
            % both lava
            if isa(op1,'lava') && isa(op2,'lava')
                
                [m1, n1, d1, w1] = size(op1);
                [m2, n2, d2, w2] = size(op2);
                if (max(m1,n1)==1) && ((m1 ~= m2) || (n1 ~= n2))
                    op1 = kron(op1,ones(size(op2)));
                    [m1, n1, d1, w1] = size(op1);
                elseif (max(m2,n2)==1) && ((m1 ~= m2) || (n1 ~= n2))
                    op2 = kron(op2,ones(size(op1)));
                    [m2, n2, d2, w2] = size(op2);
                elseif m1~=m2 || n1~=n2
                    error('wrong size');
                end
                
                opVar1 = permute(op1.opVar,[3 1 2 4]);
                coeff1 = permute(op1.coeff,[3 1 2]);
                opVar2 = permute(op2.opVar,[3 1 2 4]);
                coeff2 = permute(op2.coeff,[3 1 2]);
                
                opVar1 = kron(reshape(opVar1,m1*n1*d1,w1),ones(d2,1));
                coeff1 = kron(reshape(coeff1,m1*n1*d1,1),ones(d2,1));
                t = kron(ones(d1,1),reshape(opVar2,m2*n2*d2,w2));
                s = kron(ones(d1,1),reshape(coeff2,m2*n2*d2,1));
                opVar2 = reshape(permute(reshape(t,[d2,m2*n2,d1,w2]), [1 3 2 4]), numel(t)/w2,w2);
                coeff2 = reshape(permute(reshape(s,[d2,m2*n2,d1]),[1 3 2]), numel(s),1);
                % sort! commute!
                opVar3 = permute(reshape(sort([opVar1, opVar2],2),d1*d2,m1,n1,w1+w2),[2,3,1,4]);
                coeff3 = permute(reshape(coeff1.*coeff2, d1*d2,m1,n1),[2 3 1]);
                opOut = lava(opVar3,coeff3);
                
            % one of them is double
            % simply .* the coeffs with the double and keep opVar
            elseif isnumeric(op1) && isa(op2,'lava')
                opOut = lava(op2.opVar, op1.*op2.coeff);
            elseif isa(op1,'lava') && isnumeric(op2)
                opOut = lava(op1.opVar, op2.*op1.coeff);
            end
            
            % Avoid too big representations
            opOut = simplify(opOut, true);
        end
        
        % minus
        function opOut = minus(op1,op2)
            opOut = op1 + -1.*op2;
        end
        
        % uminus
        function opOut = uminus(op1)
            opOut = -1.*op1;
        end
        
        % mtimes
        function opOut = mtimes(op1,op2)  
            % matrix multiplication between
            % a double matrix with a lava matrix
            % a lava matrix with a double matrix
            % or two lava matrices
            if isnumeric(op1) && isa(op2,'lava')
                op1 = lava.num2lava(op1);
                opOut = mtimes(op1,op2);
            elseif isa(op1,'lava') && isnumeric(op2)
                op2 = lava.num2lava(op2);
                opOut = mtimes(op1,op2);
            elseif isa(op1,'lava') && isa(op2,'lava')
                % both are lava objects
                [m1,n1,d1,w1] = size(op1);
                [m2,n2,d2,w2] = size(op2);
                if max(m1,n1)==1 || max(m2,n2)==1
                    % if one of them is a scalar
                     opOut = op1.*op2;
                elseif n1~=m2
                    error('inner dimension must match')
                else
                    % load opVar and coeff
                    opVar1 = op1.opVar; opVar2 = op2.opVar;
                    coeff1 = op1.coeff; coeff2 = op2.coeff;
                    % kronecker, match d1 with d2, and m1 with n2
                    opVar1 = kron(ones(n2*d2,1),reshape(opVar1,numel(opVar1),1));
                    coeff1 = kron(ones(n2*d2,1),reshape(coeff1,numel(coeff1),1));
                    opVar2 = kron(ones(m1*d1,1),reshape(opVar2,numel(opVar2),1));
                    coeff2 = kron(ones(m1*d1,1),reshape(coeff2,numel(coeff2),1));
                    % correct the order, n1=m2
                    opVar1 = reshape(permute(reshape(opVar1,m1,n1,d1,w1,n2,d2),[6 3 2 1 5 4]),numel(opVar1)/w1,w1);
                    coeff1 = reshape(permute(reshape(coeff1,m1,n1,d1,n2,d2),[5 3 2 1 4]),numel(coeff1),1);
                    opVar2 = reshape(permute(reshape(opVar2,m2,n2,d2,w2,m1,d1),[3 6 1 5 2 4]),numel(opVar2)/w2,w2);
                    coeff2 = reshape(permute(reshape(coeff2,m2,n2,d2,m1,d1),[3 5 1 4 2]),numel(coeff2),1);
                    % put it back in a lava object
                    % sort since commuting
                    opVar3 = permute(reshape(sort([opVar1 opVar2],2),(d1*d2)*n1,m1,n2,w1+w2),[2 3 1 4]);
                    coeff3 = permute(reshape(coeff1.*coeff2, (d1*d2)*n1,m1,n2),[2 3 1]);
                    opOut = lava(opVar3,coeff3);
                end                
            end
            
            % Avoid too big representations
            opOut = simplify(opOut, true);
        end
        
        % rdivide
        function opOut = rdivide(op1,op2)
            % Element-wise division. We only support specific cases.
            
            if ~isa(op2, 'lava') || isequal(op2.uniqueVar.opVar, 0)
                % division by scalars
                if isa(op2, 'lava')
                    opOut = lava(op1.opVar, op1.coeff./op2.coeff);
                else
                    opOut = lava(op1.opVar, op1.coeff./op2);
                end
            else
                error('Unsupported arguments');
            end
        end
        
        % mrdivide
        function opOut = mrdivide(op1,op2)
            % Matrix division, only division by a scalar is supported
            
            if numel(op2) == 1
                opOut = rdivide(op1,op2);
            else
                error('Unsupported arguments');
            end
        end
        
        % element-wise power
        function opOut = power(op1,exponent)
            if ~isnumeric(exponent) || (numel(exponent) ~= 1) || ~isequal(round(exponent), exponent) || (exponent < 0)
                error('Elemement-wise power is only supported with simple exponents.');
            end
            
            opOut = op1;
            for i = 2:exponent
                opOut = opOut.*op1;
            end
        end
        
        % matrix power
        function opOut = mpower(op1,exponent)
            % We only support matrix powers ... for scalars ;-)
            if numel(op1) ~= 1
                error('Operation only implemented on scalars');
            end
            
            opOut = power(op1,exponent);
        end
        
        % kronecker
        function opOut = kron(varargin)
            % kronecker product between
            % a double matrix with a lava matrix
            % a lava matrix with a double matrix
            % or two lava matrices
            % or more than two objects
            
            if nargin < 2
                error('Not enough arguments');
            elseif nargin > 2
                opOut = varargin{1};
                for i = 2:nargin
                    opOut = kron(opOut, varargin{i});
                end
                return;
            end
            
            assert(nargin == 2);
            % If we are here we have exactly two arguments
            op1 = varargin{1};
            op2 = varargin{2};
            
            if isnumeric(op1) && isa(op2,'lava')
                op1 = lava.num2lava(op1);
                opOut = kron(op1,op2);
            elseif isa(op1,'lava') && isnumeric(op2)
                op2 = lava.num2lava(op2);
                opOut = kron(op1,op2);
            elseif isa(op1,'lava') && isa(op2,'lava')
                [m1, n1, d1, w1] = size(op1);
                [m2, n2, d2, w2] = size(op2);
                % load opVar and coeff
                opVar1 = op1.opVar; opVar2 = op2.opVar;
                coeff1 = op1.coeff; coeff2 = op2.coeff;
                % kronecker, match d1 with d2, and m1 with n2
                opVar1 = kron(ones(m2*n2*d2,1),reshape(opVar1,numel(opVar1),1));
                coeff1 = kron(ones(m2*n2*d2,1),reshape(coeff1,numel(coeff1),1));
                opVar2 = kron(ones(m1*n1*d1,1),reshape(opVar2,numel(opVar2),1));
                coeff2 = kron(ones(m1*n1*d1,1),reshape(coeff2,numel(coeff2),1));
                % correct the order, d2 d1 m2 m1 n2 n1 (w1,w2) 
                opVar1 = reshape(permute(reshape(opVar1,m1,n1,d1,w1,m2,n2,d2),[7 3 5 1 6 2 4]),numel(opVar1)/w1,w1);
                coeff1 = reshape(permute(reshape(coeff1,m1,n1,d1,m2,n2,d2),[6 3 4 1 5 2]),numel(coeff1),1);
                opVar2 = reshape(permute(reshape(opVar2,m2,n2,d2,w2,m1,n1,d1),[3 7 1 5 2 6 4]),numel(opVar2)/w2,w2);
                coeff2 = reshape(permute(reshape(coeff2,m2,n2,d2,m1,n1,d1),[3 6 1 4 2 5]),numel(coeff2),1);
                % put it back in a lava object
                % sort since commuting
                opVar3 = permute(reshape(sort([opVar1 opVar2],2),d1*d2,m1*m2,n1*n2,w1+w2),[2 3 1 4]);
                coeff3 = permute(reshape(coeff1.*coeff2,d1*d2,m1*m2,n1*n2),[2 3 1]);
                opOut = lava(opVar3,coeff3);
            end
            
            % Avoid too big representations
            opOut = simplify(opOut, true);
        end
 
        % conjugate transposition
        function opOut = ctranspose(op1)
            opVar1 = op1.opVar;
            coeff1 = op1.coeff;
            opVar1 = permute(opVar1,[2 1 3 4]);
            coeff1 = conj(permute(coeff1,[2 1 3 4])); % hermitian conjugation
            opOut = lava(opVar1,coeff1);
        end
        
        % transposition
        function opOut = transpose(op1)
         opVar1 = op1.opVar;
            coeff1 = op1.coeff;
            opVar1 = permute(opVar1,[2 1 3 4]);
            coeff1 = permute(coeff1,[2 1 3 4]); % only transposition
            opOut = lava(opVar1,coeff1);
        end
        
        % trace
        function opOut = trace(op1)
            if size(op1,1) ~= size(op1,2)
                error('Matrix must be square.');
            end
            
            opOut = 0;
            for i = 1:size(op1,1)
                vect = [zeros(1,i-1) 1 zeros(1, size(op1,1)-i)]';
                opOut = opOut + vect'*op1*vect;
            end
            
            % Avoid too big representations
            opOut = simplify(opOut, true);
        end
        
        % partial transposition
        function opOut = Tx(op1, sys, dim)
            assert(isa(op1, 'lava'));
            assert(min(sys) >= 1);
            assert(max(sys) <= length(dim));
            assert(isequal(prod(dim)*[1 1], size(op1)));
            
            opVar1 = op1.opVar;
            coeff1 = op1.coeff;
            [m1,n1,d1,w1] = size(opVar1);
            
            nbSys = length(dim);
            permutation = (1:2*nbSys+2);
            permutation([nbSys+1-sys, 2*nbSys+1-sys]) = permutation([2*nbSys+1-sys,nbSys+1-sys]);
            opVar1 = reshape(permute(reshape(opVar1, [dim(end:-1:1), dim(end:-1:1), d1, w1]), permutation), [m1,n1,d1,w1]);
            coeff1 = reshape(permute(reshape(coeff1, [dim(end:-1:1), dim(end:-1:1), d1]), permutation(1:end-1)), [m1,n1,d1]);
            
            opOut = lava(opVar1, coeff1);
        end
        
        % partial trace
        function opOut = TrX(op1, sys, dim)
            error('TODO');
        end
        
        % sum
        function opOut = sum(op1, dim)
            %   sum(a) : normal sum (column-wise for matrices)
            %   sum(a, dim) : sum along the dimension dim
            %   sum(a, 'all') : sum all elements
            if nargin < 2
                if size(op1,1) ~= 1
                    dim = 1;
                else
                    dim = 2;
                end
            end

            % Now we call the relevant sum procedure. Since the function creates a
            % new object with the result, we keep the corresponding handle...
            switch dim
                case 'all'
                    op1 = op1.subsref(struct('type','()','subs',{{':'}}));
                    opOut = ones(1,numel(op1))*op1;
                case 1
                    opOut = ones(1,size(op1,1))*op1;
                case 2
                    opOut = op1*ones(size(op1,2),1);
                otherwise
                    error('Unexpected argument in sum');
            end
            
            % Avoid too big representations
            opOut = simplify(opOut, true);
        end
        
        % prod
        function opOut = prod(op1, dim)
            %   prod(a) : normal product (column-wise for matrices)
            %   prod(a, dim) : product along the dimension dim
            %   prod(a, 'all') : product of all elements
            if nargin < 2
                if size(op1,1) ~= 1
                    dim = 1;
                else
                    dim = 2;
                end
            end

            % Now we call the relevant sum procedure. Since the function creates a
            % new object with the result, we keep the corresponding handle...
            switch dim
                case 'all'
                    op1 = op1.subsref(struct('type','()','subs',{{':'}}));
                    opOut = prod(op1);
                case 1
                    if size(op1,1) > 1
                        op11 = op1.subsref(struct('type','()','subs',{{1,':'}}));
                        op12 = op1.subsref(struct('type','()','subs',{{2,':'}}));
                        op1rest = op1.subsref(struct('type','()','subs',{{3:size(op1,1),':'}}));
                        opOut = prod([op11.*op12; op1rest], 1);
                    else
                        % Avoid too big representations
                        opOut = simplify(op1, true);
                    end
                case 2
                    if size(op1,2) > 1
                        op11 = op1.subsref(struct('type','()','subs',{{':',1}}));
                        op12 = op1.subsref(struct('type','()','subs',{{':',2}}));
                        op1rest = op1.subsref(struct('type','()','subs',{{':',3:size(op1,2)}}));
                        opOut = prod([op11.*op12, op1rest], 2);
                    else
                        % Avoid too big representations
                        opOut = simplify(op1, true);
                    end
                otherwise
                    error('Unexpected argument in prod');
            end
        end
        
        % round
        function opOut = round(op1, nbDigits)
            % Rounding off coefficients
            %
            % round(op1) rounds the coefficients of op1 to the nearest
            %   integer
            % round(op1, nbDigits) rounds the coefficients to the nbDigits
            %   most significant digits
            %
            % See also: round
            
            opOut = op1;
            if nargin == 1
                opOut.coeff = round(opOut.coeff);
            else
                opOut.coeff = round(opOut.coeff, nbDigits);
            end
            opOut = simplify(opOut, true);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    basic tests
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % test if object is real
        function result = isreal(op1)
            result = isreal(op1.coeff);
        end
        
        % test if the object is a scalar
        function result = isscalar(op1)
            % Returns true iff op1 is 1x1 and a constant (i.e. involves no
            % variable).
            
            % First we check the dimension
            if numel(op1) ~= 1
                result = false;
                return;
            end
            
            % Now we have exactly one element. We check if it involves any
            % variables.
            op1 = simplify(op1);
            result = (op1.opVar == 0);
        end
        
        % test of symmetry
        function result = issymmetric(op1)
            % First we simplify the object
            op1 = simplify(op1);
            
            % Now we check the structure
            result = isequal(op1, op1.transpose);
        end
        
        % test of hermiticity
        function result = ishermitian(op1)
            % First we simplify the object
            op1 = simplify(op1);
            op2 = simplify(op1.ctranspose); % for matlab, -1i <= 1 <= 1i, so simplifying again is crucial here
            
            % Now we check the structure
            result = isequal(op1, op2);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    basic matrix manipulations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % reshape
        function result = reshape(op1, varargin)
            if (length(varargin)==1) && (length(varargin{1}) == 2)
                s1 = varargin{1}(1);
                s2 = varargin{1}(2);
            elseif length(varargin)==2
                s1 = varargin{1};
                s2 = varargin{2};
            else
                error('Invalid shape argument. reshape(a,2,2) or reshape(a,[2 2])')
            end
            opVar1 = op1.opVar;
            coeff1 = op1.coeff;
            opVar1 = reshape(opVar1,s1,s2,size(opVar1,3),size(opVar1,4));
            coeff1 = reshape(coeff1,s1,s2,size(coeff1,3),size(coeff1,4));
            result = lava(opVar1,coeff1);
        end
        
        % match sizes
        function [op1, op2] = matchSize(op1,op2)
            [m1,n1,d1,w1] = size(op1);
            [m2,n2,d2,w2] = size(op2);
            opVar1 = op1.opVar; coeff1 = op1.coeff;
            opVar2 = op2.opVar; coeff2 = op2.coeff;
            if d1<d2
                opVar1 = reshape(permute(reshape(opVar1,m1,n1,d1,w1),[3 1 2 4]),d1,m1*n1*w1);
                opVar1 = [zeros(d2-d1,m1*n1*w1); opVar1];
                opVar1 = reshape(permute(reshape(opVar1,d2,m1,n1, w1),[2 3 1 4]),m1,n1,d2,w1);
                coeff1 = reshape(permute(reshape(coeff1,m1,n1,d1),[3 1 2]),d1,m1*n1);
                coeff1 = [zeros(d2-d1,m1*n1); coeff1];
                coeff1 = reshape(permute(reshape(coeff1,d2,m1,n1),[2 3 1]),m1,n1,d2);
                d1 = d2;
            elseif d2<d1
                opVar2 = reshape(permute(reshape(opVar2,m2,n2,d2,w2),[3 1 2 4]),d2,m2*n2*w2);
                opVar2 = [zeros(d1-d2,m2*n2*w2); opVar2];
                opVar2 = reshape(permute(reshape(opVar2,d1,m2,n2, w2),[2 3 1 4]),m2,n2,d1,w2);
                coeff2 = reshape(permute(reshape(coeff2,m2,n2,d2),[3 1 2]),d2,m2*n2);
                coeff2 = [zeros(d1-d2,m2*n2); coeff2];
                coeff2 = reshape(permute(reshape(coeff2,d1,m2,n2),[2 3 1]),m2,n2,d1);
                d2 = d1;
            end
            if w1<w2
                opVar1 = reshape(opVar1,m1*n1*d1,w1);
                opVar1 = [zeros(m1*n1*d1,w2-w1) opVar1];
                opVar1 = reshape(opVar1,m1,n1,d1,w2);
            elseif w2<w1
                opVar2 = reshape(opVar2,m2*n2*d2,w2);
                opVar2 = [zeros(m2*n2*d2,w1-w2) opVar2];
                opVar2 = reshape(opVar2,m2,n2,d2,w1);
            end
            op1 = lava(opVar1,coeff1);
            op2 = lava(opVar2,coeff2);
        end
        
        % concatenation
        function result = horzcat(varargin)
            if nargin>2
                result = horzcat(varargin{1},horzcat(varargin{2:end}));
            elseif nargin==2
                % Make sure both object are lava objects
                op1 = varargin{1};
                op2 = varargin{2};
                if ~isa(op1, 'lava')
                    op1 = lava.num2lava(op1);
                end
                if ~isa(op2, 'lava')
                    op2 = lava.num2lava(op2);
                end

                % 0x0 or 0xd objects are the same
                if (size(op1,1) == 0) && (size(op1,2) == 0) && (size(op2,1) ~= 0)
                    % op1 is empty, matching size with op2
                    op1 = lava(zeros(size(op2,1), 0), zeros(size(op2,1), 0));
                end
                if (size(op2,1) == 0) && (size(op2,2) == 0) && (size(op1,1) ~= 0)
                    % op2 is empty, matching size with op1
                    op2 = lava(zeros(size(op1,2), 0), zeros(size(op1,2), 0));
                end
                
                % need to standardize the width and depth
                [op1, op2] = matchSize(op1,op2);
                s.type = '.';
                s.subs = 'opVar';
                opVar1 = subsref(op1,s);
                opVar2 = subsref(op2,s);    
                s.subs = 'coeff';
                coeff1 = subsref(op1,s);
                coeff2 = subsref(op2,s);
                try
                    opVar3 = [opVar1 opVar2];
                    coeff3 = [coeff1 coeff2];
                    result = lava(opVar3,coeff3);
                catch ME
                    error('Dimensions of arrays being concatenated are not consistent.')
                end
            else
                result = varargin{1};
            end
        end
        
        function result = vertcat(varargin)
            if nargin>2
                result = vertcat(varargin{1},vertcat(varargin{2:end}));
            elseif nargin==2
                % Make sure both object are lava objects
                op1 = varargin{1};
                op2 = varargin{2};
                if ~isa(op1, 'lava')
                    op1 = lava.num2lava(op1);
                end
                if ~isa(op2, 'lava')
                    op2 = lava.num2lava(op2);
                end
                
                % 0x0 or 0xd objects are the same
                if (size(op1,1) == 0) && (size(op1,2) == 0) && (size(op2,2) ~= 0)
                    % op1 is empty, matching size with op2
                    op1 = lava(zeros(0, size(op2,2)), zeros(0, size(op2,2)));
                end
                if (size(op2,1) == 0) && (size(op2,2) == 0) && (size(op1,2) ~= 0)
                    % op2 is empty, matching size with op1
                    op2 = lava(zeros(0, size(op1,2)), zeros(0, size(op1,2)));
                end
                
                % need to standardize the width and depth
                [op1, op2] = matchSize(op1,op2);
                s.type = '.';
                s.subs = 'opVar';
                opVar1 = subsref(op1,s);
                opVar2 = subsref(op2,s);
                s.subs = 'coeff';
                coeff1 = subsref(op1,s);
                coeff2 = subsref(op2,s);
                try
                    opVar3 = [opVar1; opVar2];
                    coeff3 = [coeff1; coeff2];
                    result = lava(opVar3,coeff3);
                catch ME
                    error('Dimensions of arrays being concatenated are not consistent.')
                end
            else
                result = varargin{1};
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    display functions
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % toStr
        function result = toStr(op1)
            % This function returns a string description of the polynomial matrix
            result = cell(size(op1));
            
            for i = 1:size(op1,1)
                for j = 1:size(op1,2)
                    text = ' + ';
                    for k = 1:size(op1.opVar,3)
                        if (op1.coeff(i,j,k) ~= 0)
                            coeffTxt = strrep(strrep(num2str(op1.coeff(i,j,k)), '+', ' + '), '-', ' - ');
                            % Deal with exponential notation
                            coeffTxt = strrep(strrep(coeffTxt, 'e + ', 'e+'), 'e - ', 'e-');
                            if (length(coeffTxt) > 2) && isequal(coeffTxt(1:2), '0 ')
                                % purely imaginary term, remove the zero
                                coeffTxt = coeffTxt(2:end);
                            end
                            if sum(op1.opVar(i,j,k,:)) == 0
                                % just a constant term
                                if (length(coeffTxt) > 2) && (coeffTxt(2) == '-')
                                    text = [text(1:end-3), coeffTxt, ' + '];
                                else
                                    text = [text, coeffTxt, ' + '];
                                end
                            else
                                % some operators are involved
                                if isreal(op1.coeff(i,j,k)) || isreal(1i*op1.coeff(i,j,k))
                                    % either real of imaginary coefficient
                                    % we first adjust the sign
                                    start = 1;
                                    if (length(coeffTxt) > 2) && (coeffTxt(2) == '-')
                                        text(end-1) = '-';
                                        start = 4;
                                    elseif (length(coeffTxt) > 2) && isreal(1i*op1.coeff(i,j,k))
                                        start = 4;
                                    end
                                    if isreal(op1.coeff(i,j,k))
                                        % purely real coefficient
                                        if abs(op1.coeff(i,j,k)) ~= 1
                                            text = [text, coeffTxt(start:end), '*'];
                                        end
                                    else
                                        % purely imaginary coefficient
                                        if abs(op1.coeff(i,j,k)) == 1
                                            % we write just 1i is unity
                                            text = [text, '1i*'];
                                        else
                                            text = [text, coeffTxt(start:end), '*'];
                                        end
                                    end
                                else
                                    % both real and imaginary coefficients
                                    if (length(coeffTxt) > 2) && (coeffTxt(2) == '-')
                                        text = [text, '(- ', coeffTxt(4:end), ')*'];
                                    else
                                        text = [text, '(', coeffTxt, ')*'];
                                    end
                                end
                                
                                ops = [0 squeeze(op1.opVar(i,j,k,:))' 0];
                                operatorStart = 1;
                                for l = 2:length(ops)
                                    if ops(l) ~= ops(l-1)
                                        % found a new operator, we:
                                        % 1. finish writing the previous
                                        %    operator
                                        power = l-operatorStart;
                                        if (ops(l-1) ~= 0)
                                            if (power ~= 1)
                                                text = [text, '^', num2str(power), '*'];
                                            else
                                                text = [text, '*'];
                                            end
                                        end
                                        % 2. write the new operator
                                        if l-1 <= size(op1.opVar,4)
                                            text = [text, 'x', num2str(ops(l))];
                                        end
                                        operatorStart = l;
                                    end
                                end
                                text = [text(1:end-1), ' + '];
                            end
                        end
                    end
                    if length(text) > 6
                        if isequal(text(2), '+')
                            % We don't print the initial plus sign
                            result{i,j} = text(4:end-3);
                        else
                            % Don't print the space before the first minus
                            % sign
                            result{i,j} = text(2:end-3);
                        end
                    else
                        result{i,j} = '0';
                    end
                end
            end
        end
        
        % disp
        function disp(op1,str1)
            if (nargin >= 2) && strcmpi(str1, 'mathematica')
                forMathematica = true;
                str1 = 'full';
            else
                forMathematica = false;
            end
            
            if (nargin >= 2) && strcmpi(str1,'full')
                niceDisplay = true;
                shortenDisplay = false;
            else
                niceDisplay = (size(op1.opVar,2)*size(op1.opVar,3)*size(op1.opVar,4) <= 1000);
                shortenDisplay = niceDisplay && (size(op1,1) > 30);
            end
            
            if forMathematica
                shortenDisplay = false;
            end
            
            if (numel(op1) == 0) || ~niceDisplay
                disp(['  ', num2str(size(op1,1)), 'x', num2str(size(op1,2)), ' lava array (polynomials of degree ', num2str(size(op1.opVar,4)), ' with up to ', num2str(size(op1.opVar,3)), ' terms)']);
            else
                disp(['  ', num2str(size(op1,1)), 'x', num2str(size(op1,2)), ' lava array (polynomials of degree ', num2str(size(op1.opVar,4)), ' with up to ', num2str(size(op1.opVar,3)), ' terms):']);
            end
            disp(' ');

            if niceDisplay
                % possibly display more info
                
                % identify lines to display
                if shortenDisplay
                    LinesToDisplay = [1:24, size(op1,1)-5:size(op1,1)];
                    % Get text description of all elements to be printed
                    cellDescription = op1.subsref(struct('type','()','subs',{{LinesToDisplay,':'}})).toStr;
                else
                    LinesToDisplay = 1:size(op1,1);
                    % Get text description of all elements to be printed
                    cellDescription = op1.toStr;
                end

                % Finds the longest string in each column
                colLen = zeros(1, size(op1,2));
                for i = 1:size(cellDescription,1)
                    for j = 1:size(op1,2)
                        colLen(j) = max(colLen(j), length(cellDescription{i,j}));
                    end
                end
                colLen = colLen + 2;
                
                % print each element
                if forMathematica
                    for i = 1:size(cellDescription,1)
                        if i == 1
                            text = '    {{';
                        else
                            text = '    {';
                        end
                        for j = 1:size(op1,2)
                            text = [text, cellDescription{i,j}];
                            if j < size(op1,2)
                                text = [text, ','];
                            end
                            text = [text, char(kron(' ', ones(1, colLen(j)-length(cellDescription{i,j})-1)))];
                            if j == size(op1,2)
                                if i < size(op1,1)
                                    text = [text, '},'];
                                else
                                    text = [text, '}}'];
                                end
                            end
                        end
                        disp(text);
                    end
                else
                    co = 0;
                    for i = LinesToDisplay
                        co = co + 1;
                        if shortenDisplay && (i == size(op1,1)-5)
                            disp('   (...)');
                        end
                        text = '    ';
                        for j = 1:size(op1,2)
                            text = [text, cellDescription{co,j}];
                            text = [text, char(kron(' ', ones(1, colLen(j)-length(cellDescription{co,j}))))];
                        end
                        disp(text);
                    end
                end
            end
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    Hierarchy relaxation operations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [opOut, selectedVariables] = localize(op1, list)
            % This function returns a matrix that localizes the matrix or
            % polynomial defined in 'op1'. The basis used to define the
            % localizing matrix is formed by the constant '1' and the
            % variables present in 'list'. These variables are ordered and
            % duplicates removed. The second output 'selectedVariables'
            % returns these variables in the order used.
            %
            % Note: The order of the element
            %
            % Example:
            %     x1 = lava(1);
            %     op1 = x1*x1;
            %     list = [x1, x1*x1, x1*x1*x1]
            %     op1.localize(list)
            %
            % See also:
            %     lava.uniqueVar
            
            % Make sure both inputs are lava object
            if ~isa(op1, 'lava')
                op1 = lava.num2lava(op1);
            end
            if ~isa(list, 'lava')
                list = lava.num2lava(list);
            end

            selectedVariables = uniqueVar(lava(0), list);
            base = selectedVariables*selectedVariables';
            if isequal(op1, lava(0))
                % simple case: localizing 1 gives the standard SDP matrix
                opOut = base;
            else
                opOut = kron(base, op1);
            end
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    other operations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % simplify
        % slow: elementwise operation!
        % try arrayfun?
        function opOut = simplify(op1, quick)
            % If the (optional) quick parameter is true, no element-wise
            % operation is performed, only global ones.
            
            if nargin < 2
                quick = false;
            end
            
            [m1, n1, d1, w1] = size(op1);
            opVar1 = op1.opVar;
            coeff1 = op1.coeff;
            
            % First, we simplify multiplications assuming that the multiplication commutes:
            opVar1 = reshape(sort(reshape(opVar1,m1*n1*d1,w1),2),m1,n1,d1,w1);

            % Now we sort the addition dimension, d
            % elementwise
            opVar1 = reshape(permute(reshape(opVar1,m1,n1,d1,w1),[3 4 1 2]),d1,numel(opVar1)/d1);
            coeff1 = reshape(permute(reshape(coeff1,m1,n1,d1), [3 1 2]),d1,numel(coeff1)/d1);
            
            if ~quick
                % element-wise merging of identical terms
                for ii=1:m1*n1
                    tmpV = opVar1(:,(1:w1)+(ii-1)*w1);
                    tmpC = coeff1(:,ii);
                    [uniqV, ~, labels] = unique(tmpV,'rows');
                    if size(uniqV,1)<size(tmpV,1)
                        % identify repeated terms only
                        repeated = find(sum(sparse(1:length(labels), labels, 1)) > 1);
                        for jj = repeated
                            % collapse the coeff into the last one
                            % and set opVar to all 0 after the coeff is removed
                            fidx = find(labels == jj);
                            assert(length(fidx) > 1); % We should have more than one term if we arrive here
                            tmpC(fidx(end)) = sum(tmpC(fidx));
                            tmpC(fidx(1:end-1)) = 0;
                            tmpV(fidx(1:end-1),:) = 0;
                        end
                    end
                    % any variable with a zero coefficient is removed
                    tmpV(tmpC==0,:) = 0;
                    % sort according to opVar
                    [tmpV, idx] = sortrows(tmpV); 
                    tmpC = tmpC(idx); % the coeffs follows
                    % the constant is moved
                    positionOfConstant = find((tmpC~=0).*prod(tmpV==0,2));
                    targetPosition = sum(sum(tmpV,2)==0);
                    if ~isempty(positionOfConstant)
                        tmpC([positionOfConstant, targetPosition]) = tmpC([targetPosition, positionOfConstant]);
                    end
                    % save the result
                    opVar1(:,(1:w1)+(ii-1)*w1) = tmpV;
                    coeff1(:,ii) = tmpC;
                end
            end
            
            % trim the top zeros, if any
            % reduce maximum depth
            idx = prod(coeff1==0,2) == 1;
            if sum(idx)==d1
                % We keep always at least depth of 1
                opVar1 = 0*opVar1(1,:);
                coeff1 = 0*coeff1(1,:);
            else
                opVar1(idx,:) = [];
                coeff1(idx,:) = [];
            end
            % new maxDepth
            d1 = size(opVar1,1);
            
            % trim the left most zeros, if any
            % reduce maximal width
            opVar1 = sort(reshape(permute(reshape(opVar1,d1,w1,n1,m1),[3 4 1 2]),numel(opVar1)/w1,w1),2);
            coeff1 = reshape(permute(reshape(coeff1,d1,n1,m1),[2 3 1]),m1,n1,d1);
            idx = prod(opVar1==0,1) == 1;
            % We always keep at least a width of 1
            if isequal(idx, ones(1,w1))
                idx(1) = 0;
            end
            opVar1(:,idx) = [];
            % new maxWidth
            w1 = size(opVar1,2);
            % put it back
            opOut = lava(reshape(opVar1,m1,n1,d1,w1),coeff1);
        end
        
        function opOut = substitute(op1, op2, value)
            assert((prod(size(op2)) == 1) && isscalar(value), 'Only assignment of a single operator is currently supported');
            
            assert(~isa(value, 'lava'), 'Only scalar assignments is currently supported')
            
            % We make sure lava objects are in their simplified form
            op1 = simplify(op1);
            op2 = simplify(op2);
            assert(size(op2,3) == 1, 'Only monomial assignment is currently available');
            
            [m1, n1, d1, w1] = size(op1);
            opVar1 = op1.opVar;
            coeff1 = op1.coeff;
            [m2, n2, d2, w2] = size(op2);
            opVar2 = op2.opVar;
            coeff2 = op2.coeff;
            
            if w1 < w2
                % Nothing to be done
                opOut = op1;
                return;
            end
            
            % We simply check each element
            for i = 1:m1
                for j = 1:n1
                    for k = 1:d1
                        sel = [];
                        co = 1;
                        for l = 1:w1
                            if opVar1(i,j,k,l) == opVar2(1,1,1,co)
                                sel = [sel l];
                                co = co + 1;
                                if co > w2
                                    co = 1;
                                end
                            end
                        end
                        if length(sel) >= w2
                            % We found the monomial we were looking for at
                            % least once, let's substitute it
                            
                            % We remove partial identifications
                            nbFound = floor(length(sel)/w2);
                            sel = sel(1:w2*nbFound);
                            
                            % We remove the monomial
                            opVar1(i,j,k,sel) = 0;
                            
                            % And adjust the coefficient
                            coeff1(i,j,k) = coeff1(i,j,k)*(value/coeff2(1,1,1)).^nbFound;
                        end
                    end
                end
            end
            
            % Create and slightly simplify the result
            opOut = lava(reshape(opVar1,m1,n1,d1,w1),coeff1);
            simplify(opOut, true);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    basic variable-related operations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % TODO IN THIS SECTION : Ignore variables which come with a
        %                        coefficient 0!
        
        function [out, IA] = unique(varargin)
            % from several lava objects
            % returns the list of unique elements
            % example: a = lava([1 2 3])
            %          out = unique(a,kron(a,a))
            %          [out, IA] = unique(a,kron(a,a))
            
            % we concatenate all objects together
            op = [];
            allElements = struct('type', '()', 'subs', {{':'}});
            for ii = 1:nargin
                assert(isa(varargin{ii},'lava'));
                op = [op; subsref(varargin{ii}, allElements)];
            end
            
            % make sure all elements are in simplified form
            op = simplify(op);
            
            % and finally keep only one copy of each unique element
            s = size(op.opVar);
            if length(s) < 4
                s(end+1:4) = 1;
            end
            table = [reshape(op.opVar, s(1)*s(2), s(3)*s(4)), reshape(op.coeff, s(1)*s(2), s(3))];
            [uniques, IA] = unique(table, 'rows');
            
            % Creates the corresponding lava object
            s2 = s;
            s2(1) = size(uniques,1);
            opVar1 = reshape(uniques(:,1:s2(3)*s2(4)), s2(1), s2(2), s2(3), s2(4));
            coeff1 = reshape(uniques(:,s2(3)*s2(4)+1:end), s2(1)*s2(2), s2(3));
            out = lava(opVar1, coeff1);
            
            % eventually, adjust size
            if (nargin == 1) && (size(varargin{1},1) == 1) && (size(varargin{1},2) > 1)
                out = out.';
            end
        end
        
        function out = uniqueVar(varargin)
            % from several lava objects
            % returns the list of unique variables
            % either in a lava objects (default)
            % or a unique list (n-by-width) array
            % example: a = lava([1 2 3])
            %          out = uniqueVar(a,kron(a,a));
            %          out = uniqueVar(a,kron(a,a),'array');
            
            if strcmp(varargin{end},'array')
                maxWidth = 1;
                list = zeros(0,1);
                for ii=1:nargin-1
                    % Make sure this is a lava object
                    if ~isa(varargin{ii}, 'lava')
                        varargin{ii} = lava.num2lava(varargin{ii});
                    end
                    [~,~,~,w] = size(varargin{ii});
                    opVar1 = varargin{ii}.opVar;
                    tmp = reshape(opVar1,numel(opVar1)/w,w);
                    % We remove the variables which come with a coefficient
                    % of zero
                    tmp = tmp(varargin{ii}.coeff ~= 0,:);
                    if w<maxWidth
                        tmp = [zeros(size(tmp,1),maxWidth-w) tmp];
                    elseif w>maxWidth
                        list = [zeros(size(list,1),w-maxWidth) list];
                        maxWidth = w;
                    end
                    list = unique([list; tmp],'rows');
                end
                out = list;
            else
                % We construct the output from the unique array
                list = uniqueVar(varargin{:},'array');
                out = simplify(lava(mat2cell(list,ones(1,size(list,1)),size(list,2))));
            end
        end
        
        % read out to variable index
        function [maxWidth, maxBase, uniqueIdx, varargout] = assignVarIdx(varargin)
            % [maxWidth, maxBase, uniqueIdx, varIdxCell] = assignVarIdx(lava1, lava2, {lava3, lava4},...)
            
            % in case some argin are cells of lavas
            co=1;
            for ii=1:nargin
                if iscell(varargin{ii})
                    for jj = 1:numel(varargin{ii})
                        input{co} = varargin{ii}{jj};
                        co = co+1;
                    end
                else
                    input{co} = varargin{ii};
                    co = co+1;
                end
            end
            
            % Make sure all inputs are lava object
            for i = 1:length(input)
                if ~isa(input{i}, 'lava')
                    input{i} = lava.num2lava(input{i});
                end
            end
            
            maxWidth = -1;
            maxBase = -1;
            uniqueIdx = [];
            for ii=1:numel(input)
                % find the maximum width and maximum bases needed
                maxWidth = max([maxWidth,size(input{ii},4)]);
                maxBase = max([maxBase, max(input{ii}.opVar,[],'all')]);
            end
            % the largest number needed
            maxBase = maxBase+1;
            
            % compute output
            output = cell(size(input));
            for co = 1:length(input)
                [m, n, d, w] = size(input{co});
                tmp = reshape(input{co}.opVar,[m*n*d, w]);
                vvarIdx = reshape(lava.severalBases2dec(tmp,maxBase),m,n,d);
                output{co}.varIdx = vvarIdx;
                % just keeping a copy of the coefficients, for the
                % corresponding varIdx
                % maybe we don't need it
                output{co}.coeff = input{co}.coeff; 
                % wrap vvarIdx into column vector to make the concatenation
                uniqueIdx = unique([uniqueIdx; vvarIdx(:)]);
            end

            % put the output in the same way as input
            % distinguish the case of lava and cell of lavas
            varIdxCell = cell(1,nargin);
            co = 1;
            for ii=1:nargin
                if iscell(varargin{ii})
                    for jj = 1:numel(varargin{ii})
                        varIdxCell{ii}{jj} = output{co};
                        co = co+1;
                    end
                else
                    varIdxCell{ii} = output{co};
                    co = co+1;
                end
            end
            varargout = varIdxCell;
        end
        
        % extracts the decomposition in terms of monomials
        function [M, monomials] = decompose(op1)
            % [M, monomials] = decompose(op1)
            % 
            % Decomposes all elements of a lava object in terms of its
            % monomials. The output is always a 2-D matrix. For a n x 1
            % object op1, the output satisfies
            %    op1 = M*monomials'
            %
            %
            % Example:
            % 
            % polys = kron(lava([1 2 3]'), lava(round(rand(3,1))))
            % [M, monomials] = decompose(polys)
            
            monomials = uniqueVar(op1).';

            [maxWidth, maxBase, uniqueIdx, result] = assignVarIdx(op1(:));

            tmp = inversePerm(uniqueIdx+1);
            indices = squeeze(result.varIdx)+1;
            coeffs = squeeze(result.coeff);

            M = full(sparse([1:size(indices,1)]'*ones(1,size(indices,2)), tmp(indices), coeffs));
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    Crater bridge operations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        function [varargout] = assignSdpVar(varargin)
            % return sdp matrices, also the sparse variable vSp
            % in the same structure as varargin
            % example: [ma, mb, mC, vSp] = assignSdpVar(a,b,C), where C is
            % a cell of lava objects. Then mC will be a cell of sdp matrices
            
            [~, ~, uniqueIdx, varIdxCell{1:nargin}] = assignVarIdx(varargin{:});
            nbVar = length(uniqueIdx);
            fprintf('  Number of sdpvar needed: %d', nbVar)
            fprintf(newline);fprintf(newline)
            v = sdpvar(1,nbVar);
            
            % sparse assignment
            % need to add 1, because we can't assign vSp(0)
            vSp = sparse(uniqueIdx(:)+1,ones(nbVar,1),v);
            % lava(0) = 1:
            vSp(1) = 1;
            
            varargout = cell(1,nargin);
            for ii=1:nargin
                if iscell(varargin{ii})
                    % if it is a cell
                    for kk = 1:numel(varargin{ii})
                        % operate on each cell element
                        fprintf('Assigning cell %d matrix %d', ii, kk)
                        fprintf(newline)
                        tmpVarIdx = varIdxCell{ii}{kk}.varIdx;
                        tmpCoeff = varIdxCell{ii}{kk}.coeff;
                        
                        [m1,n1,d1,~] = size(varargin{ii}{kk});
                    
                        tmp = reshape(vSp(tmpVarIdx+1),m1,n1,d1);
                        varargout{ii}{kk} = sum(tmp.*tmpCoeff,3);
                    end
                else
                    fprintf('Assigning matrix %d', ii)
                    fprintf(newline)
                    tmpVarIdx = varIdxCell{ii}.varIdx;
                    tmpCoeff = varIdxCell{ii}.coeff;

                    % make the cell into a big matrix
                    % assign sdpvar, only once
                    % then do matrix dot product
                    % finally make it into the original shape
                    [m1,n1,d1,~] = size(varargin{ii});
                    
                    tmp = reshape(vSp(tmpVarIdx+1),m1,n1,d1);
                    varargout{ii} = sum(tmp.*tmpCoeff,3);
                end
            end
            varargout = [varargout {vSp}];
        end
        
        % Equality constraint
        function constraint = eq(op1, op2)
            % Make sure both object are lava objects
            if ~isa(op1, 'lava')
                op1 = lava.num2lava(op1);
            end
            if ~isa(op2, 'lava')
                op2 = lava.num2lava(op2);
            end
            
            constraint = crater(op1, '==', op2);
        end
        
        % <= constraint
        function constraint = le(op1, op2)
            % Make sure both object are lava objects
            if ~isa(op1, 'lava')
                op1 = lava.num2lava(op1);
            end
            if ~isa(op2, 'lava')
                op2 = lava.num2lava(op2);
            end
            
            constraint = crater(op1, '<=', op2);
        end
        
        % >= constraint
        function constraint = ge(op1, op2)
            % Make sure both object are lava objects
            if ~isa(op1, 'lava')
                op1 = lava.num2lava(op1);
            end
            if ~isa(op2, 'lava')
                op2 = lava.num2lava(op2);
            end
            
            constraint = crater(op1, '>=', op2);
        end
        
        % < constraint
        function constraint = lt(op1, op2)
            constraint = le(op1, op2);
        end
        
        % > constraint
        function constraint = gt(op1, op2)
            constraint = ge(op1, op2);
        end
        
     end
     
     methods(Static)
         
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         %    Some more constructors
         %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
         
         function opOut = symmetric(dim, var)
            % Constructs a symmetric real matrix of size dim x dim, starting
            % with variable var. var can also be a vector of variables
            % (with dim*(dim+1)/2 elements). In this case, the first dim
            % elements are assigned to the diagonal
            %
            % Example:
            %   lava.symmetric(4, 1)
            %   lava.symmetric(3, [0 1 2 0 3 0])
            
            % Check input
            if numel(var) > 1
                if (numel(var) ~= dim*(dim+1)/2)
                    error('Var should be a scalar, or a vector of length dim*(dim+1)/2');
                end
                diagIndices = cumsum([1 (dim:-1:2)]);
                vars = zeros(1,dim*(dim+1)/2);
                vars(diagIndices) = var(1:dim);
                vars(setdiff(1:dim*(dim+1)/2,diagIndices)) = var(dim+1:end);
            else
                vars = var + [0:dim*(dim+1)/2-1];
            end
            
            % We identify diagonal and lower-diagonal elements
            d = 1:dim+1:dim^2;
            down = (tril(ones(dim, dim),-1) ~= 0);
            
            M = zeros(dim, dim);
            M(d) = 1;
            M(down) = 1;
            
            M(M~=0) = vars;
            Mdown = M;
            Mdown(d) = 0;
            M = M + Mdown.';
            
            opOut = lava(M);
        end
        
        function opOut = antiSymmetric(dim, var)
            % Constructs an antisymmetric real matrix of size dim x dim, starting
            % with variable var. Var can also be a vector of variables
            % (with dim*(dim-1)/2 elements).
            %
            % Example:
            %   lava.antiSymmetric(4, 1)
            %   lava.antiSymmetric(3, [1 2 3])
            
            % Check input
            if numel(var) > 1
                if (numel(var) ~= dim*(dim-1)/2)
                    error('Var should be a scalar, or a vector of length dim*(dim+1)/2');
                end
                vars = var(:);
            else
                vars = var + [0:dim*(dim-1)/2-1];
            end
            
            % We identify lower-diagonal elements
            down = (tril(ones(dim, dim),-1) ~= 0);
            
            M = zeros(dim, dim);
            M(down) = 1;
            
            M(M~=0) = vars;
            Mdown = M;
            M = M + Mdown.';
            
            % The coefficients are not all the same
            coeffs = zeros(dim, dim);
            coeffs(down) = -1;
            coeffs = coeffs - coeffs.';
            
            opOut = lava(M, coeffs);
        end
        
        function opOut = hermitian(dim, var)
            % Constructs a hermitian matrix of size dim x dim, starting
            % with variable var. var can also be a vector of variables
            % (with dim^2 elements)
            %
            % Example:
            %   lava.hermitian(3, 1)
            
            % Check input
            if numel(var) > 1
                if (numel(var) ~= dim^2)
                    error('Var should be a scalar, or a vector of length dim*(dim+1)/2');
                end
                vars = var(:);
            else
                vars = var + [0:dim^2-1];
            end
            
            % The real part is symmetric
            if numel(var) > 1
                % Keep the ordering given by the user
                re = lava.symmetric(dim, var(1:dim*(dim+1)/2));
            else
                re = lava.symmetric(dim, vars(1));
            end
            
            % The imaginary part is antisymmetric
            im = lava.antiSymmetric(dim, vars(dim*(dim+1)/2+1:end));
            
            opOut = re + 1i*im;
        end

                 
        function out = num2lava(num)
            % converts a double matrix into a lava objects
            % opVar: [0] means constant 1
            if isnumeric(num) && ismatrix(num)
                out = lava(zeros(size(num)),num);
            else
                error('Input must be numeric matrix.')
            end
        end
        
        % a function converts numbers in some bases to decimal
        function [vDeci] = severalBases2dec(vIn,basesIn)
            % converts a vector vIn
            % denoting a number in a basis specified by basesIn
            % to  a single number in decimal
            % inspired by J.D. Bancal
            % by Cai Yu, 2019-11-20


            if length(basesIn)==1
                basesIn = ones(1,size(vIn,2))*basesIn;
            end

            if size(vIn,2)~=length(basesIn)
                error('length mismatch');
            end

            vDeci = zeros(size(vIn,1),1);
            for jj=1:size(vIn,2)
                % we read from left to right, but lowest power is on the right
                ii = size(vIn,2)-jj + 1;
                if vIn(:,ii)>=basesIn(ii)
                    error('number larger than basis')
                end
                vDeci = vDeci + vIn(:,ii).*basesIn(ii)^(jj-1);
            end

            end
        
         % static methods...
%          function M = buildLasserreMatrix(basisOp, coreOp)
%             % build the lava object representing the 
%             % localizing matrix of core localized by basis
%             fprintf('to be done')
%          end
        
     end
end


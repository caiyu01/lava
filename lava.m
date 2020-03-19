classdef lava
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
% x = lava({[1 2] 2; 3 4},[1 1;-1 1i]); % specifies the variables
% and coefficients
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
        % a = lava; % empty
        % a = lava([1 2;3 4]); % specifies the variables, but coeff=1
        % a = lava({[1 2] 2; 3 4},[1 1;-1 1i]); % specifies the variables
        % and coefficients
            switch nargin
                case 0
                    % if coeff1 is not specified
                    opOut = lava.empty; % size is zero
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
                        % match every cell with the maximum width and depth
                        for ii=1:numel(opVar)
                            if (size(opVar{ii},1)~=size(coeff{ii},1) || size(coeff{ii},2)~=1)... % if size does not match
                                    && ~isempty(opVar{ii}) && ~isempty(coeff{ii}) % and it is not empty
                                error('Wrong size. coeff must be an n-by-1 vector, where n is the height of opVar')
                            end
                        end
                            maxDepth = max(cellfun(@(x) size(x,1), opVar),[],'all');
                            maxWidth = max(cellfun(@(x) size(x,2), opVar),[],'all');
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
                        % sort, since commuting!
                        opVar = reshape(sort(reshape(opVar,numel(opVar)/maxWidth,maxWidth),2),m1,n1,maxDepth,maxWidth);
                        coeff = permute(reshape(cell2mat(coeff),maxDepth,m1,n1),[2,3,1]);
                        opOut.opVar = opVar;
                        opOut.coeff = coeff;
                    elseif isa(opVar,'double') && isa(coeff,'double')
                        % should allow construction from 4-D opVar with 3-D coeff
                        if ismatrix(opVar) && ismatrix(coeff)
                            % both are 2-D matrix, 
                            % interpreted depth=1 and width = 1
                            if size(opVar,1)~=size(coeff,1) || size(opVar,2)~=size(coeff,2)
                                error('Wrong size. opVar and coeff must have same size.')
                            end
                            opOut.opVar(:,:,1,1) = opVar;
                            opOut.coeff(:,:,1,1) = coeff;
                        elseif ~ismatrix(opVar) % && ~ismatrix(coeff)
                            % in case some dimensions has length 1
                            if size(opVar,1)==1 || size(opVar,2)==1
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
            if isempty(op)
                [varargout{1:nargout}] = size([],varargin{:});
            elseif nargin==1 && nargout<=1
                % output a vector, of only the matrix dimension
                % not depth and width of opVar
                opVar1 = op.opVar;
                varargout{1} = [size(opVar1,1),size(opVar1,2)];
            else
                [varargout{1:nargout}] = size(op.opVar,varargin{:});
            end
        end
        
        function n = numel(op)
            % only the first two matrix dimensions
            n = size(op,1)*size(op,2);
        end
        
        function l = length(op)
            [m,n] = size(op);
            if m>=n
                l = m;
            else
                l = n;
            end
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    subscripted referecing and assignment
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % Squeezed operators and coefficients access
        function result = sOpVar(op1)
            result = squeeze(op1.opVar);
        end
        
        function result = sCoeff(op1)
            result = squeeze(op1.coeff);
        end
        
        % subscripted referencing
        function varargout = subsref(op1,s)
            % Supported indexing
            % a.opVar, a.coeff
            % a.opVar(1), a.opVar(1,1) (not squeezed, use a(1,1).sOpVar for squeezed opVars!)
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
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    basic arithmetic operations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function opOut = plus(op1,op2)
        % addition
            [m1, n1, d1, w1] = size(op1);
            [m2, n2, d2, w2] = size(op2);
            % in case one of them is empty
            if isempty(op1)
                opOut = op2;
                return;
            elseif isempty(op2)
                opOut = op1;
                return;
            elseif max(m1,n1)==1
                % scalar with matrix
                op1 = simplify(kron(op1,ones(size(op2))));
                [m1, n1, d1, w1] = size(op1);
            elseif max(m2,n2)==1
                % matrix with scalar
                op2 = simplify(kron(op2,ones(size(op1))));
                [m2, n2, d2, w2] = size(op2);
            elseif m1~=m2 || n1~=n2
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
                if max(m1,n1)==1 
                    op1 = simplify(kron(op1,ones(size(op2))));
                    [m1, n1, d1, w1] = size(op1);
                elseif max(m2,n2)==1
                    op2 = simplify(kron(op2,ones(size(op1))));
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
            elseif isa(op1,'double') && isa(op2,'lava')
                opOut = lava(op2.opVar, op1.*op2.coeff);       
            elseif isa(op1,'lava') && isa(op2,'double')
                opOut = lava(op1.opVar, op2.*op1.coeff);
            end
        end
                
        % minus
        function opOut = minus(op1,op2)
            opOut = op1 + -1.*op2;
        end
        
        % mtimes
        function opOut = mtimes(op1,op2)  
            % matrix multiplication between
            % a double matrix with a lava matrix
            % a lava matrix with a double matrix
            % or two lava matrices
            if isa(op1,'double') && isa(op2,'lava')
                op1 = lava(zeros(size(op1)),op1);
                opOut = mtimes(op1,op2);
            elseif isa(op1,'lava') && isa(op2,'double')
                op2 = lava(zeros(size(op2)),op2);
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
            if ~isa(exponent, 'double') || (numel(exponent) ~= 1) || ~isequal(round(exponent), exponent) || (exponent < 0)
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
        function opOut = kron(op1,op2)
            % kronecker product between
            % a double matrix with a lava matrix
            % a lava matrix with a double matrix
            % or two lava matrices
            if isa(op1,'double') && isa(op2,'lava')
                op1 = lava(zeros(size(op1)),op1);
                opOut = kron(op1,op2);
            elseif isa(op1,'lava') && isa(op2,'double')
                op2 = lava(zeros(size(op2)),op2);
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
        
        % test if the object is a scalar
        function result = isscalar(op1)
            % Returns true iff op1 is 1x1 and a constant (i.e. involves no
            % variable).
            
            % First we check the dimension
            if numel(op1) > 1
                result = false;
                return;
            end
            
            % Now we have only one element. We check if it involves any
            % variables.
            op1 = op1.simplify;
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
            
            % Now we check the structure
            result = isequal(op1, op1.ctranspose);
        end
        
        % reshape
        function result = reshape(op1, varargin)
            if length(varargin)==1
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
                % need to standardize the width and depth
                [op1, op2] = matchSize(varargin{1},varargin{2});
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
                result = horzcat(varargin{1},horzcat(varargin{2:end}));
            elseif nargin==2
                % need to standardize the width and depth
                [op1, op2] = matchSize(varargin{1},varargin{2});
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
        
%         function opOut = setKron(op1,op2)
%             % set kronecker, for generating higher order relaxations
%             % only concerns the opVar, coeff are set to 1!
%             % disregards the size of inputs, only the unique list of opVars
%             % returns a row lava
%             % op = setKron(lava1,lava2)
%             op3 = kron(op1,op2);
%             opVar1 = uniquecell([op3.opVar]);
%             opOut = lava(opVar1(:))'; % all coefficients are real, so ' is ok
%         end
        
        function out = uniqueVar(varargin)
            % from several lava objects
            % returns the list of unique variables
            % either in a lava objects (default)
            % or a unique list (n-by-width) array
            % example: a = lava([1 2 3])
            %          out = uniqueVar(a,kron(a,a));
            %          out = uniqueVar(a,kron(a,a),'array');
            
            maxWidth = 1;
            list = [];
            if strcmp(varargin{end},'array')
                for ii=1:nargin-1
                    [~,~,~,w] = size(varargin{ii});
                    opVar1 = varargin{ii}.opVar;
                    tmp = reshape(opVar1,numel(opVar1)/w,w);
                    if w<maxWidth
                        tmp = [zeros(size(tmp,1),maxWidth-w) tmp];
                    elseif w>maxWidth
                        list = [zeros(size(list,1),w-maxWidth) list];
                        maxWidth = w;
                    end
                    list = unique([list; tmp],'rows');
                end
                out = list;
                return;
            else
                for ii=1:nargin
                    [~,~,~,w] = size(varargin{ii});
                    opVar1 = varargin{ii}.opVar;
                    tmp = reshape(opVar1,numel(opVar1)/w,w);
                    if w<maxWidth
                        tmp = [zeros(size(tmp,1),maxWidth-w) tmp];
                    elseif w>maxWidth
                        list = [zeros(size(list,1),w-maxWidth) list];
                        maxWidth = w;
                    end
                    list = unique([list; tmp],'rows');
                end
                out = lava(mat2cell(list,ones(1,size(list,1)),w));
            end
        end
        
        % toStr
        function result = toStr(op1)
            % This function returns a string description of the polynomial matrix
            result = cell(size(op1));
            
            for i = 1:size(op1,1)
                for j = 1:size(op1,2)
                    text = ' + ';
                    for k = 1:size(op1.opVar,3)
                        if (op1.coeff(i,j,k) ~= 0)
                            if (isreal(op1.coeff(i,j,k)) && (op1.coeff(i,j,k) < 0)) || ...
                               (isreal(1i*op1.coeff(i,j,k)) && (-1i*op1.coeff(i,j,k) < 0))
                                % correct the sign
                                text = [text(1:end-2), '- '];
                            end
                            if sum(op1.opVar(i,j,k,:)) == 0
                                % just a constant term
                                text = [text, num2str(abs(op1.coeff(i,j,k))), ' + '];
                            else
                                % some operators are involved
                                if isreal(op1.coeff(i,j,k))
                                    % purely real coefficient
                                    if abs(op1.coeff(i,j,k) ~= 1)
                                        % write the coefficient if not unity
                                        text = [text, num2str(abs(op1.coeff(i,j,k))), '*'];
                                    end
                                elseif isreal(1i*op1.coeff(i,j,k))
                                    % purely imaginary coefficient
                                    if abs(op1.coeff(i,j,k) ~= 1)
                                        % write the coefficient if not unity
                                        text = [text, num2str(abs(op1.coeff(i,j,k))), 'i*'];
                                    end
                                else
                                    % both real and imaginary
                                    % coefficients
                                    text = [text, '(', num2str(op1.coeff(i,j,k)), ')*'];
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
                        result{i,j} = text(4:end-3);
                    else
                        result{i,j} = '0';
                    end
                end
            end
        end
        
        % display
        function disp(op1,str1)
            if nargin==2 && strcmp(str1,'full')
                % possibly better display
                ...
            else
                % just the default display
                builtin('disp',op1);
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
            %
            % See also:
            %     lava.uniqueVar
            
            selectedVariables = uniqueVar(lava(0), list);
            base = selectedVariables*selectedVariables';
            opOut = kron(base, op1);
        end
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %    other operations
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        % simplify
        % slow: elementwise operation!
        % try arrayfun?
        function opOut = simplify(op1)
            [m1, n1, d1, w1] = size(op1);
            opVar1 = op1.opVar;
            coeff1 = op1.coeff;
            
            % sort the addition dimension, d
            % elementwise
            opVar1 = reshape(permute(reshape(opVar1,m1,n1,d1,w1),[3 4 2 1]),d1,numel(opVar1)/d1);
            coeff1 = reshape(permute(reshape(coeff1,m1,n1,d1), [3 2 1]),d1,numel(coeff1)/d1);
            
            for ii=1:m1*n1
                tmpV = opVar1(:,(1:w1)+(ii-1)*w1);
                tmpC = coeff1(:,ii);
                % sort according to opVar
                [tmpV, idx] = sortrows(tmpV); 
                tmpC = tmpC(idx); % the coeffs follows
                uniqV=  unique(tmpV,'rows');
                if size(uniqV,1)<size(tmpV,1)
                    for jj=1:size(uniqV,1)
                        % collapse the coeff into the last one
                        % and set opVar to all 0 after the coeff is removed
                        idx = prod(tmpV==uniqV(jj,:),2);
                        nbRemove = sum(idx)-1;
                        if nbRemove>=1
                            tmpC(find(idx,1,'last')) = sum(tmpC(logical(idx)));
                            tmpC(find(idx,nbRemove)) = 0;
                            tmpV(find(idx,nbRemove),:) = zeros(nbRemove,w1);
                        end
                    end
                end
                % move the 0 coefficient
                [tmpC, idx] = sort(tmpC);
                tmpV = tmpV(idx,:);
                opVar1(:,(1:w1)+(ii-1)*w1) = tmpV;
                coeff1(:,ii) = tmpC;
            end
            
            % trim the top zeros, if any
            % reduce maximum depth
%             idx = sum(abs(coeff1),2)==0;
            idx = prod(coeff1==0,2) == 1;
            opVar1(idx,:) = [];
            coeff1(idx,:) = [];
            % new maxDepth
            d1 = size(opVar1,1);
            
            % trim the left most zeros, if any
            % reduce maximal width
            opVar1 = sort(reshape(permute(reshape(opVar1,d1,w1,n1,m1),[4 3 1 2]),numel(opVar1)/w1,w1),2);
            coeff1 = reshape(permute(reshape(coeff1,d1,n1,m1),[3 2 1]),m1,n1,d1);
%             idx = sum(abs(opVar1),1)==0;
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
        
    end
    
     methods 
         
        % read out to variable index
        function [maxWidth, maxBase, uniqueIdx, varargout] = assignVarIdx(varargin)
            % [maxWidth, maxBase, uniqueIdx, varIdxCell] = assignVarIdx(lava1, lava2, {lava3, lava4},...)
            
            % in case some argin are cells of lavas
            co=1;
            for ii=1:nargin
                if iscell(varargin{ii})
                    for jj = 1:numel(varargin{ii})
                        if isa(varargin{ii}{jj},'lava')
                            input{co} = varargin{ii}{jj};
                            co = co+1;
                        else
                            error('Invalid input. Must be lava or cell of lavas')
                        end
                    end
                elseif isa(varargin{ii},'lava')
                    input{co} = varargin{ii};
                    co = co+1;
                else
                    error('Invalid input. Must be lava or cell of lavas')
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
            
            % varargout is a structure, with varIdx and coeff
            s1.type = '.';
            s1.subs = 'opVar';
            s2.type = '.';
            s2.subs = 'coeff';
            
            % put it in the output the same way as input
            % distinguish the case of lava and cell of lavas
            varIdxCell = cell(1,nargin);
            co = 1;
            for ii=1:nargin
                if iscell(varargin{ii})
                    for jj = 1:numel(varargin{ii})
                        [m, n, d, w] = size(input{co});
                        tmp = reshape(input{co}.opVar,[m*n*d, w]);
                        vvarIdx = reshape(lava.severalBases2dec(tmp,maxBase),m,n,d);
                        varIdxCell{ii}{jj}.varIdx = vvarIdx;
                        % just keeping a copy of the coefficients, for the
                        % corresponding varIdx
                        % maybe we don't need it
                        varIdxCell{ii}{jj}.coeff = input{co}.coeff; 
                        co = co+1;
                        % wrap vvarIdx into column vector to make the concatenation
                        uniqueIdx = unique([uniqueIdx; vvarIdx(:)]);
                    end
                else
                    [m, n, d, w] = size(input{co});
                    tmp = reshape(input{co}.opVar,[m*n*d, w]);
                    vvarIdx = reshape(lava.severalBases2dec(tmp,maxBase),m,n,d);
                    varIdxCell{ii}.varIdx = vvarIdx;
                    % just keeping a copy of the coefficients, for the
                    % corresponding varIdx
                    % maybe we don't need it
                    varIdxCell{ii}.coeff = input{co}.coeff; 
                    co = co+1;
                    % wrap vvarIdx into column vector to make the concatenation
                    uniqueIdx = unique([uniqueIdx; vvarIdx(:)]);
%                     vvarIdx = cellfun( @(x) (severalBases2dec(x,maxBase)), subsref(input{co},s1),'Uni',0);
%                     varIdxCell{ii}.varIdx = vvarIdx;
%                     % just keeping a copy of the coefficients, for the
%                     % corresponding varIdx
%                     % maybe we don't need it
%                     varIdxCell{ii}.coeff = subsref(input{co},s2);
%                     co = co+1;
%                     % wrap vvarIdx into column vector to make the concatenation
%                     % of uniqueIdx easier
%                     if iscell(vvarIdx)
%                         vvarIdx = cellfun(@(x) x', vvarIdx,'Uni',0);
%                         vvarIdx = [vvarIdx{:}];
%                     else
%                         vvarIdx = vvarIdx';
%                     end
%                     uniqueIdx = unique([uniqueIdx unique(vvarIdx(:))']);
                end
            end
            varargout = varIdxCell;
        end
        
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
            % lava(0) = 1;
            vSp(1) = 1;
            
            varargout = cell(1,nargin);
            for ii=1:nargin
                if iscell(varargin{ii})
                    % if it is a cell
                    for kk = 1:numel(varargin{ii})
                        % operate on each cell element
                        fprintf('Assigning cell %d matrix %d', ii, kk)
                        fprintf(newline)
                        clear m
                        tmpVarIdx = varIdxCell{ii}{kk}.varIdx;
                        tmpCoeff = varIdxCell{ii}{kk}.coeff;
                        
                        [m1,n1,d1,~] = size(varargin{ii}{kk});
                    
                        tmp = reshape(vSp(tmpVarIdx+1),m1,n1,d1);
                        varargout{ii}{kk} = sum(tmp.*tmpCoeff,3);
                        
                    end
                else
                    fprintf('Assigning matrix %d', ii)
                    fprintf(newline)
                    clear m
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
     end
     
     methods(Static)
                 
        function out = num2lava(num)
            % converts a double matrix into a lava objects
            % opVar: [0] means constant 1
            if isa(num,'double') && ismatrix(num)
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


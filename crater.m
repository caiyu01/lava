classdef crater
% Lava lakes are contained in volcanic craters.
% This object describes constraints on LAsserre VAriables

    properties
        cons % constraints
    end
    
    methods
        % Constructor for a single constaint
        function result = crater(lhs, relation, rhs)
            if ~isa(lhs, 'lava')
                lhs = lava.num2lava(lhs);
            end
            if ~isa(rhs, 'lava')
                rhs = lava.num2lava(rhs);
            end
            assert(isequal(relation, '<=') || isequal(relation, '==') || isequal(relation, '>='));
            
            result.cons = {struct('lhs', lhs, 'rel', relation, 'rhs', rhs)};
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Concatenation
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % concatenation of multiple constraints
        function result = horzcat(varargin)
            result = crater(0,'==',0); % we want to return a crater object
            result.cons = {};          % reinitialization
            for i = 1:nargin
                if ~isempty(varargin{i})
                    assert(isa(varargin{i}, 'crater'));
                    result.cons = [result.cons, varargin{i}.cons];
                end
            end
        end
        
        % vertical concatenation is the same
        function result = vertcat(varargin)
            result = horzcat(varargin);
        end
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Display
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        % disp
        function disp(cons1)
            if numel(cons1.cons) == 1
                disp('1 Polynomial constraint:')
            else
                disp([num2str(numel(cons1.cons)), ' Polynomial constraints:']);
            end
            
            for i = 1:numel(cons1.cons)
                dim1 = max(size(cons1.cons{i}.lhs,1), size(cons1.cons{i}.rhs,1));
                dim2 = max(size(cons1.cons{i}.lhs,2), size(cons1.cons{i}.rhs,2));
                text = ['  ', num2str(i), ' : ', num2str(dim1), 'x', num2str(dim2), ' '];
                switch cons1.cons{i}.rel
                    case '<='
                        text = [text, 'inequality <='];
                    case '=='
                        text = [text, 'equality'];
                    case '>='
                        text = [text, 'inequality >='];
                    otherwise
                        error('Unknown relation');
                end
                disp(text);
            end
        end

        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        % Yalmip interface
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        function [sol, Fv, objv, lhsv, rhsv] = optimize(varargin)
            t = tic;
            assert(nargin >= 2);
            if isa(varargin{2}, 'crater')
                error('Constraints should be put as the first parameter of solvesdp.');
            end
            assert(isa(varargin{1}, 'crater'));
            assert(isa(varargin{2}, 'lava'));
            F = varargin{1};
            obj = varargin{2};
            options = [];
            if nargin >= 3
                options = varargin{3};
            end
            
            % First, we translate all lava objects into SDP variables
            nbConstr = length(F.cons);
            lhs = cell(1,nbConstr);
            rhs = cell(1,nbConstr);
            for i = 1:nbConstr
                lhs{i} = F.cons{i}.lhs;
                rhs{i} = F.cons{i}.rhs;
            end
            
            [objv, lhsv, rhsv] = assignSdpVar(obj, lhs, rhs);

            Fv = [];
            for i = 1:nbConstr
                switch F.cons{i}.rel
                    case '<='
                        Fv = [Fv, lhsv{i} <= rhsv{i}];
                    case '=='
                        Fv = [Fv, lhsv{i} == rhsv{i}];
                    case '>='
                        Fv = [Fv, lhsv{i} >= rhsv{i}];
                    otherwise
                        error('Unknown relation');
                end
            end
            tlava = toc(t);
            
            % Solve the problem through yalmip
            sol = optimize(Fv, objv, options);
            
            sol.lavatime = tlava;
        end
        
    end

end

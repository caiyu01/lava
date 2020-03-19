classdef lavacons
% Object describing a set of Constraints on LAserre VAriables

    properties
        cons % constraints
    end
    
    methods
        % Constructor for a single constaint
        function result = lavacons(lhs, relation, rhs)
            if ~isa(lhs, 'lava')
                lhs = lava.num2lava(lhs);
            end
            if ~isa(rhs, 'lava')
                rhs = lava.num2lava(rhs);
            end
            assert(isequal(relation, '<=') || isequal(relation, '==') || isequal(relation, '>='));
            
            result.cons = {struct('lhs', lhs, 'rel', relation, 'rhs', rhs)};
        end
        
        % concatenation of multiple constraints
        function result = horzcat(varargin)
            result = varargin{1}; % we want to return a lavacons object
            result.cons = {};     % reinitialization
            for i = 1:nargin
                assert(isa(varargin{i}, 'lavacons'));
                result.cons = [result.cons, varargin{i}.cons];
            end
        end
        
        % vertical concatenation is the same
        function result = vertcat(varargin)
            result = horzcat(varargin);
        end

        % disp
        function disp(cons1)
            if numel(cons1.cons) == 1
                disp('Polynomial constraint:')
            else
                disp('Polynomial constraints:')
            end
            
            for i = 1:numel(cons1.cons)
                dim1 = max(size(cons1.cons{i}.lhs,1), size(cons1.cons{i}.rhs,1));
                dim2 = max(size(cons1.cons{i}.lhs,2), size(cons1.cons{i}.rhs,2));
                text = [num2str(dim1), 'x', num2str(dim2), ' '];
                switch cons1.cons{i}.rel
                    case '<='
                        text = [text, 'inequality <='];
                    case '=='
                        text = [text, 'equality'];
                    case '>='
                        text = [text, 'inequality >='];
                end
                disp(text);
            end
        end
           
    end

end

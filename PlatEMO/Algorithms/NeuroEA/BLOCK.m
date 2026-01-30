classdef BLOCK < handle & matlab.mixin.Heterogeneous
%BLOCK - The superclass of blocks.

%------------------------------- Copyright --------------------------------
% Copyright (c) 2026 BIMK Group. You are free to use the PlatEMO for
% research purposes. All publications which use this platform or any code
% in the platform should acknowledge the use of "PlatEMO" and reference "Ye
% Tian, Ran Cheng, Xingyi Zhang, and Yaochu Jin, PlatEMO: A MATLAB platform
% for evolutionary multi-objective optimization [educational forum], IEEE
% Computational Intelligence Magazine, 2017, 12(4): 73-87".
%--------------------------------------------------------------------------

    properties(SetAccess = protected)
        parameter;      % Parameters in the block
        lower;          % Lower bound of each parameter
        upper;          % Upper bound of each parameter
        output;         % Current output of the block
        nextOut = 1;	% Index of next output solution
    end
    properties
        trainTime = 0;	% Number of training times
    end
    methods
        function Main(obj,Problem,Precursors,Ratio)
        %Main - Main procedure of the block.
        %
        %   obj.Main(Pro,Pre,Ratio) performs the main procedure of block
        %   obj. Pro is a PROBLEM object, Pre are multiple BLOCK objects,
        %   and Ratio are the ratios of solutions gathered from the outputs
        %   of blocks Pre. The output of the main procedure is stored in
        %   obj.output.
        %
        %   Example:
        %       Blocks(i).Main(Problem,Blocks(logical(Graph(:,i))),Graph(:,i));
        end
        function ParameterAssign(obj)
        %ParameterAssign - Assign parameters to block-specific variables.
        %
        %   This function is automatically called when the value of
        %   obj.parameter is changed.
        end
    end
	methods(Sealed = true)     
        function ParameterSet(obj,value)
        %ParameterSet - Set the parameters of multiple blocks.
        %
        %   obj.ParameterSet(Par) sets obj.parameter to Par. Here obj can
        %   be multiple BLOCK objects, where Par is a vector concatenating
        %   the parameters of all the objects.
        %
        %   obj.ParameterAssign() is automatically called after this
        %   function.
        %
        %   Example:
        %       Blocks(1:5).ParameterSet(Par);
        
            k = 1;
            for i = 1 : length(obj)
                obj(i).parameter = min(obj(i).upper,max(obj(i).lower,value(k:k-1+length(obj(i).parameter))));
                obj(i).ParameterAssign();
                obj(i).trainTime = obj(i).trainTime + 1;
                k = k + length(obj(i).parameter);
            end
        end
        function value = parameters(obj)
        %parameters - Get the parameters of multiple blocks.
        %
        %   Par = obj.parameters returns a vector concatenating the
        %   parameters of multiple blocks obj.
        
            value = cat(2,obj.parameter);
        end
        function value = lowers(obj)
        %lowers - Get the lower bounds of the parameters of multiple blocks.
        %
        %   Lower = obj.lowers returns a vector concatenating the lower
        %   bounds of the parameters of multiple blocks obj.
        
            value = cat(2,obj.lower);
        end
        function value = uppers(obj)
        %uppers - Get the upper bounds of the parameters of multiple blocks.
        %
        %   Upper = obj.uppers returns a vector concatenating the upper
        %   bounds of the parameters of multiple blocks obj.
        
            value = cat(2,obj.upper);
        end
        function Output = Gather(obj,Problem,Predecessors,Ratio,inType,multiple)
        %Gather - Gather the output from multiple precursors
        %
        %   Output = obj.Gather(Pro,Pre,Ratio,inType,mul) gathers outputs
        %   from multiple blocks for block obj. Pro is a PROBLEM object,
        %   Pre are multiple BLOCK objects, Ratio are the ratios of
        %   solutions gathered from the outputs of blocks Pre, inType is
        %   the type of gathered solutions (1. SOLUTION objects, 2. decision
        %   matrix), and mul indicates that the number of gathered
        %   solutions should be a multiple of mul.
        %
        %   This function is usually called at the beginning of obj.Main.
        %
        %   Example:
        %       Population = obj.Gather(Problem,Predecessors,Ratio,1,1);
        %       ParentDec = obj.Gather(Problem,Predecessors,Ratio,2,2);
        
            Ratio  = Ratio(Ratio>0);
            Output = [];
            Lens   = [];
            for i = 1 : length(Predecessors)
                if inType == 1
                    % Get SOLUTION objects from the outputs of predecessors
                    if ~isa(Predecessors(i).output,'SOLUTION')
                        Predecessors(i).output = Problem.Evaluation(Predecessors(i).output);
                    end
                    Out    = Predecessors(i).output;
                    Index  = mod(Predecessors(i).nextOut-1:Predecessors(i).nextOut+floor(length(Out)*Ratio(i))-2,length(Out)) + 1;
                    Output = [Output,Out(Index)];
                    Lens   = [Lens,length(Index)];
                    Predecessors(i).nextOut = mod(Index(end),length(Out)) + 1;
                else
                    % Get decision matrix from the outputs of predecessors
                    if ~isa(Predecessors(i).output,'SOLUTION')
                        Out = Predecessors(i).output;
                    else
                        Out = Predecessors(i).output.decs;
                    end
                    Index  = mod(Predecessors(i).nextOut-1:Predecessors(i).nextOut+floor(size(Out,1)*Ratio(i))-2,size(Out,1)) + 1;
                    Output = [Output;Out(Index,:)];
                    Lens   = [Lens,length(Index)];
                    Predecessors(i).nextOut = mod(Index(end),size(Out,1)) + 1;
                end
            end
            % Interleave the outputs from multiple precursors
            Index = repmat(1:max(Lens),length(Lens),1) + repmat([0,cumsum(Lens(1:end-1))]',1,max(Lens));
            Index = min(Index,sum(Lens));
            Index = unique(Index(:),'stable');
            if inType == 1
                Output = Output(Index(1:floor(end/multiple)*multiple));
            else
                Output = Output(Index(1:floor(end/multiple)*multiple),:);
            end
            % Reset the next output solution of the current block
            obj.nextOut = 1;
        end
        function Validity(obj,Graph)
        %Validity - Check the validity of an algorithm with mulitple blocks
        %
        %   obj.Validity(G) throws an error if the algorithm is invalid.
        %   obj is an array of BLOCK objects constituting the algorithm and
        %   G is the adjacency matrix. After the error err is caught, use
        %   err.identifier to determine the error type and use
        %   str2num(err.cause{1}.message) to obtain the indexes of invalid
        %   objects. Besides, nothing happens if the algorithm is valid.
        %
        %   Example:
        %       try
        %           Blocks.Validity(Graph);
        %       catch err
        %           switch err.identifier
        %               case 'BLOCK:NoInput'
        %                   str2num(err.cause{1}.message)
        %               case 'BLOCK:NoOutput'
        %                   str2num(err.cause{1}.message)
        %               case ...
        %                   ...
        %           end
        %       end

            try
                type = strrep(arrayfun(@class,obj,'UniformOutput',false),'Block_','');
                G    = digraph(Graph);
                invalidStr = '';
                assert(strcmp('Population',type{1}),'BLOCK:NoPopulation','the first block is not a population.',invalidStr);
                assert(any(ismember({'Crossover','Exchange','Mutation'},type)),'BLOCK:NoOperator','the algorithm does not contain variation operator.',invalidStr);
                invalidStr = num2str(find(indegree(G)==0)');
                assert(isempty(invalidStr),'BLOCK:NoInput','the block #%s have no predecessor.',invalidStr);
                invalidStr = num2str(find(outdegree(G)==0)');
                assert(isempty(invalidStr),'BLOCK:NoOutput','the block #%s have no successor.',invalidStr);
                invalidStr = num2str(find(conncomp(G,'Type','weak')>1));
                assert(isempty(invalidStr),'BLOCK:Isolation','the block #%s are isolated.',invalidStr);
                invalidStr = num2str(find(diag(Graph))');
                assert(isempty(invalidStr),'BLOCK:SelfLoop','the block #%s have self-loop.',invalidStr);
                if ismethod(G,'allcycles')
                    cycles     = allcycles(G);
                    invalidStr = num2str(cell2mat(cycles(find(cellfun(@(a)~ismember('Population',type(a)),cycles),1))));
                    assert(isempty(invalidStr),'BLOCK:InfLoop','the cycle #%s have no population.',invalidStr);
                end
            catch err
                err = addCause(err,MException('',invalidStr));
                rethrow(err);
            end
        end
    end
end
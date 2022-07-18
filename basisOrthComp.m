% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %
%                                                                           %
%    basisOrthComp                                                          %
%                                                                           %
%                                                                           %
% OUTPUT: Returns details about the network and its stoichiometric          %
%    subspace, its orthogonal complement of the stoichiometric subspace,    %
%    and the basis for the orthogonal complement of the stoichiometric      %
%    subspace with labels for the nonpivot species. The output variables    %
%    'model', 'basis', and 'nonpivot' allow the user to view the following, %
%    respectively:                                                          %
%       - Complete network with all the species listed in the 'species'     %
%            field of the structure 'model'                                 %
%       - Table showing the basis for the orthogonal complement of the      %
%            stoichiometric subspace                                        %
%       - Species representing nonpivot columns                             %
%                                                                           %
% INPUT: model: a structure, representing the chemical reaction network     %
%    (CRN) (see README.txt for details on how to fill out the structure)    %
%                                                                           %
% Note: The idea for the first part of the algorithm comes from Soranzo and %
%    Altafini (2009).                                                       %
%                                                                           %
% Reference: Soranzo N, Altafini C (2009) ERNEST: a toolbox for chemical    %
%    reaction network theory. Bioinform 25(21):2853–2854.                   %
%    https://doi.org/10.1093/bioinformatics/btp513                          %
%                                                                           %
% Created: 15 July 2022                                                     %
% Last Modified: 18 July 2022                                               %
%                                                                           %
% % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % % %



function [model, basis, nonpivot] = basisOrthComp(model)

    %
    % Step 1: Add to 'model.species' all species indicated in the reactions
    %
    
    % Initialize model.species
    model.species = { };
    
    % Get all species from reactants
    for i = 1:numel(model.reaction)
        for j = 1:numel(model.reaction(i).reactant)
            model.species{end+1} = model.reaction(i).reactant(j).species;
        end
    end
    
    % Get species from products
    for i = 1:numel(model.reaction)
        for j = 1:numel(model.reaction(i).product)
            model.species{end+1} = model.reaction(i).product(j).species;
        end
    end
    
    % Get only unique species
    model.species = unique(model.species);
    
    % Count the number of species
    m = numel(model.species);
    
    
    
    %
    % Step 2: Get the matrix of reaction vectors of the network
    %
    
    % Initialize the matrix of reactant complexes
    reactant_complex = [ ];
    
    % Initialize the matrix of product complexes
    product_complex = [ ];
    
    % Initialize the stoichiometric matrix
    N = [ ];
    
    % For each reaction in the model
    for i = 1:numel(model.reaction)
      
        % Initialize the vector for the reaction's reactant complex
        reactant_complex(:, end+1) = zeros(m, 1);
        
        % Fill it out with the stoichiometric coefficients of the species in the reactant complex
        for j = 1:numel(model.reaction(i).reactant)
            reactant_complex(find(strcmp(model.reaction(i).reactant(j).species, model.species), 1), end) = model.reaction(i).reactant(j).stoichiometry;
        end
        
        % Initialize the vector for the reaction's product complex
        product_complex(:, end+1) = zeros(m, 1);
        
        % Fill it out with the stoichiometric coefficients of the species in the product complex
        for j = 1:numel(model.reaction(i).product)
            product_complex(find(strcmp(model.reaction(i).product(j).species, model.species), 1), end) = model.reaction(i).product(j).stoichiometry;
        end
        
        % Create a vector for the stoichiometric matrix: Difference between the two previous vectors
        N(:, end+1) = product_complex(:, end) - reactant_complex(:, end);
        
        % If the reaction is reversible
        if model.reaction(i).reversible
          
            % Insert a new vector for the reactant complex: make it same as the product complex
            reactant_complex(:, end+1) = product_complex(:, end);
            
            % Insert a new vector for the product complex: make it the same as the reactant complex
            product_complex(:, end+1) = reactant_complex(:, end-1);
            
            % Insert a new vector in the stoichiometric matrix: make it the additive inverse of the vector formed earlier
            N(:, end+1) = -N(:, end);
        end
    end
    
    % Matrix of reaction vectors
    R = N';
    
    
    
    %
    % Step 3: Get the basis for the orthogonal complement of the stoichiometric subspace
    %
    
    % Form basis_R for the rowspace of R (so we'll work with R = span(basis_R))
    % Write R in reduced row echelon form: the transpose of R is used so 'pivot_row' will give the pivot rows of R
    %    - 'A' is W in reduced row echelon form
    %    - 'pivot_row' gives the row numbers of R which form a basis for the rowspace of R
    [~, pivot_row] = rref(R');
    
    % Form the basis
    basis_R = R(pivot_row, :);
    
    % Write the basis for R in rref: This is equivalent to solving x in basis_R*x = 0
    [B, basis_R_rref_pivot_row] = rref(basis_R);
    
    % Get the nonpivots
    nonpivot = setdiff(1:size(R, 2), basis_R_rref_pivot_row);
    
    % Prepare the basis vectors of the orthogonal complement of R: This is where we'll put the solution to basis_R*x = 0
    basis = zeros(numel(nonpivot), size(R, 2));
    
    % Go through each row of B
    for i = 1:size(B, 1)
        
        % Get the pivot of the row
        pivot = find(B(i, :), 1, 'first');
        
        % Check only the nonpivots
        for j = 1:numel(nonpivot)
            
            % If the nonpivot column in the row is nonzero
            if B(i, nonpivot(j)) ~= 0
                
                % The nonpivot is part of the linear combo of the pivot
                basis(j, pivot) = -B(i, nonpivot(j));
            end
        end
    end
    
    % The nonpivots are linear combos of themselves
    for i = 1:numel(nonpivot)
        basis(i, nonpivot(i)) = 1;
    end
    
    % Basis for the orthogonal complement
    basis = basis';
    
    
    
    %
    % Step 4: Display result
    %
    
    disp(['We have a vector subspace of R', num2str(size(R, 2))]);
    disp(['Rank of the stoichiometric subspace: ', num2str(rank(R))]);
    disp(['Rank of its orthogonal complement: ', num2str(rank(basis))]);
    
    % Get list of species and make it a table
    species = table(model.species');
    
    % Convert matrix of basis vectors to table
    basis = array2table(basis);
    
    % Combine the two tables
    T = [species, basis];
    
    % The header of the basis vectors should be the nonpivot species
    names = ["Species", string(model.species(nonpivot))];
    
    % Rename the header of the table
    T = renamevars(T, 1:width(T), names);
    
    % Display the table of basis vectors
    fprintf('\n')
    disp(T)

end
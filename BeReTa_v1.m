function BeReTaSolution = BeReTa_v1(model, reg_net, exp_comp, product_rxn)
% Beneficial Regulator Targeting (BeReTa) algorithm [Minsuk Kim et al, Bioinformatics, Accepted]
%
%INPUTS
% model             COBRA model structure
% reg_net           Transcriptional regulatory network structure
%       regulator              Cell containing regulator names
%       target                 Cell containing target names
%       sign                   Vector containing sign of regulation (activator, 1; repressor, -1; unknown, 0) (optional)
%       regulator, target, and sign should have same length
% exp_comp          Gene expression compendium data structure
%       genes                   Cell containing gene names
%       expression              Gene expression levels
%       genes and expression should have same length
% product_rxn       Exchange reaction for target product
%
%OUTPUTS
% BeReTaSolution    BeReTa solution structure (target regulator, beneficial score, p-value)

% Generate regulatory strength matrix (RegStr) and flux slope vector (q_slope)
[RegStr, RegStr_gene] = generateRegStrMatrix(reg_net, exp_comp, model);
q_slope = calculateFluxSlope(model, product_rxn);

% Trimming
reg_list = RegStr(2:end,1);
RegStr = cell2mat(RegStr(2:end,2:end));
q_slope = q_slope(1:size(RegStr,2));

% Calculate beneficial scores for all transcriptional regulators (TRs)
beneficial_score = RegStr * q_slope;

% Perform permutation test
beneficial_score_pval = permutationTestBeReTa(model, RegStr, q_slope, beneficial_score, reg_list);

% Apply target criteria to select BeReTa targets
BeReTaSolution = selectBeReTaTargets(model, RegStr, RegStr_gene, q_slope, beneficial_score, beneficial_score_pval, reg_list);

end


function [RegStr, RegStr_gene] = generateRegStrMatrix(reg_net, exp_comp, model)
%% Generate regulatory strength matrix (RS)
regulator = reg_net.regulator;
target = reg_net.target;
if length(fields(reg_net)) == 3
    sign = reg_net.sign;
else
    sign = zeros(size(regulator));
end
expression = exp_comp.expression;
expressionid = exp_comp.genes;
model_genes = model.genes;
model_rxns = model.rxns;

% Discard regulatory interactions for which correlation cannot be computed.
regulator_temp = {};
target_temp = {};
sign_temp = [];
for i = 1:length(regulator)
    regulator_id = find(strcmp(regulator(i), expressionid));
    target_id = find(strcmp(target(i), expressionid));
    if (length(regulator_id) == 1) && (length(target_id) == 1)
        regulator_temp = [regulator_temp; regulator(i)];
        target_temp = [target_temp; target(i)];
        sign_temp = [sign_temp; sign(i)];
    end
end
regulator = regulator_temp;
target = target_temp;
sign = sign_temp;

% Generate regulatory strength (RS) matrix for TR-gene interactions
reg_list = unique(regulator);
corr_mat1 = zeros(length(reg_list), length(model_genes));
corr_mat1_row = reg_list;
corr_mat1_col = model_genes';
for i = 1:length(regulator)
    row_id = find(strcmp(regulator(i), corr_mat1_row));
    col_id = find(strcmp(target(i), corr_mat1_col));
    x1 = expression(find(strcmp(regulator(i), expressionid)),:);
    x2 = expression(find(strcmp(target(i), expressionid)),:);
    corr_val = corr2(x1, x2);
    if sign(i) == 0
        corr_mat1(row_id, col_id) = corr_val;
    elseif sign(i) == 1
        corr_mat1(row_id, col_id) = abs(corr_val);
    elseif sign(i) == -1
        corr_mat1(row_id, col_id) = - abs(corr_val);
    end
end
% Remove zero columns in the matrix
row_index = find(sum(abs(corr_mat1')));
col_index = find(sum(abs(corr_mat1)));
corr_mat2 = corr_mat1(row_index,col_index);
corr_mat2_row = corr_mat1_row(row_index);
corr_mat2_col = corr_mat1_col(col_index);
RegStr_gene = [['RSmat_gene', corr_mat2_col]; [corr_mat2_row, num2cell(corr_mat2)]];

% Generate regulatory strength (RS) matrix for TR-rxn interactions
% GPR mapping using average function
corr_mat3 = zeros(length(corr_mat2_row), length(model_rxns));
corr_mat3_row = corr_mat2_row;
corr_mat3_col = model_rxns';
for i = 1:size(corr_mat2,1)
    gpr_genes = corr_mat2_col';
    gpr_levels = corr_mat2(i,:)';
    gpr_genes2 = gpr_genes(find(gpr_levels));
    gpr_levels2 = gpr_levels(find(gpr_levels));
    levels = gene_to_reaction_levels(model, gpr_genes2, gpr_levels2, @(x,y)((x+y)/2), @(x,y)((x+y)/2));
    corr_mat3(i,find(~isnan(levels))) = levels(find(~isnan(levels)));
end
RegStr = [['RSmat', corr_mat3_col]; [corr_mat3_row, num2cell(corr_mat3)]];

end
function reaction_levels = gene_to_reaction_levels( model, genes, levels, f_and, f_or )
% Original authors: Daniel Machado and Markus Herrgard, PLoS Computational Biology, 2014
% Code obtained from https://github.com/cdanielmachado/transcript2flux
%
% Convert gene expression levels to reaction levels using GPR associations.
% Level is NaN if there is no GPR for the reaction or no measured genes.
%
% INPUTS
%       model - cobra model
%       genes - gene names
%       levels - gene expression levels
%       f_and - function to replace AND
%       f_or - function to replace OR
%
% OUTPUTS
%       reaction_levels - reaction expression levels
%
% Author: Daniel Machado, 2013

    reaction_levels = zeros(length(model.rxns), 1);

    for i = 1:length(model.rxns)
        level = eval_gpr(model.grRules{i}, genes, levels, f_and, f_or);
        reaction_levels(i) = level;
    end

end
function [result, status] = eval_gpr(rule, genes, levels, f_and, f_or)
% Original authors: Daniel Machado and Markus Herrgard, PLoS Computational Biology, 2014
% Code obtained from https://github.com/cdanielmachado/transcript2flux
%
% Evaluate the expression level for a single reaction using the GPRs.
% Note: Computes the expression level even if there are missing measured
% values for the given rule. This implementation is a modified version of
% an implementation provided in [Lee et al, BMC Sys Biol, 2012]

    EVAL_OK = 1;
    PARTIAL_MEASUREMENTS = 0;
    NO_GPR_ERROR = -1;
    NO_MEASUREMENTS = -2;
    MAX_EVALS_EXCEEDED = -3;

    MAX_EVALS = 1000;
    NONETYPE = 'NaN';

    NUMBER = '[0-9\.\-e]+';
    MAYBE_NUMBER = [NUMBER '|' NONETYPE];

    expression = rule;
    result = NaN;
    status = EVAL_OK;

    if isempty(expression)
        status = NO_GPR_ERROR;
    else
        rule_genes = setdiff(regexp(expression,'\<(\w|\-)+\>','match'), {'and', 'or'});
        
        total_measured = 0;
        
        for i = 1:length(rule_genes)
            j = find(strcmp(rule_genes{i}, genes));
            if isempty(j)
                level = NONETYPE;
            else
                level = num2str(levels(j));
                total_measured = total_measured + 1;
            end
            expression = regexprep(expression, ['\<', rule_genes{i}, '\>'], level );
        end
        
        
        if total_measured == 0
            status = NO_MEASUREMENTS;
        else
            if total_measured < length(rule_genes)
                status = PARTIAL_MEASUREMENTS;
            end
            
            maybe_and = @(a,b)maybe_functor(f_and, a, b);
            maybe_or = @(a,b)maybe_functor(f_or, a, b); 
            str_wrapper = @(f, a, b)num2str(f(str2double(a), str2double(b)));

            counter = 0;
            
            while isnan(result)

                counter = counter + 1;
                if counter > MAX_EVALS
                    status = MAX_EVALS_EXCEEDED;
                    break
                end

                try 
                    result = eval(expression);            
                catch e   
                    paren_expr = ['\(\s*(', MAYBE_NUMBER,')\s*\)'];
                    and_expr = ['(',MAYBE_NUMBER,')\s+and\s+(',MAYBE_NUMBER,')'];
                    or_expr = ['(',MAYBE_NUMBER,')\s+or\s+(',MAYBE_NUMBER,')'];

                    expression = regexprep(expression, paren_expr, '$1');
                    expression = regexprep(expression, and_expr, '${str_wrapper(maybe_and, $1, $2)}');
                    expression = regexprep(expression, or_expr, '${str_wrapper(maybe_or, $1, $2)}');
                end
            end
            
        end
    end

end
function c = maybe_functor(f, a, b)
% Original authors: Daniel Machado and Markus Herrgard, PLoS Computational Biology, 2014
% Code obtained from https://github.com/cdanielmachado/transcript2flux
    
    if isnan(a) && isnan(b)
        c = nan;
    elseif ~isnan(a) && isnan(b)
        c = a;
    elseif isnan(a) && ~isnan(b)
        c = b;
    else 
        c = f(a,b);
    end
end

function q_slope = calculateFluxSlope(model, product_rxn)
%% Calculate flux slope vector (q_slope)
% Calculate initial and maximum flux of product Rxn.
initial_sol = optimizeCbModel(model);
product_rxn_id = ismember(model.rxns, product_rxn);
initial_product_flux = initial_sol.x(product_rxn_id);
model_product = changeObjective(model, product_rxn);
maximum_sol = optimizeCbModel(model_product);
maximum_product_flux = maximum_sol.x(product_rxn_id);

% Flux scanning while increasing the product flux
num_steps = 20;
flux_values = zeros(length(model.rxns), num_steps);
for i = 1:num_steps
    product_flux = initial_product_flux + (maximum_product_flux - initial_product_flux) * ((i-1)/num_steps);
    model_enforced = changeRxnBounds(model, product_rxn, product_flux, 'b');
    enforced_sol = optimizeCbModel(model_enforced, 'max', 'one');
    flux_values(:,i) = abs(enforced_sol.x);
end

% Use linear regression to estimate flux slopes
flux_values_product = flux_values(product_rxn_id,:);
q_slope = zeros(size(model.rxns));
for i = 1:length(model.rxns)
    linear_regression = polyfit(flux_values_product,flux_values(i,:),1);
    q_slope(i,1) = linear_regression(1);
end
q_slope(q_slope<0) = 0;

end

function beneficial_score_pval = permutationTestBeReTa(model, RegStr, q_slope, beneficial_score, reg_list)
%% Calculate p-values for beneficial scores
% Permute only for gene-associated rxns
selExc = findExcRxns(model, 1);
exchanges = model.rxns(selExc);
orphans = findOrphanRxns(model);
gene_associated = setdiff(model.rxns, exchanges);
gene_associated = setdiff(gene_associated, orphans);
rand_rxn_ids = find(ismember(model.rxns, gene_associated));

% Calculate permuted beneficial scores
rand_num = 10000;
beneficial_score_rand = zeros(length(reg_list), rand_num);
for i = 1:rand_num
    [~, rand_order] = sort(rand(size(rand_rxn_ids)));
    q_slope_rand = q_slope;
    q_slope_rand(rand_rxn_ids(rand_order)) = q_slope(rand_rxn_ids);
    beneficial_score_rand(:,i) = RegStr * q_slope_rand;
end

% Calculate p-values
beneficial_score_pval = ones(size(reg_list));
for i = 1:length(reg_list)
    if beneficial_score(i) > 0 
        beneficial_score_pval(i) = length(find(beneficial_score_rand(i,:) > beneficial_score(i))) / rand_num;
    elseif beneficial_score(i) < 0
        beneficial_score_pval(i) = length(find(beneficial_score_rand(i,:) < beneficial_score(i))) / rand_num;
    end
end

end

function BeReTaSolution = selectBeReTaTargets(model, RegStr, RegStr_gene, q_slope, beneficial_score, beneficial_score_pval, reg_list)
%% Apply target criteria to select BeReTa targets
% Define target criteria   (1) TR should have non-zero beneficial score.
pval_cut = 0.05;         % (2) The p-value of the beneficial score should be less than 0.05.
n_target_cut = 2;        % (3) TR should have two or more effective gene/reaction targets.
f_target_cut = 0.1;      % (4) At least 10% of target metabolic genes of TR should be beneficial, i.e. have positive flux slopes.

% Calculate the metrics
n_effective_genes = zeros(size(reg_list));
n_effective_rxns = zeros(size(reg_list));
f_effective_genes = zeros(size(reg_list));
for i = 1:length(reg_list)
    regulated_genes = RegStr_gene(1,(find(cell2mat(RegStr_gene(i+1,2:end)))+1))';
    effective_rxns = model.rxns(find(q_slope .* RegStr(i,:)'));
    n_effective_rxns(i) = length(effective_rxns);
    if length(effective_rxns) > 0
        effective_genes_temp = findGenesFromRxns(model, effective_rxns);
        effective_genes = {};
        for j = 1:length(effective_genes_temp)
            effective_genes = [effective_genes; effective_genes_temp{j}];
        end
        effective_genes = unique(effective_genes);
        effective_genes = intersect(effective_genes, regulated_genes);
        n_effective_genes(i) = length(effective_genes);
        f_effective_genes(i) = length(effective_genes)/length(regulated_genes);
    else
        n_effective_genes(i) = 0;
        f_effective_genes(i) = 0;
    end
end

% Select BeReTa targets
BeReTa_metrics = [beneficial_score, beneficial_score_pval, n_effective_genes, n_effective_rxns, f_effective_genes];
[~, score_order] = sort(abs(beneficial_score), 'descend');
reg_list = reg_list(score_order);
BeReTa_metrics = BeReTa_metrics(score_order,:);
BeReTaSolution_reg_list = {};
BeReTaSolution_metrics = [];
for i = 1:length(reg_list)
    if (BeReTa_metrics(i,2) < pval_cut) && (BeReTa_metrics(i,3) >= n_target_cut) && (BeReTa_metrics(i,4) >= n_target_cut) && (BeReTa_metrics(i,5) >= f_target_cut)
        BeReTaSolution_reg_list = [BeReTaSolution_reg_list; reg_list(i)];
        BeReTaSolution_metrics = [BeReTaSolution_metrics; BeReTa_metrics(i,:)];
    end
end
BeReTaSolution = [BeReTaSolution_reg_list, num2cell(BeReTaSolution_metrics(:,1:2))];

end

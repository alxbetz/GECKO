
initCobraToolbox
changeCobraSolver('ibm_cplex','all')


creMixo = readCbModel('/Users/alx/data/creinhardtii/metabolism/gmm/iCre1355_mixo.xml');
creMixo = changeObjective(creMixo,'Biomass_Chlamy_mixo')
%add pseudo reactions for DNA, RNA, Lipid(BB and Chain) Protein and Carbon


[ecModel,ecModel_batch] = enhanceGEM(creMixo)

modelRaven = ravenCobraWrapper(creMixo);
org_name = 'chlamydomonas reinhardtii';
name = 'creMixoAdapted'
version = '1.0'
%preprocess model
cd change_model
[model,name,version] = preprocessModel(modelRaven,name,version);

%Retrieve kcats & MWs for each rxn in model:
cd ../get_enzyme_data
model_data = getEnzymeCodes(model);
kcats      = matchKcats(model_data,org_name);

%Integrate enzymes in the model:
cd ../change_model
ecModel                 = readKcatData(model_data,kcats);
%[ecModel,modifications] = manualModifications(ecModel);
%GAM and NGAM values from Imam / Schaeuble 2015
%GAM – 92.43 mmol ATP gDW1, NGAM – 2.85 mmol ATP gDW1 h1)

sigma    = 0.5;      %Optimized for glucose
Ptot     = 0.303;      %Assumed constant ; totale measured protein content in g/gDW mixo : 0.303, auto: 0.261, hetero: 0.222 taken from Boyle & Morgan 2009
gR_exp   = 0.066 ;     %[g/gDw h] mixotrophic from Boyle & Morgan 2009
GAM = 92.43;
%c_source = 'D-glucose exchange (reversible)'; %Rxn name for the glucose uptake reaction



%5.4967e-02
%EX_co2_e

c_%source = 'EX_co2_e';
c_source = 'CO2 exchange';

cd ../limit_proteins

%infIDX = find(~isinf(table2array(proteinAbundanceTest(:,2))));
%proteinAbundanceTest1 = proteinAbundanceTest(infIDX,:)
%import data via UI

pIDs = table2array(proteinAbundanceTest(:,1));
pData = table2array(proteinAbundanceTest(:,2));
%check biomass sum function is working correctly
[X,P,C,R,D,L] = sumBioMass(ecModel)


[ecModel_batch,OptSigma] = getConstrainedModel(ecModel,c_source,sigma,Ptot,gR_exp);
%disp(['Sigma factor (fitted for growth on glucose): ' num2str(OptSigma)])
    %Get f (estimated mass fraction of enzymes in model)
    [f,~] = measureAbundance(ecModel.enzymes);
    %Change media to batch conditions:
    cd ../kcat_sensitivity_analysis
    ecModel = changeMedia_batch(ecModel,c_source);
    cd ../limit_proteins
    %Get a preliminary enzyme constrained model for performing the Kcats
    %sensitivity analysis
    %constrainEnzymes(model,Ptot,sigma,f,GAM,pIDs,data,gRate,c_UptakeExp,c_source)
    %[ecModel_batch,~,~] = constrainEnzymes(ecModel,Ptot,sigma,f);
    %[ecModel_batch,~,~] = constrainEnzymes(ecModel,Ptot,sigma,0.4,GAM);
    %%contrain model without experiemntal data
    % constrain model with experimental data
    [ecModel_batch_exp,~,~] = constrainEnzymes(ecModel,Ptot,sigma,0.4,GAM,pIDs,pData);
    
    
    
	solution            = solveLP(ecModel_batch,1);
    if ~isempty(solution.f)
        %Set the media according to the experimental conditions
        cd ../kcat_sensitivity_analysis
        ObjIndex = find(ecModel_batch.c);
        % If the model is overconstrained
        if (gRate-solution.x(ObjIndex))>0 
            fprintf('\n')
            disp('***************************************************************')
            disp('                The ECmodel is overconstrained                 ')
            %Perform a sensitivity analysis on the objective function with 
            %respect to the individual Kcat coefficients, the algorithm will 
            %iterate replacing the top limiting value according to the maximum 
            %value available in BRENDA for the same EC number until the objective
            %is no longer underpredicted 
            ecModel_batch = modifyKcats(ecModel_batch,gRate,modifications,name);
        else
            fprintf('\n')
            disp('***************************************************************')
            disp('              The ECmodel is not overconstrained               ')
        end    
        %The sigma factor is reffited for the specified conditions (constraints in the model)
        fprintf('\n')
        disp('***************************************************************')
        disp('        Fitting the average enzymes saturation factor          ')
        OptSigma          = sigmaFitter(ecModel_batch,Ptot,gRate,f);
        enzymePos         = strcmp(ecModel_batch.rxns,'prot_pool_exchange');
        currentEnzymeUB   = ecModel_batch.ub(enzymePos);
        ecModel_batch     = setParam(ecModel_batch,'ub','prot_pool_exchange', ...
                                     currentEnzymeUB*OptSigma/sigma);
        
        %Simulate growth on minimal media and export the top ten used 
        %enzymes to the file "topUsedEnzymes.txt" in the containing folder
        solution          = solveLP(ecModel_batch,1);
        topUsedEnzymes(solution.x,ecModel_batch,{'Min_glucose'},name);
        cd ../limit_proteins
    else
        disp('ecModel with enzymes pool constraint is not feasible')
    end

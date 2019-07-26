function [modelBMSplit] = splitBM(model,bmReaction)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here
bmIDX = find(strcmp('Biomass_Chlamy_mixo',creMixo.rxns))
%bmEductIDX = find(creMixo.S(:,bmIDX) < 0)
%bmProductIDX = find(creMixo.S(:,bmIDX) > 0)
bmStoich = nonzeros(creMixo.S(:,bmIDX))

bmMets =findMetsFromRxns(creMixo,bmIDX )
%check if comps and reaction names and order are equal
isequal(bmMets,comps(:,1))

cMetIDX = find(strcmp('C',comps(:,3)));
lMetIDX = find(strcmp('L',comps(:,3)));
pMetIDX = find(strcmp('P',comps(:,3)));
dMetIDX = find(strcmp('D',comps(:,3)));
rMetIDX = find(strcmp('R',comps(:,3)));
nMetIDX = find(strcmp('N',comps(:,3)));

cMet = cat(1,comps(cMetIDX,1), { 'carbon_pseudo[c]'})
lMet = cat(1,comps(lMetIDX,1), { 'lipid_pseudo[c]'} )
pMet = cat(1,comps(pMetIDX,1), { 'protein_pseudo[c]'} )
dMet = cat(1,comps(dMetIDX,1), { 'dna_pseudo[c]'} )
rMet = cat(1,comps(rMetIDX,1), { 'rna_pseudo[c]'} )
nMet = cat(1,comps(nMetIDX,1), { 'carbon_pseudo[c]'}, { 'lipid_pseudo[c]'}, { 'protein_pseudo[c]'}, { 'dna_pseudo[c]'}, { 'rna_pseudo[c]'} )

cStoi = cat(1,bmStoich(cMetIDX), 1.0)
lStoi = cat(1,bmStoich(lMetIDX), 1.0)
pStoi = cat(1,bmStoich(pMetIDX), 1.0)
dStoi = cat(1,bmStoich(dMetIDX), 1.0)
rStoi = cat(1,bmStoich(rMetIDX), 1.0)
nStoi = cat(1,bmStoich(nMetIDX), -1.0, -1.0, -1.0, -1.0, -1.0)


creMixoMod = addReaction(creMixo,'carbohydrate pseudoreaction','metaboliteList',cMet,'stoichCoeffList',cStoi, 'reversible',false);
creMixoMod = addReaction(creMixoMod,'lipid pseudoreaction','metaboliteList',lMet,'stoichCoeffList',lStoi, 'reversible',false);
creMixoMod = addReaction(creMixoMod,'protein pseudoreaction','metaboliteList',pMet,'stoichCoeffList',pStoi, 'reversible',false);
creMixoMod = addReaction(creMixoMod,'dna pseudoreaction','metaboliteList',dMet,'stoichCoeffList',dStoi, 'reversible',false);
creMixoMod = addReaction(creMixoMod,'rna pseudoreaction','metaboliteList',rMet,'stoichCoeffList',rStoi, 'reversible',false);
creMixoMod = addReaction(creMixoMod,'biomass pseudoreaction','metaboliteList',nMet,'stoichCoeffList',nStoi, 'reversible',false);

creMixoMod = changeObjective(creMixoMod, 'biomass pseudoreaction')
normGrowth = optimizeCbModel(creMixo);
modGrowth = optimizeCbModel(creMixoMod);

normGrowth.f 
modGrowth.f

modelBMSplit = inputArg2;
end


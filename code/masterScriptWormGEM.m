%
% FILE NAME:    masterScriptWormGEM.m
%
%
% PURPOSE: This script is for reconstruction of the Worm-GEM, by using
%          the Human-GEM as template and taking into account specie-specific
%          pathways/reactions.
%
%


%% Load Human-GEM as template
load('Human-GEM.mat');

% convert gene identifiers from Ensembl ids to gene symbols
[grRules,genes,rxnGeneMat] = translateGrRules(ihuman.grRules,'Name','ENSG');
ihuman.grRules    = grRules;
ihuman.genes      = genes;
ihuman.rxnGeneMat = rxnGeneMat;



%% Use MA reactions identifiers 

% load reaction annotaiton files
rxnAssoc = jsondecode(fileread('humanGEMRxnAssoc.JSON'));

%replace reaction identifiers with MA ids if available
ind = getNonEmptyList(rxnAssoc.rxnMAID);
ihuman.rxns(ind) = rxnAssoc.rxnMAID(ind);



%% Generate Worm-GEM by using Human-GEM as template

% get ortholog pairs from human to worm
wormOrthologPairs = extractAllianceGenomeOrthologs('human2WormOrthologs.json');
wormGEM = getModelFromOrthology(ihuman, wormOrthologPairs);



%% Incorporate species-specific reactions

% get metabolic networks based on the KEGG annoation using RAVEN function
KEGG_human=getKEGGModelForOrganism('hsa');
KEGG_worm=getKEGGModelForOrganism('cel');

% remove reactions shared with human
WormSpecificRxns=setdiff(KEGG_worm.rxns,KEGG_human.rxns);

% remove reactions included in Human-GEM
WormSpecificRxns=setdiff(WormSpecificRxns,rxnAssoc.rxnKEGGID);


% get species-specific network for manual inspection and then
% organize species-specific pathways into two tsv files:
wormSpecificNetwork=removeReactions(KEGG_worm,...
    setdiff(KEGG_worm.rxns,WormSpecificRxns), true, true, true);

% "wormSpecificMets.tsv" contains species-specific metabolites
metsToAdd = importTsvFile('wormSpecificMets.tsv');

% "wormSpecificRxns.tsv" contains species-specific reactions
rxnsToAdd = importTsvFile('wormSpecificRxns.tsv');
rxnsToAdd.subSystems = cellfun(@(s) {{s}}, rxnsToAdd.subSystems);

% integrate worm-specific metabolic network
[wormGEM, modelChanges] = addMetabolicNetwork(wormGEM, rxnsToAdd, metsToAdd);


%% Gap-filling for biomass formation
[wormGEM, gapfillNetwork]=gapfill4EssentialTasks(wormGEM,ihuman);
% Added 36 reactions for gap-filling


%% Save the model into mat, yml, and xml

wormGEM.id = 'Worm-GEM';
save('../model/Worm-GEM.mat', 'wormGEM');
writeHumanYaml(wormGEM, '../model/Worm-GEM.yml');
exportModel(wormGEM, '../model/Worm-GEM.xml');


from .calculateGMCSwithCycles import getGPRDict, stackLayersWithCycles, checkRegulatoryNetwork, bidict
from .MinimalCutSetClass import tqdm
from .cobraTools import re, initiateLogger
from .ProblemDefinitions import readjustmentProblem
from copy import deepcopy
from datetime import datetime
from .Utilities import createLog
from .calculateGMCSwithCycles import boolTransformer, mpbn

def calculateReadjustment(gMISList, model, regulatory_dataframe, num_layers=2, solver='gurobi', **kwargs):
    # Name of the analysis
    name = (
        "gMCSpy_Readjustment_"
        + model.id
        + "_"
        + datetime.now().strftime("%H-%M-%S--%d-%m-%Y")
    )
    
    # Save all the parameters in a log file    
    pathN = "logs/" + name + "/"
    path = kwargs.get('path', pathN)
    
    createLog(path)
    configName = ("configuration")
    logging = initiateLogger(configName + ".log", path + configName + ".log")
    # Save the parameters in the log file
    logging.info(f"solver: {solver}, num_layers: {num_layers}, solver: {solver}")
    logging.info(f"Model: {model.id}")
    logging.info(f"gMIS to check: {gMISList}")
    # close logger
    handler = logging.handlers[0]
    handler.close()
    logging.removeHandler(handler)
    
    
    ### Results logger
    
    if kwargs.get('path', None) is None:
        path = "logs/" + name
        name = ("results")
        handler = name + ".log"
        createLog(path + "/results/" + name + ".log")
        logging = initiateLogger(name + ".log", path + "/results/" + name + ".log")
    else:
        name = ("results.log")
        createLog(path + '/' + name)
        logging = initiateLogger(name, path + name)
        
    kwargs.pop('path', None)
    
    # Prepare the regulatory network
    regulatory_dict = regulatory_dataframe.to_dict('list')
    # Check if the regulatory network has valid sources and targets
    regulatory_dict, referenceDict = checkRegulatoryNetwork(regulatory_dict)
    # create a dictionary with the gpr rules associated to each reaction to get a unique list of gprs
    gprDict = getGPRDict(model)
    
    # Expand all the rules to the desired number of layers
    rulesDict = {}
    i = 0
    for key, value in tqdm(gprDict.items()):
        GPR = key, value
        rules, genes = calculateRules(GPR, regulatory_dict, num_layers, solver, **kwargs)
        rulesDict[i] = {'genes': genes, 'rules': rules}
        i += 1
    
    print('checking gMIS')
    for gMIS in tqdm(gMISList):
        if len(gMIS) == 1:
            continue
        calculateReadjustmentSingleGMIs(gMIS, rulesDict, referenceDict, logging, **kwargs)


def calculateReadjustmentSingleGMIs(gMIS, rulesDict, referenceDict, logging=None, **kwargs):        
        precomputedDict = {}
        # check if the gMIS is in the reference dictionary
        gMISDict = bidict()
        encodedgMIS = []
        for gene in gMIS:
            if '_KO' in gene:
                geneT = gene.split('_KO')[0]
            else:
                geneT = gene.split('_KI')[0]
            if geneT in referenceDict.values():
                encodedgMIS.append(referenceDict.inverse[geneT])
                gMISDict[gene] = referenceDict.inverse[geneT]
            else:
                encodedgMIS.append(geneT)
                gMISDict[gene] = geneT
        
        filteredRules = {}
        for gene in encodedgMIS:
            for key, value in rulesDict.items():
                ruleGenes = value['genes']
                if set([gene]).issubset(set(ruleGenes)):
                    filteredRules[key] = value
        
        #print(gMIS)
        for key, value in filteredRules.items():
            results = calculateReadjustmentInRule(key, value['rules'], gMISDict, precomputedDict)
            if results is not None:
                for res in results:
                    if res['readjustment']:
                        logging.info(res)
                        break
                    

            
def calculateRules(GPR, regulatory_dict, num_layers=2, solver='gurobi', **kwargs):
    key, value = GPR
    expandedRules = stackLayersWithCycles(value['genes'], regulatory_dict, num_layers, solver)
    boolform = boolTransformer(expandedRules, value)
    genes = list(boolform.keys())
    return boolform, genes

def calculateReadjustmentInRule(keyGPR, rules, gMISDict, precomputedDict):
    readaptList = []    
    for key, items in gMISDict.items():
        # First check if the "problem" gene is in the rules
        if items not in rules:
            continue        
        # Check if the gene is knocked out or knocked in, to define intervention
        interaction = 0 if key.endswith('_KO') else 1
        # Get the acompanying genes
        provitionalList = list(gMISDict.keys())
        # Remove the gene "problem" from the list
        provitionalList.remove(key)
        # Calculate the value of the set without the "problem" gene
        val = calculateSetValues(provitionalList, rules)
        # If val is none, only the "problem" gene is in the set, continue
        if val is None:
            continue
        # Perform the knockout to the target gene
        network, reachable_from = boolean_Knockout(items, interaction, rules)
        # Prepare the constraints to obtain the desire attractors
        constraints = {}
        for gene in provitionalList:
            if gene.endswith('_KI'):
                gene_without_lastname = gene.split('_KI')[0]
                tInt = 1
            else:
                gene_without_lastname = gene.split('_KO')[0]
                tInt = 0
            constraints[gene_without_lastname] = tInt
        # Implement a cache strategy to avoid recalculating the same attractors
        problemGeneDef = f"i_{items}_{interaction}"
        identifierSet = set([problemGeneDef, keyGPR])
        for key, value in constraints.items():
            identifierSet.add(f"c_{key}_{value}")
            
        identifier = frozenset(identifierSet)
        if identifier in precomputedDict:
            pass
        else:
            act = network.attractors(reachable_from=reachable_from, constraints=constraints, limit=1)
            precomputedDict[identifier] = list(act)
        
        try:
            attractors = precomputedDict[identifier][0]
        except:
            readaptList.append({'readjustment': True, 'gene': key, 'gMISDict': list(gMISDict.keys()), 'key': keyGPR})
        
    return readaptList
    '''
    
    names = bidict()
    resultsDict = {}
    for name in gMISDict.keys():
        gene_without_lastname = name.split('_KI')[0] if name.endswith('_KI') else name.split('_KO')[0]
        names[name] = gene_without_lastname
        resultsDict[name] = []
        break
        # Check if configuration exists
    
    try:
        attractors = list(act)[0]
    except:
        gMIS = list(gMISDict.keys())
        return {'readjustment': True, 'gMIS': gMIS, 'constraints': constraints}
    '''

def calculateSetValues(provitionalList, rules):
    val = None
    for gene in provitionalList:
        if gene.endswith('_KI'):
            gene_without_lastname = gene.split('_KI')[0]
            tInt = 1
        else:
            gene_without_lastname = gene.split('_KO')[0]
            tInt = 0
        if gene_without_lastname in rules:
            if val is None:
                val = tInt
            else:
                val = val + tInt
    return val
    
def boolean_Knockout(target:str, interventionVal:int, rules:dict):
    '''
    This function will create a new set of rules with the target gene knocked out
    '''
    newRules = deepcopy(rules)
    newRules.pop(target)
    # Creating a boolean network
    network = mpbn.MPBooleanNetwork(newRules)
    # Setting the target to 0
    reachable_from = {target: interventionVal}
    
    return network, reachable_from
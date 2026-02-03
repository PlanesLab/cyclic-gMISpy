from .ProblemDefinitions import cobra
from .ProblemInterpreters import np

from .OptimizationProblem import OptimizationProblem

from .ProblemInterpreters import beartype
from .ProblemDefinitions import List
from .ProblemDefinitions import scipy
from .ProblemDefinitions import buildDictionaryDossageNetwork, buildDictionaryMCSGeneProblem


from .calculateMCS import calculateMCS

from .cobraTools import createSparseMatrix, parseGPRToModel, simplifyGMatrix, relatedRows, mergeIsforms, initiateLogger

from collections import defaultdict
import warnings
import time

from .MinimalCutSetClass import MinimalCutSetProblem
from .MinimalCutSetClass import uuid
from .MinimalCutSetClass import tqdm
from .MinimalCutSetClass import logging
from .MinimalCutSetClass import prepareModelNutrientGeneMCS

from datetime import datetime
import pandas as pd
from .Utilities import createLog, re

from bidict import bidict
import concurrent.futures
import subprocess
import mpbn
import sys
import os
from contextlib import contextmanager
from copy import deepcopy
import shutil
from pathlib import Path
import networkx as nx


@contextmanager
def suppress_output():
    # Save current stdout and stderr
    stdout = sys.stdout
    stderr = sys.stderr
    # Redirect stdout and stderr to devnull (a special file that discards all data written to it)
    sys.stdout = open(os.devnull, 'w')
    sys.stderr = open(os.devnull, 'w')
    try:
        yield
    finally:
        # Restore stdout and stderr
        sys.stdout = stdout
        sys.stderr = stderr

def calculateParallelGMIS(cobraModel, regulatory_dataframe, num_layers=2, **kwargs):
    '''
    '''
    # Name of the analysis
    name = (
        "gMCSpy_GMIS_"
        + cobraModel.id
        + "_"
        + datetime.now().strftime("%H-%M--%d_%m_%Y")
    )
    
    # Prefered solver is gurobi due to performance experience during benchmarking
    solver = kwargs.get("solver", "gurobi")
    maxKOLength = kwargs.get('maxKOLength', 3)
    timeLimit = kwargs.get("timeLimit", 1e4)
    numWorkers = kwargs.get("numWorkers", 1)
    maxNumberGMCS = kwargs.get('maxNumberGMCS', 1e6)
    verbose = kwargs.get('verbose', 0)
    earlyStop = kwargs.get("earlyStop", True)
    removeRelated = kwargs.get("removeRelated", False)
    isoformSeparator = kwargs.get("isoformSeparator", None)
    targetKOs = kwargs.get("targetKOs", None)
    geneSubset = kwargs.get("geneSubset", None)
    isNutrient = kwargs.get("isNutrient", False)
    exchangeReactions = kwargs.get("exchangeReactions", None)
    saveGMatrix = kwargs.get("saveGMatrix", False)
    checkSolutions = kwargs.get("checkSolutions", False)
    saveSolutions =  kwargs.get('saveSolutions', True)
    onlyGenes = kwargs.get('onlyGenes', False)
    if kwargs.get("forceLength") is None:
        if solver == "gurobi":
            forceLength = False
        elif solver == "cplex":
            forceLength = True
        elif solver == "scip":
            forceLength = True   
    else:
        forceLength = kwargs.get("forceLength", False)
        
        
    # Save all the parameters in a log file    
    path = "logs/" + name + "/"
    createLog(path)
    configName = ("configuration")
    logging = initiateLogger(configName + ".log", path + configName + ".log")
    # Save the parameters in the log file
    logging.info(f"solver: {solver}, maxKOLength: {maxKOLength}, timeLimit: {timeLimit}, numWorkers: {numWorkers}, maxNumberGMCS: {maxNumberGMCS}, num_layers: {num_layers}, verbose: {verbose}, forceLength: {forceLength}, earlyStop: {earlyStop}, removeRelated: {removeRelated}, isoformSeparator: {isoformSeparator}, targetKOs: {targetKOs}, geneSubset: {geneSubset}, isNutrient: {isNutrient}, exchangeReactions: {exchangeReactions}, saveGMatrix: {saveGMatrix}, checkSolutions: {checkSolutions}, saveSolutions: {saveSolutions}, forceLength: {forceLength}, onlyGenes: {onlyGenes}")
    logging.info(f"Model: {cobraModel.id}")
    # close logger
    handler = logging.handlers[0]
    handler.close()
    logging.removeHandler(handler)
    
    # Pop the parameters that are not used downstream
    kwargs.pop("maxKOLength", None)
    kwargs.pop("maxNumberGMCS", None)
    kwargs.pop("verbose", None)
    kwargs.pop("earlyStop", None)
    kwargs.pop("removeRelated", None)
    kwargs.pop("isoformSeparator", None)
    kwargs.pop("targetKOs", None)
    kwargs.pop("geneSubset", None)
    kwargs.pop("isNutrient", None)
    kwargs.pop("exchangeReactions", None)
    kwargs.pop("saveGMatrix", None)
    kwargs.pop("checkMCS", None)
    kwargs.pop('saveSolutions', None)
    kwargs.pop("forceLength", None)
    kwargs.pop("timeLimit", None)
    kwargs.pop("OnlyGenes", None)
    path = "logs/" + name + "/process/"
    createLog(path)
    path = "logs/" + name
    name = ("process")
    handler = name + ".log"
    logging = initiateLogger(name + ".log", path + "/process/" + name + ".log")
    allTime = time.time()

    
       
    if isNutrient:
        cobraModel = prepareModelNutrientGeneMCS(cobraModel, exchangeReactions)

    # Prepare the regulatory network
    regulatory_dict = regulatory_dataframe.to_dict('list')
    # Check if the regulatory network has valid sources and targets
    regulatory_dict, referenceDict = checkRegulatoryNetwork(regulatory_dict)
    # Check if the regulatory network has the correct format, should have source_ENSEMBL, target_ENSEMBL, and interaction
    expectedColumns = ['source_ENSEMBL', 'target_ENSEMBL', 'interaction']
    for column in expectedColumns:
        if column not in regulatory_dict:
            raise ValueError(f"Expected column {column} not found in the regulatory network")
    
    # Filter the regulatory network to keep only the relationships that are genes
    if onlyGenes:
        # Drop the rows that are not genes
        regulatory_dataframe = regulatory_dataframe[regulatory_dataframe['source_ENSEMBL'].str.contains('ENSG')]
        regulatory_dataframe = regulatory_dataframe[regulatory_dataframe['target_ENSEMBL'].str.contains('ENSG')]
    
    #print(regulatory_dataframe.shape)    
    

    startTime = time.time()
    gObject = calculateRegNetGMatrix(cobraModel, regulatory_dict, num_layers, solver=solver, maxKOLength=maxKOLength, path=path, onlyGenes=onlyGenes, numWorkers=numWorkers)
    endTime = time.time()
    logging.info("GMatrix: " + str(endTime - startTime))

    
    gDict = gObject["gDict"]
    gMatrix = gObject["gMatrix"]
    relationships = gObject["relationships"]
    numberNewGenesByKO = gObject["numberNewGenesByKO"]
    
    
    # Filter the G matrix to keep only genes that are in the geneSubset
    if geneSubset:
        #Make a set of the geneSubset
        geneSet = set(geneSubset)
        # Extract all possible interventions
        gDictKeys = list(gDict.keys())
        # Remove the interventions that are not in the geneSubset
        for key in gDictKeys:
            if not key.issubset(geneSet):
                gDict.pop(key)

        gMatrix = createSparseMatrix(gDict, cobraModel.reactions)
        [relationships, numberNewGenesByKO] = relatedRows(gDict, mergeIsforms(isoformSeparator))
   
    if saveGMatrix:        
        gMatrix_Strings = []
        for ko, reactions  in gDict.items():
            string_ko = "@@".join(list(ko))
            string_reactions = "@@".join(cobraModel.reactions[x].id for x in reactions)
            row = string_ko + " -> " + string_reactions
            gMatrix_Strings.append(row)
        filename = f"Gmat_python_{cobraModel.id}_{solver}.csv"
        gPath = path + "/GMatrix/" 
        createLog(gPath)
        completePath = gPath + filename
        df = pd.DataFrame({'gMatrix': gMatrix_Strings})
        df.to_csv(completePath)


    # Initialize the genes to be knocked out, e.g. genes that you want to be part of the gMCS
    if targetKOs is None:
        targetKOs = []

    # Make sure that targetKOs is a list
    if not isinstance(targetKOs, list):
        raise ValueError("targetKOs must be a list")

    # Initialize the subset of genes to be studied , e.g. the space of genes to calculate the gMCS
    if geneSubset is None:
        geneSubset = []

    # Make sure that geneSubset is a list
    if not isinstance(geneSubset, list):
        raise ValueError("geneSubset must be a list")

    # List of possible interventions
    gDictKeys = list(gDict.keys())

    startTime = time.time()
    problem = buildDictionaryMCSGeneProblem(
        cobraModel,
        gDict,
        gMatrix,
        relationships=relationships,
        maxKOLength=maxKOLength,
        forceLength=forceLength,
        numberNewGenesByKO=numberNewGenesByKO,
        verbose=verbose,
        logPath=path
    )
    endTime = time.time()
    logging.info("BuildProblem: " + str(endTime - startTime))

        
    gurobi_default_params = {
        "MIPFocus": 3,
        "TimeLimit": max(1, timeLimit),
        "Threads": numWorkers,
        "PoolSolutions": 10000,
        "PoolGap": 0.01,
        "PoolSearchMode": 1,
        "Cutoff": maxKOLength,
        "Presolve": 1,
        "PreSOS1Encoding": 2,
        "Cuts": 2,
        #"OptimalityTol": 1e-9,
    }

    if forceLength:
        gurobi_default_params['MIPFocus'] = 3
        gurobi_default_params['TimeLimit'] = 500

    problem.setCharacteristicParameter(maxKOLength)
    problem.setModelName(cobraModel.id)
    problem.setLogPath(path)
    # Set the solver
    
    problem.setSolver(solver)

    cplex_default_params = {
        "mip.tolerances.integrality": 1e-5,
        "emphasis.mip": 4,
        "timelimit": max(10, timeLimit),
        "threads": numWorkers,
        "preprocessing.presolve": 1,
        "mip.limits.populate": 1000,
        "mip.pool.relgap": 0.1,
        "mip.tolerances.uppercutoff": maxKOLength,
        "mip.pool.intensity": 3,
        "mip.strategy.probe": 2,
        "mip.strategy.variableselect": 4,
        "preprocessing.sos1reform": 1,
    }


    if solver == "gurobi":
        params = gurobi_default_params
    elif solver == "cplex":
        params = cplex_default_params
    else:
        params = {}
    # Interpret the problem, depending on the solver each problem has to be interpreted differently into the solver's interface
    optimization = problem.interpretProblem(verbose=2, parameters=params)

    # Import the common problem methods
    problemMethods = MinimalCutSetProblem()

    # Get the objective variables, they are the variables that we are interested in
    knockouts = list(gDict.keys())
    zpVariables = [
        key for key, data in problem.getVariables().items() if "zp_" in data["name"]
    ]
    solutionVariables = {"indices": zpVariables, "names": knockouts}
    # Set the bounds for the solver, e.g. maximum number of solutions, maximum numbers in a solution, or if we want to force the step by step calculation etc.
    bounds = {
        "MAX_SOLUTION_LENGTH": maxKOLength,
        "MAX_SOLUTIONS": maxNumberGMCS,
        "forceLength": forceLength,
    }

    # Solve the problem iteratively
    solutionDict = problemMethods.iterateOverSolver(
        problem=problem,
        problemInterpreted=optimization,
        solver=solver,
        solutionVariables=solutionVariables,
        bounds=bounds,
        mergeSolutions=True,
        earlyStop=earlyStop,
        removeRelated=removeRelated,
        verbose=2,
        handler=handler,
        iterationProgress=True,
    )
    allEndTime = time.time()

    logging.info("totalTime: " + str(allEndTime - allTime))
        
    if saveSolutions:
        solutionPath = path + '/solutions/'
        solutionFileName = 'solutions.log'
        createLog(solutionPath)
        loggingSolutions = initiateLogger(solutionFileName, solutionPath + solutionFileName)
        loggingSolutions.info('gene,reactions')
        geneSolutionList = []
        for order, ko in solutionDict.items():
            reactions_log_list = []
            for gene in ko['solution']:
                reactions = gDict[frozenset([gene])]
                reactions_log_list.append([cobraModel.reactions[x].id for x in reactions])
            reactions = ','.join(list(set().union(*reactions_log_list)))
            ko_log = ','.join(list(ko['solution']))  
            ko_log_decode = []
            for gene in ko['solution']:
                  if gene.startswith('M_I_'):
                      ko_log_decode.append(decode_string(gene, referenceDict))
                  else:
                      ko_log_decode.append(gene)                
            loggingSolutions.info(f'{ko_log_decode}:[{reactions}]')

        solutionHandler = loggingSolutions.handlers[0]
        loggingSolutions.removeHandler(solutionHandler)
        solutionHandler.close()
               
    handler = logging.handlers[0]
    logging.removeHandler(handler)
    handler.close()
    
    return solutionDict
    

def analysisGPRs(model:cobra.Model,
                 regulatory_dict: dict, 
                 gDict: dict, 
                 num_layers: int,
                 path: str, 
                 solver='gurobi', 
                 numWorkers=1, 
                 maxKOLength=3,
                 gprDict=None, 
                 listProblemGPRs=None):
    if gprDict is None: 
        gprDict = getGPRDict(model)
  
    # some times inconsistencies in the regulatory network can cause unreprogrammable gprs, and we need to reduce the number of layers
    # This check is useful when we have a list of gprs that we want to re analyze. 
    checkFlag = False
    if listProblemGPRs is not None:
        print(f"Analyzing again the gpr {gprDict[listProblemGPRs[0]]['id']} with the reduced number of layers {num_layers}")
        checkFlag =True
    parallelGPRDict = []
    expandedRulesDict = {}
    
    
    for key, value in tqdm(gprDict.items(), disable=checkFlag):
        if checkFlag:
            if key not in listProblemGPRs:
                continue
        expandedRulesWithCycles = stackLayersWithCycles(value['genes'], regulatory_dict, num_layers)        
        if len(expandedRulesWithCycles['sources']) == 0:
            regulation = False
            sourceGenes = []
            targetGenes = []
            # update the gprDict with the new information
            gprDict[key]['regulation'] = regulation
            gprDict[key]['source'] = sourceGenes
            gprDict[key]['target'] = targetGenes
            gprDict[key]['numLayers'] = num_layers
            
            calculateTraditionalMCS(gDict, value['gpr'], value['reactions'], solver=solver, numWorkers=numWorkers, maxKOLength=maxKOLength)  
        else:
            regulation = True
            sourceGenes = set(expandedRulesWithCycles['sources'])
            targetGenes = set(expandedRulesWithCycles['targets'])
            gprDict[key]['regulation'] = regulation
            gprDict[key]['source'] = sourceGenes
            gprDict[key]['target'] = targetGenes
            gprDict[key]['numLayers'] = num_layers
            parallelGPRDict.append(key)
            expandedRulesDict[key] = expandedRulesWithCycles
    
    return parallelGPRDict, expandedRulesDict, gprDict

# save the expanded rules and the GPRs
def preProcessing(model: cobra.Model,
            regulatory_network: dict,
            layers: int, 
            maxLen: int, 
            gDict: dict, 
            solver: dict, 
            numWorkers: int, 
            path: str=None,
            processLogger: logging.Logger=None,
            gprDict=None, 
            listProblemGPRs=None):
    if path is None:
        path = "."
    # create a new folder inside bnAnalysis 
    os.makedirs(f'{path}/bnAnalysis/', exist_ok=True)
    os.makedirs(f'{path}/bnAnalysis/bnnets', exist_ok=True)
    
    path = f'{path}/bnAnalysis/'

       
    parallelGPRList, expandedRulesDict, gprDict = analysisGPRs(model=model,
                                                      regulatory_dict=regulatory_network,
                                                      gDict=gDict,
                                                      num_layers=layers,
                                                      path=path,
                                                      solver=solver,
                                                      numWorkers=numWorkers,
                                                      maxKOLength=maxLen,
                                                      gprDict=gprDict, 
                                                      listProblemGPRs=listProblemGPRs)
    
    task_params = []
    for key in parallelGPRList:
        task_params.append((key, gprDict[key], expandedRulesDict, path))
    
    # Start time
    start = time.time()

    for parms in task_params:        
        key, value, expandedRulesDict, path = parms
        id = value['id']
        expandedRulesWithCycles = expandedRulesDict[key]
        boolForm = boolTransformer(expandedRulesWithCycles, value)
        filenamePathBn = f'{path}/bnnets/{id}-{layers}.bnet'
        boolSaver(boolForm, filenamePathBn)

    # End time
    end = time.time()
    if processLogger is not None:
        processLogger.info(f"Time to save the bnet files: {end - start} seconds")
    
    return parallelGPRList, gprDict

def find_complex_cycles(data_string):
    nodes_in_cycles = set()
    G = nx.DiGraph()
    
    # --- Parsing ---
    lines = data_string.strip().split('\n')
    for line in lines:
        if not line.strip(): continue
        
        parts = line.split(',')
        target = parts[0].strip()
        
        if len(parts) > 1:
            # The right side contains the regulators (Sources)
            regulators = parts[1].strip().split('|')
            
            for source in regulators:
                source = source.strip()
                # --- CRITICAL STEP: Ignore Self-Loops ---
                if source and source != target:
                    G.add_edge(source, target)

    # --- Detection ---
    # nx.simple_cycles finds all elementary cycles
    # Since we didn't add self-loops to G, all found cycles will have length >= 2
    cycles = list(nx.simple_cycles(G))
    
    # --- Output ---
    if not cycles:
        return list(nodes_in_cycles)  # No cycles found
    else:
        #print(f"Found {len(cycles)} complex cycle(s):")
        for i, cycle in enumerate(cycles, 1):
            # Format the output to look like a flow: A -> B -> C -> A
            cycle_path = " -> ".join(cycle) + " -> " + cycle[0]
            nodes_in_cycles.update(cycle)
    
    return list(nodes_in_cycles)

def calculateRegNetGMatrix(model, regulatory_dict, num_layers=2, **kwargs):  
    global maxKOLength  
    isoformSeparator = kwargs.get('isoformSeparator', None)
    maxKOLength = kwargs.get('maxKOLength', 3)
    dictManager = mergeIsforms(isoformSeparator)
    path = kwargs.get('path', None)
    onlyGenes = kwargs.get('onlyGenes', False)
    solver = kwargs.get('solver', 'gurobi')
    numWorkers = kwargs.get('numWorkers', 1)  
    processLogger = logging.getLogger('process.log')
    modelReactions = [reaction.id for reaction in model.reactions]
    gDict = {}
    print("Number of workers: ", numWorkers)
    # stepOne
    parallelGPRList, gprDict = preProcessing(model=model,
                              regulatory_network=regulatory_dict,
                              layers=num_layers,
                              maxLen=maxKOLength,
                              gDict=gDict,
                              solver=solver,
                              numWorkers=numWorkers,
                              path=path,
                              processLogger=processLogger)
    
    targert_path = Path(path, 'bnAnalysis/bnnets')
    resultDict = parallelAnalysisBonesis(targert_path, numWorkers=numWorkers, timeLimitBonesis=10*60,
                              processLogger=processLogger)
    
    # transform parallelGPRDict into a dictionary with key the enumeration of the gprDict and keep the same value
    enumGPRDict = {i: v for i, v in enumerate(gprDict.items()) if v[0] in parallelGPRList}
        
    targert_path = Path(path, 'bnAnalysis/completed_bnnets')
    cycles_in_enumGPRDict = {}
    for gpr, value in enumGPRDict.items():
        bnet_file = targert_path / f"{gpr}-{num_layers}.bnet"
        nodes = find_complex_cycles(bnet_file.read_text())
        if nodes:
            cycles_in_enumGPRDict[gpr] = nodes
    
    # insert the results from the bonesis analysis into the gDict
    for key, value in resultDict.items():
        if key in cycles_in_enumGPRDict:
            #print(f"Processing GPR: {key}, gpr with cycles in genes: {cycles_in_enumGPRDict[key]}")
            pass
        # if the length of the value is 1, we only have gpr=0 or empty
        if len(value) == 1:
            continue
        if value == 'empty':
            # if bonesis cannot find any solution, we calculate try to reduce the number of layers until we find a solution
            gpr = {enumGPRDict[key][0]:enumGPRDict[key][1]} 
            reactions = enumGPRDict[key][1]['reactions']
            gmis = calculate_gMIS_iterative(gDict, 
                                        model, 
                                        regulatory_dict, 
                                        path, 
                                        gpr, 
                                        gprDict,
                                        reactions, 
                                        num_layers, 
                                        numWorkers, 
                                        maxKOLength, 
                                        solver) 
                    # value is a list of dictionaries
            for resDict in gmis:
                geneSolutions = set()
                # if directly blocking the gpr is part of the solution, we skip it
                if "gpr" in list(resDict.keys()):
                    continue
                if key in cycles_in_enumGPRDict:
                    genes_in_cycle = [gene for gene in resDict.keys() if gene in cycles_in_enumGPRDict[key]] 
                    if len(genes_in_cycle) > 0:
                        print(f'cycle found in gpr {key}, genes in cycle: {genes_in_cycle}')
                        continue
                    for resDictKey, resDictVals in resDict.items():                        
                        lastName = None
                        if resDictVals == 0:
                            lastName = 'KO'
                        elif resDictVals == 1:
                            lastName = 'KI'
                        sol = resDictKey + '_' + lastName
                        geneSolutions.add(sol)
                    # jump
                    geneSolutions = frozenset(geneSolutions)
                    addToDict(gDict, geneSolutions, reactions)
                    
                else:    
                    for resDictKey, resDictVals in resDict.items():
                        lastName = None
                        if resDictVals == 0:
                            lastName = 'KO'
                        elif resDictVals == 1:
                            lastName = 'KI'
                        sol = resDictKey + '_' + lastName
                        geneSolutions.add(sol)
                    # jump
                    geneSolutions = frozenset(geneSolutions)
                    addToDict(gDict, geneSolutions, reactions)
            continue 
        # key is the enumeration of the gprDict
        reacts = enumGPRDict[key][1]['reactions']
        # value is a list of dictionaries
        for resDict in value:
            geneSolutions = set()
            # if directly blocking the gpr is part of the solution, we skip it
            if "gpr" in list(resDict.keys()):
                continue
            if key in cycles_in_enumGPRDict:
                genes_in_cycle = [gene for gene in resDict.keys() if gene in cycles_in_enumGPRDict[key]] 
                if len(genes_in_cycle) > 0:
                    print(f'cycle found in gpr {key}, genes in cycle: {genes_in_cycle}')
                    continue
                for resDictKey, resDictVals in resDict.items():                    
                    lastName = None
                    if resDictVals == 0:
                        lastName = 'KO'
                    elif resDictVals == 1:
                        lastName = 'KI'
                    sol = resDictKey + '_' + lastName
                    geneSolutions.add(sol)
                # jump
                geneSolutions = frozenset(geneSolutions)
                addToDict(gDict, geneSolutions, reacts)
            else:
                for resDictKey, resDictVals in resDict.items():
                    lastName = None
                    if resDictVals == 0:
                        lastName = 'KO'
                    elif resDictVals == 1:
                        lastName = 'KI'
                    sol = resDictKey + '_' + lastName
                    geneSolutions.add(sol)
                # jump
                geneSolutions = frozenset(geneSolutions)
                addToDict(gDict, geneSolutions, reacts)
    
    # save the gprDict
    nameGPRFile = f'{path}/GPRs2.txt'
    filename =  nameGPRFile
    # Write set to file, create file if it doesn't exist
    with open(filename, 'w') as file:
        file.write('reactions,genes,gpr,id,source,target,regulation,numLayers\n')
        for gprKey, gprVals in gprDict.items():
            #file.write(",".join(map(str,gprVals)))
            file.write(",".join(map(str,gprVals.values())))
            file.write("\n")
    
    # Make sure that each entry in the gDict is unique list of reactions
    for key, value in gDict.items():
        gDict[key] = list(set(value))
    
    # Filter out the interventions that are larger than the maxKOLength
    filteredGDict = gDict.copy()
    for key in gDict.keys():
        if len(key) > maxKOLength:
            filteredGDict.pop(key)
    
    if onlyGenes:
        # Filter out the interventions that are not composed of genes
        newDict = filteredGDict.copy()
        for key, items in filteredGDict.items():
            flag = False
            for gene in key:
                if 'ENSG' not in gene:
                    flag = True
            if flag:
                newDict.pop(key)
        filteredGDict = newDict
            
    
    newGDict = transformReactionIdsIntoIndexes(filteredGDict, modelReactions)
    
    simpG = simplifyGMatrix(newGDict, dictManager)
    
    [relationships, numberNewGenesByKO] = relatedRows(simpG, dictManager)
    
    gMatrix = createSparseMatrix(simpG, modelReactions)
    
    gObject = {
        "gMatrix": gMatrix,
        "gDict": simpG,
        "relationships": relationships,
        "numberNewGenesByKO": numberNewGenesByKO,
    }
    
    
    
    return gObject

def parallelAnalysisBonesis(target_path: Path, numWorkers=1, **kwargs):
    #Export clingo options
    os.environ['CLINGO_OPTS'] = "--sign-def=pos --vsids-acids"
    
    global timeLimitBonesis
    
    # create a new folder inside bnAnalysis 
    path = str(target_path.parent)
    os.makedirs(f'{path}/completed_bnnets', exist_ok=True)
    
    # get a list of all the files in the target path
    target_files = [f for f in target_path.iterdir() if f.is_file()]

    timeLimitBonesis = kwargs.get('timeLimitBonesis', 5*60)
    processLogger = kwargs.get('processLogger', None)
    
    # Time the execution of the subprocesses
    startime = time.time()
    
    dict_res = {}
    # while target_files is greater than 0
    while len(target_files) > 0:
        # order the files by name (0,1, ..., 101, 102, ... n)
        target_files = sorted(target_files, key=lambda x: int(x.name.split('-')[0]))

        # Using ThreadPoolExecutor to run subprocesses in parallel
        with concurrent.futures.ThreadPoolExecutor(max_workers=numWorkers) as executor:
            # Map each file to the calculate_gMIS function
            results = list(executor.map(calculate_gMIS, target_files))
        for res in results:
            if res is not None:
                dict_res.update(res)
        target_files = [f for f in target_path.iterdir() if f.is_file()]
        timeLimitBonesis += 180
    endtime = time.time()
    if processLogger is not None:
        processLogger.info(f"Total time to run the bonesis analysis: {endtime - startime}")
    return dict_res

# Function to process each target file
def calculate_gMIS(bnet):
    precessLogger = logging.getLogger('process')
    fileName = str(bnet)
    name = fileName.split("/")[-1].split("-")[0]
    try:
        result = subprocess.run(
            [
                "bonesis-reprogramming",
                bnet,
                '{"gpr":0}',
                str(maxKOLength)
            ],
            capture_output=True,
            text=True,
            timeout=timeLimitBonesis  # Timeout set to 5 minutes
        )
        
        shutil.move(fileName, fileName.replace('bnnets', 'completed_bnnets'))
        try:
            # Split the output by new lines
            lines = result.stdout.splitlines()
            # Convert the output to a list of dictionaries
            dict_list = [eval(line) for line in lines]
            empty_dict_list = [{}]
            # extract what bnt file was used from filename (x.bnet)
            fileName = fileName.replace("\\", "/")
            name = fileName.split("/")[-1].split("-")[0].replace(".bnet", "")
            name = int(name)
            # We need to calculate only the metabolic gene cut sets if bonesis cannot find any solution
            # When we encounter an empty dictionary, this means there is no possible reprogramming.
            if dict_list == empty_dict_list:
                return {name: "empty"}
            res = {name: dict_list}
            return res
        except:
            precessLogger.info("The output could not be converted to a list of dictionaries. (The result is either empty or null)")

    except subprocess.TimeoutExpired:
        print(f"The command for {bnet} took too long to run and was terminated. Will retry with {timeLimitBonesis + 180} seconds")
    except subprocess.CalledProcessError as e:
        print(f"The command for {bnet} encountered an error: {e}")

           
def calculate_gMIS_iterative(gDict, model, regulatory_dict, path, gpr, gprDict, reactions, num_layers, numWorkers, maxKOLength, solver):
    flag = True
    while flag:
        num_layers -= 1
        parallelGPRDict = preProcessing(model=model,
                                regulatory_network=regulatory_dict,
                                layers=num_layers,
                                maxLen=maxKOLength,
                                gDict=gDict,
                                numWorkers=numWorkers,
                                path=path, 
                                gprDict=gprDict,
                                listProblemGPRs=list(gpr.keys()), 
                                solver=solver)

        targert_path = Path(path, 'bnAnalysis/bnnets')
        resultDict = parallelAnalysisBonesis(targert_path, numWorkers=numWorkers, timeLimitBonesis=10*60,
                              processLogger=None)
        
        for key, value in resultDict.items():
            if value != 'empty':
                flag = False
                return value
             

def calculateTraditionalMCS(gDict, gpr, reactions, solver, numWorkers, maxKOLength):
    # get the genes from the gpr rule
    stringGPR = gpr.to_string()
    genes = [gene for gene in gpr.genes]
    if len(genes) == 1:
        # Rename the genes as they are KOs
        genes = [f'{gene}_KO' for gene in genes]
        key=frozenset(genes)
        addToDict(gDict, key, reactions)
    # Only and
    elif bool(re.search(' and ', stringGPR) and not re.search(' or ', stringGPR)):
        # Rename the genes as they are KOs
        genes = [f'{gene}_KO' for gene in genes]
        # Only one gene needs to be KO to stop the reaction
        for gene in genes:
            key=frozenset([gene])
            addToDict(gDict, key, reactions)

    # Only or
    elif bool(re.search(' or ', stringGPR) and not re.search(' and ', stringGPR)):
        # Rename the genes as they are KOs
        genes = [f'{gene}_KO' for gene in genes]
        # All genes need to be KO to stop the reaction
        key=frozenset(genes)
        addToDict(gDict, key, reactions)   # And and or
    else:
        treeGenes, model = parseGPRToModel(stringGPR, 'gpr', None)
        solutionDict = calculateMCS(model, MAX_LENGTH_MCS=maxKOLength, MAX_MCS=1e6, rxnSubset=list(treeGenes), solver=solver, numWorkers=numWorkers)
        for _, data in solutionDict.items():
            # Rename the genes as they are KOs
            genes = [f'{gene}_KO' for gene in data["solution"]]
            key = frozenset(genes)
            addToDict(gDict, key, reactions)
        
def addToDict(gDict, key, reactions):
    if key in gDict:
        gDict[key] = gDict[key] + (reactions)
    else:
        gDict[key] = reactions
            
def getGPRDict(model):
    gpr_dict = {}
    for reaction in model.reactions:
        gpr = reaction.gpr.to_string()
        if gpr in gpr_dict:
            gpr_dict[gpr]['reactions'].append(reaction.id)
        elif gpr == '':
            pass
        else:
            gpr_dict[gpr] = {'reactions': [reaction.id], 'genes': [gene.id for gene in reaction.genes], 'gpr': reaction.gpr}
    for i, gpr in enumerate(gpr_dict):
        gpr_dict[gpr]['id'] = i
        gpr_dict[gpr]['source'] = {}
        gpr_dict[gpr]['target'] = {}
        gpr_dict[gpr]['regulation'] = None
        
    return gpr_dict
    

def constructLayerWithCycles(genes, regulatory_dict):
    layer_dict = {}
    layer_dict['sources'] = []
    layer_dict['targets'] = []
    layer_dict['interaction'] = []
    layer_dict['genes'] = genes
    for gene in genes:
        if gene in regulatory_dict['target_ENSEMBL']:
            # get the index of all the target genes that match to gene
            target_indices = [i for i, x in enumerate(regulatory_dict['target_ENSEMBL']) if x == gene]
            [layer_dict['sources'].append(regulatory_dict['source_ENSEMBL'][i]) for i in target_indices]
            [layer_dict['targets'].append(regulatory_dict['target_ENSEMBL'][i]) for i in target_indices]
            [layer_dict['interaction'].append(regulatory_dict['interaction'][i]) for i in target_indices]
            layer_dict['genes'] = list(set(layer_dict['genes']).union(set(layer_dict['sources'])))
    return layer_dict
    
def stackLayersWithCycles(genes, regulatory_dict, num_layers):
    if num_layers == 0:
        return {'sources': [], 'targets': [], 'interaction': [], 'genes': []}
    for i in range(num_layers):
        layer = constructLayerWithCycles(genes, regulatory_dict)
        if layer:
            genes = layer['genes']
            safeLayer = layer.copy()
        else:
            break
    if not layer:
        if i == 0:
            safeLayer = {'sources': [], 'targets': [], 'interaction': [], 'genes': []}
        return safeLayer
    else:
        layer['layer_hasCycles'] = [i + 1, False]
    return layer


def checkRegulatoryNetwork(regulatoryDict):
    # Banned symbols
    bannedSymbols = ['*', '+', '/', ',', '~', "'", '-', ':', ';']
    # Implement a bidirectional dictionary as a more efficient way to store the reference
    referenceDict = bidict()
    # Check if the regulatory network has valid sources and targets, raise warning if not and suppress the character
    for i, zipped in enumerate(zip(regulatoryDict['source_ENSEMBL'], regulatoryDict['target_ENSEMBL'])):
        source, target = zipped
        if any(symbol in source for symbol in bannedSymbols) or source[0].isdigit() or 'and' in source or 'or' in source:
            #flag = f'{source} contains banned symbols, removing them, string has been modified'
            #warnings.warn(flag)
            # Encode the string to remove banned symbols
            if source in referenceDict.values():
                identifier = referenceDict.inverse[source]
            else:
                identifier = f'M_I_{str(uuid.uuid4())}'.replace('-', '')
            referenceDict[identifier] = source
            regulatoryDict['source_ENSEMBL'][i] = identifier
        if any(symbol in target for symbol in bannedSymbols) or target[0].isdigit() or 'and' in target or 'or' in target:
            #flag = f'{target} contains banned symbols, removing them, string has been modified'
            #warnings.warn(flag)
            # Encode the string to remove banned symbols
            if target in referenceDict.values():
                identifier = referenceDict.inverse[target]
            else:
                identifier = f'M_I_{str(uuid.uuid4())}'.replace('-', '')
               
            regulatoryDict['target_ENSEMBL'][i] = identifier
        
    return regulatoryDict, referenceDict

def transformReactionIdsIntoIndexes(gDict, reactionList):
    newGDict = defaultdict(list)
    for key, value in gDict.items():
        reactionIndexes = []
        for reaction in value:
            reactionIndexes.append(reactionList.index(reaction))
        newGDict[key] = reactionIndexes
    return newGDict

# Decode a string
def decode_string(string: str, referenceDict: dict):
    list_strings = string.split("_K")
    if len(list_strings) > 1:
        encoded_string = list_strings[0]
        last_name = list_strings[1]
        last_name = f"_K{last_name}"
    else:
        encoded_string = string
        last_name = ""
    decoded_string = referenceDict[encoded_string]
    decoded_string = decoded_string + last_name
    return decoded_string


def boolTransformer(expandedRules, value):
    newDict = {}
    # transform the grp into a boolean form
    trans_gpr = str(value['gpr'])
    trans_gpr = trans_gpr.replace(' and ', ' & ').replace(' or ', ' | ')
    newDict['gpr'] = trans_gpr

    for source, target, interaction in zip(expandedRules['sources'], expandedRules['targets'], expandedRules['interaction']):
        if interaction == -1:
            source = f'!{source}'
        if target not in newDict.keys():
            newDict[target] = source
        else:
            newDict[target] = f'{newDict[target]} | {source}'
    
    for gene in expandedRules['genes']:
        if gene not in newDict.keys():
            newDict[gene] = gene
            
    return newDict

def boolSaver(boolForm, path):
    f = mpbn.MPBooleanNetwork(boolForm)
    f.save(path)
    
    
    
def calculate_G_Matrix_GMIS(cobraModel, regulatory_dataframe, num_layers=2, **kwargs):
    '''
    '''
    # Name of the analysis
    name = (
        "gMCSpy_GMIS_"
        + cobraModel.id
    )
    
    # Prefered solver is gurobi due to performance experience during benchmarking
    solver = kwargs.get("solver", "gurobi")
    maxKOLength = kwargs.get('maxKOLength', 3)
    numWorkers = kwargs.get("numWorkers", 1)
    isoformSeparator = kwargs.get("isoformSeparator", None)
    saveGMatrix = kwargs.get("saveGMatrix", False)
    onlyGenes = kwargs.get('onlyGenes', False)
    geneSubset = kwargs.get('geneSubset', None)
    modifier = kwargs.get('modifier', None)

    # Save all the parameters in a log file    
    path = "logs/" + name + "/"
    createLog(path)
    configName = ("configuration")
    logging = initiateLogger(configName + ".log", path + configName + ".log")
    logging.info(f"Model: {cobraModel.id}")
    # close logger
    handler = logging.handlers[0]
    handler.close()
    logging.removeHandler(handler)
    
    path = "logs/" + name + "/process/"


    # Prepare the regulatory network
    regulatory_dict = regulatory_dataframe.to_dict('list')
    # Check if the regulatory network has valid sources and targets
    regulatory_dict, referenceDict = checkRegulatoryNetwork(regulatory_dict)
    # Check if the regulatory network has the correct format, should have source_ENSEMBL, target_ENSEMBL, and interaction
    expectedColumns = ['source_ENSEMBL', 'target_ENSEMBL', 'interaction']
    for column in expectedColumns:
        if column not in regulatory_dict:
            raise ValueError(f"Expected column {column} not found in the regulatory network")
    
    # Filter the regulatory network to keep only the relationships that are genes
    if onlyGenes:
        # Drop the rows that are not genes
        regulatory_dataframe = regulatory_dataframe[regulatory_dataframe['source_ENSEMBL'].str.contains('ENSG')]
        regulatory_dataframe = regulatory_dataframe[regulatory_dataframe['target_ENSEMBL'].str.contains('ENSG')]
    
    #print(regulatory_dataframe.shape)    
    

    startTime = time.time()
    gObject = calculateRegNetGMatrix(cobraModel, regulatory_dict, num_layers, solver=solver, maxKOLength=maxKOLength, path=path, onlyGenes=onlyGenes, numWorkers=numWorkers)
    endTime = time.time()
    logging.info("GMatrix: " + str(endTime - startTime))

    # Remove the folder bnAnalysis if it exists
    bnAnalysisPath = Path(path, 'bnAnalysis')
    if bnAnalysisPath.exists():
        shutil.rmtree(bnAnalysisPath)
        
    
    gDict = gObject["gDict"]

    # Filter the G matrix to keep only genes that are in the geneSubset
    if geneSubset:
        #Make a set of the geneSubset
        geneSet = set(geneSubset)
        # Extract all possible interventions
        gDictKeys = list(gDict.keys())
        # Remove the interventions that are not in the geneSubset
        for key in gDictKeys:
            if not key.issubset(geneSet):
                gDict.pop(key)
   
    if saveGMatrix:        
        gMatrix_Strings = []
        for ko, reactions  in gDict.items():
            string_ko = "@@".join(list(ko))
            string_reactions = "@@".join(cobraModel.reactions[x].id for x in reactions)
            row = string_ko + " -> " + string_reactions
            gMatrix_Strings.append(row)
        filename = f"Gmat_python_{cobraModel.id}_{solver}.csv"
        if modifier: 
            filename = filename.replace('.csv', f'_{modifier}.csv')
            
        gPath = path + "/GMatrix/" 
        createLog(gPath)
        completePath = gPath + filename
        df = pd.DataFrame({'gMatrix': gMatrix_Strings})
        df.to_csv(completePath)

    return gDict
from .ProblemInterpreters import np
from .ProblemInterpreters import beartype

from .ProblemDefinitions import cobra
from .ProblemDefinitions import List
from .ProblemDefinitions import scipy

from .calculateMCS import calculateMCS

import re
from ast import And, BoolOp, Or, Name, parse

from collections import defaultdict
import warnings

from .MinimalCutSetClass import uuid
from .MinimalCutSetClass import tqdm
from .MinimalCutSetClass import logging



def buildGMatrix(
    cobraModel: cobra.core.Model,
    maxKOLength: int = 1e6,
    name: str = None,
    isoformSeparator: str = None,
    verbose: int = 0,
    **kwargs,
):
    """

    **Build the G matrix, G dict and the relationships between perturbations**

    The G matrix is a matrix with the perturbations of the genes in the model, that render a reaction infeasible.
    Each row of the matrix corresponds to a perturbation, and each column corresponds to a reaction.
    This matrix is used to calculate the genetic minimal cut sets of the model.

    :param cobraModel: The cobra model that represents the metabolism of an organism.
    :param name: The name of the model. If None, the id of the model will be used.
    :param separateIsoform: If not None, the isoforms of the genes will merge into one gene.
    :param verbose: The level of verbosity of the function.

    :return: A List with the G matrix, the G dictionary and the relationships between perturbations.


    """
    
    maxKOLength = kwargs.get('maxKOLength', 3)
    dictManager = mergeIsforms(isoformSeparator)
        
    modelReactions = [reaction.id for reaction in cobraModel.reactions]
    # create a dictionary with the gpr rules associated to each reaction to get a unique list of gprs
    gprDict = getGPRDict(cobraModel)
    gDict = {}
    
    # Analyze the GPRs to find the perturbations that stop the reactions
    for key, value in tqdm(gprDict.items()):
            calculate_GPR_MCS(gDict, value['gpr'], value['reactions'], dictManager)


    # Filter out the interventions that are larger than the maxKOLength
    filteredGDict = gDict.copy()
    for key in gDict.keys():
        if len(key) > maxKOLength:
            filteredGDict.pop(key)
            
    newGDict = transformReactionIdsIntoIndexes(filteredGDict, modelReactions)
    
            
    simpG = simplifyGMatrix(newGDict, dictManager)

    [relationships, numberNewGenesByKO] = relatedRows(simpG, dictManager)
    
    gMatrix = createSparseMatrix(simpG, modelReactions)

    # scipy.sparse.save_npz(cobraModelName + '_gMatrix.npz', gMatrix)
    gObject = {
        "gMatrix": gMatrix,
        "gDict": simpG,
        "relationships": relationships,
        "numberNewGenesByKO": numberNewGenesByKO,
    }
    return gObject

def calculate_GPR_MCS(gDict, gpr, reactions, addToDict):
    # get the genes from the gpr rule
    stringGPR = gpr.to_string()
    genes = [gene for gene in gpr.genes]
    if len(genes) == 1:
        key=frozenset(genes)
        [addToDict(gDict, key, reaction) for reaction in reactions]
    # Only and
    elif bool(re.search(' and ', stringGPR) and not re.search(' or ', stringGPR)):
        # Only one gene needs to be KO to stop the reaction
        for gene in genes:
            key=frozenset([gene])
            [addToDict(gDict, key, reaction) for reaction in reactions]

    # Only or
    elif bool(re.search(' or ', stringGPR) and not re.search(' and ', stringGPR)):
        # All genes need to be KO to stop the reaction
        key=frozenset(genes)
        [addToDict(gDict, key, reaction) for reaction in reactions]  
    
    # And and or
    else:
        treeGenes, model = parseGPRToModel(stringGPR, 'gpr', None)
        solutionDict = calculateMCS(model, MAX_LENGTH_MCS=5, MAX_MCS=1e6, rxnSubset=list(treeGenes))
        for _, data in solutionDict.items():
            key = frozenset(data["solution"])
            [addToDict(gDict, key, reaction) for reaction in reactions]
                   
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
    return gpr_dict

def simplifyGMatrix(gMatrixDict: defaultdict, __addToDict) -> defaultdict:
    """**Simplify the G matrix**

    Analyze the G matrix to find the rows that are related to each other.
    If a common gene is more perturbations, then the perturbations are related to each other.
    And should be an individual perturbation.

    :param gMatrixDict: The dictionary with the perturbations.

    :return: The simplified G matrix.
    """
    
    # Empty dictionary to store the simplified G matrix
    simpG = gMatrixDict.copy()
    
    # Get the interventions compossed of one gene
    singleGeneInterventions = [gene for gene in gMatrixDict.keys() if len(gene) == 1]
    
    # Get the interventions compossed of more than one gene
    multipleGeneInterventions = [gene for gene in gMatrixDict.keys() if len(gene) > 1]
    
    # Intersect each intervention with the other interventions to find common elements
    for interventionA in multipleGeneInterventions:
        for interventionB in multipleGeneInterventions:
            if interventionA != interventionB:
                intersection = interventionA.intersection(interventionB)
                if intersection:
                    # If the intersection is not empty, then the interventions are related
                    # Check if the intersection is already a single gene intervention
                    if len(intersection) == 1:
                       if intersection not in singleGeneInterventions and intersection not in gMatrixDict.keys():
                           # Intersection must not be a key of the dictionary
                            if intersection not in simpG.keys():
                                 # Add the intersection as a key of the dictionary
                                __addToDict(simpG, intersection, [])
                                simpG[intersection] = []
                    else:
                        for element in intersection:
                            if element not in singleGeneInterventions and element not in gMatrixDict.keys():
                                # Intersection must not be a key of the dictionary
                                if intersection not in simpG.keys():
                                     # Add the intersection as a key of the dictionary
                                    __addToDict(simpG, intersection, [])
                                    simpG[intersection] = []
                                

    return simpG

def createSparseMatrix(gMatrixDict, modelReactions):
    """**Build the G matrix**

    Create sparse matrix from the dictionary with the perturbations.

    :param gMatrixDict: The dictionary with the perturbations.
    :param modelReactions: The reactions of the model.

    :return: The G matrix.
    """
    gMatrix = scipy.sparse.lil_matrix((len(gMatrixDict), len(modelReactions)))
    for i, items in enumerate(gMatrixDict.items()):
        gMatrix[i, list(items[1])] = 1
    gMatrix = gMatrix.tocsc()
    return gMatrix

def relatedRows(gMatrixDict: defaultdict, __addToDict) -> List:
    """**Find the rows that are related to each other**

    If a row is a subset of another row,
    then that means that when the perturbation of the superset is activated, the subset is also affected.
    i. e. gene A and gene B inactivates reaction1, and gene A inactivates reaction2. Knocking out gene A and gene B will affect both reactions.

    :param gMatrixDict: The dictionary with the perturbations.

    :return: A list with the rows that are related to each other.
    """
    relationships = defaultdict(list)
    numberNewGenesByKO = defaultdict(int)
    sortedDict = sorted(gMatrixDict.items(), key=lambda x: len(x[0]))
    for key, value in gMatrixDict.items():
        numberNewGenesByKO[key] = len(key)
        if len(key) < 2:
            continue
        keysConstained = []
        for key2, value2 in sortedDict:
            if len(key2) >= len(key):
                break
            if key2.issubset(key):
                keysConstained.append(key2)
                value.extend(value2)
                __addToDict(relationships, key, key2)
        if len(keysConstained) > 0:
            key2 = frozenset.union(*keysConstained)
            numberNewGenesByKO[key] = len(key) - len(key2)
        gMatrixDict[key] = list(set(value))
    return [relationships, numberNewGenesByKO]

def calculateMCSForGPRs(
    briefing: dict, isoformSeparator: str, verbose: int = 0, **kwargs
):
    """**Calculate the Minimal Cut Sets for GPRs with both AND and OR relationships**

    Given a complex GPR, we parse de GPR to a cobra model and calculate the MCSs for the model.

    :param briefing: The dictionary with the information of the GPRs.
    :param verbose: The verbosity level.

    :return: A List with the List of GPRs and a List with the MCSs.
    """
    gpr = []
    mcs = []
    parmMCS = {
        "MAX_MCS": 10000,
        "targetB": 1e-3,
        "timelimit": 60,
        "forceLength": True,
        "numWorkers": 0,
        "verbose": verbose,
        "useIndicators": True,
        "solver": "gurobi",
    }
    for key, value in kwargs.items():
        if key in parmMCS:
            parmMCS[key] = value
        else:
            warnings.warn("Parameter " + key + " is not defined in calculateMCS.")

    for i, stringGPR in enumerate(briefing["gprsANDOR"]):
        treeGenes, model = parseGPRToModel(stringGPR, i, isoformSeparator)
        if verbose > 0:
            for i in model.reactions:
                print(i)

        solutionDict = calculateMCS(
            model, MAX_LENGTH_MCS=len(treeGenes), rxnSubset=list(treeGenes), **parmMCS
        )
        for key, data in solutionDict.items():
            gpr.append(stringGPR)
            mcs.append(data["solution"])
    return gpr, mcs

@beartype
def parseGPRToModel(stringGPR: str, reactionName, isoformSeparator):
    """**Parse a GPR string to a cobra model**


    The model is used to calculate the MCS to calculate the combinations of genes that are required for the reaction.

    :param stringGPR: The GPR string to parse. (e.g. "gen1 and gen2 or gen3")
    :param reactionName: The name of the reaction to parse. Just to identify the model.

    :return:
        - **treeGenes** - A list of genes that participate in the GPR.
        - **model** - A cobra model with the GPR parsed to a metabolic model.
    """
    if isoformSeparator is not None:
        # Remove the isoform separator and the following string and numbers
        stringGPR = re.sub(f"([{isoformSeparator}][0-9]+)", "", stringGPR)
    [gprTree, treeGenes] = cobra.core.gene.parse_gpr(stringGPR)
    parentMetabolite = cobra.core.Metabolite(f"r{0}", name=f"r{0}")
    model = cobra.core.Model(f"reaction_{reactionName}")
    objectiveReaction = cobra.core.Reaction("objective")
    objectiveReaction.name = "objective"
    objectiveReaction.lower_bound = 0
    objectiveReaction.upper_bound = 1000
    objectiveReaction.add_metabolites({parentMetabolite: -1.0})
    model.add_reactions([objectiveReaction])
    model.objective = "objective"
    level = 0
    model = buildReactionFromBranch(expr=gprTree.body, model=model, level=level)

    return [treeGenes, model]

@beartype
def buildReactionFromBranch(
    expr: BoolOp,
    model: cobra.core.Model,
    level: int = 0,
    parentMetabolite: cobra.core.Metabolite = None,
) -> cobra.core.Model:
    """**Recursive function to generate metabolic model from GPR.**

    Each GPR can be expressed as a metabolic model. To that end, a boolean tree can be build, where a branch can be either an AND, an OR or a gene.
    We move through the tree, and recursively call this function on each node. Until all genes are explored (end nodes).
    As we move through the tree, we add reactions to the model. The reactions are named with a unique identifier, to avoid name conflicts.
    An example building a metabolic model from a GPR is shown below.

    .. figure:: ./figs/binaryTreeExample.png
        :height: 350px
        :width: 450 px
        :align: center

        Figure 1. Example of a binary tree.

    In the example above, the GPR is (g2 and (g5 or g6)) or g7. The function will be called recursively on each node.

    :param expr: A BoolOp object from the ast module to be explored. i.e. a tree.  (ast.BoolOp)
    :param model: A cobra model to build add the tree.
    :param level: The level of the tree, i.e. the depth of the current node.  (int)
    :param parentMetabolite: The metabolite that is the parent of the current node.  (cobra.core.Metabolite)

    :returns:
            - **model**: A cobra model with the tree added as reactions, genes have become metabolites, and the objective function set to the root node, the reaction of interest.
    """
    if isinstance(expr, BoolOp):
        branch = expr.op
        if isinstance(branch, Or):
            if parentMetabolite is None:
                parentMetabolite = cobra.core.Metabolite(
                    f"r{level}", name=f"reaction_{level}"
                )
            reaction = cobra.core.Reaction(str(uuid.uuid4()))
            reaction.name = f"OR_{level}_{parentMetabolite.id}"
            reaction.lower_bound = 0
            reaction.upper_bound = 1000
            reaction.add_metabolites({parentMetabolite: 1.0})

            for i in expr.values:
                if isinstance(i, Name):
                    level = level + 1
                    reaction = cobra.core.Reaction(str(uuid.uuid4()))
                    reaction.name = f"OR_{level}_{parentMetabolite.id}"
                    reaction.lower_bound = 0
                    reaction.upper_bound = 1000
                    reaction.add_metabolites({parentMetabolite: 1.0})
                    childMetabolite = cobra.core.Metabolite(i.id, name=i.id)
                    reaction.add_metabolites({childMetabolite: -1.0})
                    model.add_reactions([reaction])
                    if i.id in model.reactions:
                        continue
                    reaction = cobra.core.Reaction(i.id)
                    reaction.name = i.id
                    reaction.lower_bound = 0
                    reaction.upper_bound = 1000
                    reaction.add_metabolites({childMetabolite: 1.0})
                    model.add_reactions([reaction])
                else:
                    level = level + 1
                    reaction = cobra.core.Reaction(str(uuid.uuid4()))
                    reaction.name = f"OR_{level}_{parentMetabolite.id}"
                    reaction.lower_bound = 0
                    reaction.upper_bound = 1000
                    reaction.add_metabolites({parentMetabolite: 1.0})
                    childMetabolite = cobra.core.Metabolite(
                        f"branch_OR_met{level}_{parentMetabolite.id}",
                        name=f"reaction_{level}",
                    )
                    reaction.add_metabolites({childMetabolite: -1.0})
                    model.add_reactions([reaction])
                    buildReactionFromBranch(i, model, level, childMetabolite)
        elif isinstance(branch, And):
            if parentMetabolite is None:
                parentMetabolite = cobra.core.Metabolite(
                    f"r{level}", name=f"reaction_{level}"
                )
            reaction = cobra.core.Reaction(str(uuid.uuid4()))
            reaction.name = f"AND_{level}_{parentMetabolite.id}"
            reaction.lower_bound = 0
            reaction.upper_bound = 1000
            reaction.add_metabolites({parentMetabolite: 1.0})
            for i in expr.values:
                if isinstance(i, Name):
                    childMetabolite = cobra.core.Metabolite(i.id, name=i.id)
                    reaction.add_metabolites({childMetabolite: -1.0})
                    if i.id in model.reactions:
                        continue
                    newReaction = cobra.core.Reaction(i.id)
                    newReaction.name = i.id
                    newReaction.lower_bound = 0
                    newReaction.upper_bound = 1000
                    newReaction.add_metabolites({childMetabolite: 1.0})
                    model.add_reactions([newReaction])
                else:
                    level = level + 1
                    childMetabolite = cobra.core.Metabolite(
                        f"branch_AND_met{level}_{parentMetabolite.id}",
                        name=f"reaction_{level}",
                    )
                    reaction.add_metabolites({childMetabolite: -1.0})
                    buildReactionFromBranch(i, model, level, childMetabolite)
            model.add_reactions([reaction])

    return model

def mergeIsforms(isoformSeparator):
    if isoformSeparator is None:

        def __addToDict(dict: defaultdict, key: List, value) -> None:
            """**Add to dictionary**

            Internal function to add a value to a dictionary with a key. If the key does not exist, it is created,
            and the value is added to the list of values.
            If the key exists, the value is appended to the list of values.

            Dictionary mut be a `defaultdict <https://docs.python.org/3/library/collections.html#collections.defaultdict>`_.

            :param dict: Dictionary to add gene reaction interaction. (defaultdict)
            :param key: List of genes that render the reaction inactivate. (list), if is only one gene, it must be a list with one element.
            :param value: Any value can be stored in the dictionary, we store the corresponding index of the reactions perturbed. (any)

            """
            #key = frozenset(key)
            if key in dict:
                dict[key].append(value)
            else:
                dict[key] = [value]

    elif isinstance(isoformSeparator, str):

        def __addToDict(dict: defaultdict, key: List, value) -> None:
            key = [name.split(isoformSeparator)[0] for name in key]
            key = frozenset(key)
            if key in dict:
               dict[key].append(value)
            else:
                dict[key] = [value]

    else:
        raise TypeError("isoformSeparator must be a string or None")

    return __addToDict

@beartype
def buildRxnGeneMat(GPRs: List, reactions: List, genes: List):
    """**Build a reaction-gene matrix from a list of GPRs.**

    The binary matrix is build from each reaction GPR,
    each row corresponds to a reaction and each column corresponds to a gene.
    The matrix is sparse and is stored in the CSR format.

    :param GPRs: List of GPRs from a cobra model.
    :param reactions: List of reactions from a cobra model.
    :param genes: List of genes from a cobra model.

    :returns:
        - **rxnGeneMatrix:** A sparse matrix with the number of reactions in rows and the number of genes in columns. (scipy.sparse.csr_matrix)
        - **summary:** A dictionary with the number of reactions classified by the number of genes in their GPR. (dict)
    """
    rxnGeneMatrix = scipy.sparse.lil_matrix((len(reactions), len(genes)))

    # Build reaction-gene matrix
    for position, reaction in enumerate(reactions):
        GPR = reaction.gpr.genes
        for gene in GPR:
            rxnGeneMatrix[position, genes.index(gene)] = 1

    rxnGeneMatrix = rxnGeneMatrix.tocsr()

    rxnZeroGenes = np.where(rxnGeneMatrix.sum(axis=1) == 0)[0]
    rxnOneGene = np.where(rxnGeneMatrix.sum(axis=1) == 1)[0]

    # indices of reactions not in rxnZeroGenes and rxnOneGene
    indexRxnMoreThanOneGene = np.setdiff1d(
        np.arange(len(reactions)), np.concatenate((rxnZeroGenes, rxnOneGene))
    )

    # Keep only GPRs that include two or more genes
    GPRStrings = []
    for gpr in GPRs:
        if len(gpr.genes) >= 2:
            GPRStrings.append(gpr.to_string())

    # identify reactions with "or" and "and" relationships between genes
    rxnOr = [bool(re.search(" or ", i)) for i in GPRStrings]
    rxnAnd = [bool(re.search(" and ", i)) for i in GPRStrings]

    rxnOnlyOr = np.logical_and(rxnOr, np.logical_not(rxnAnd))
    rxnOnlyAnd = np.logical_and(rxnAnd, np.logical_not(rxnOr))
    rxnOrAnd = np.logical_and(rxnOr, rxnAnd)

    dictGPRstring = {
        position: reactions[position].gpr.to_string()
        for position in indexRxnMoreThanOneGene
    }

    gprOnlyOr = np.where(rxnOnlyOr)[0]
    indexRxnOnlyOr = []
    for gpr in gprOnlyOr:
        uniqueGPR = GPRStrings[gpr]
        indexRxnOnlyOr.extend(
            [key for key, value in dictGPRstring.items() if value == uniqueGPR]
        )

    gprRxnOnlyAnd = np.where(rxnOnlyAnd)[0]
    indexRxnOnlyAnd = []
    for gpr in gprRxnOnlyAnd:
        uniqueGPR = GPRStrings[gpr]
        indexRxnOnlyAnd.extend(
            [key for key, value in dictGPRstring.items() if value == uniqueGPR]
        )

    gprRxnOrAnd = np.where(rxnOrAnd)[0]
    indexRxnOrAnd = []
    for gpr in gprRxnOrAnd:
        uniqueGPR = GPRStrings[gpr]
        indexRxnOrAnd.extend(
            [key for key, value in dictGPRstring.items() if value == uniqueGPR]
        )

    # count number of reactions in each category
    nRxnZeroGenes = len(rxnZeroGenes)
    nRxnOneGene = len(rxnOneGene)
    nRxnOnlyOr = len(indexRxnOnlyOr)
    nRxnOnlyAnd = len(indexRxnOnlyAnd)
    nRxnOrAnd = len(indexRxnOrAnd)

    nRxnTotal = nRxnZeroGenes + nRxnOneGene + nRxnOnlyOr + nRxnOnlyAnd + nRxnOrAnd

    # Keep only GPRs that will be used in the analysis of minimal cut sets
    UniqueGPRStringsANDOR = [GPRStrings[i] for i in gprRxnOrAnd]

    # summarize results as a dictionary
    summary = {
        "nRxnZeroGenes": nRxnZeroGenes,
        "indexRxnZeroGenes": rxnZeroGenes,
        "nRxnOneGene": nRxnOneGene,
        "indexRxnOneGene": rxnOneGene,
        "nRxnOnlyOr": nRxnOnlyOr,
        "indexRxnOnlyOr": indexRxnOnlyOr,
        "nRxnOnlyAnd": nRxnOnlyAnd,
        "indexRxnOnlyAnd": indexRxnOnlyAnd,
        "nRxnOrAnd": nRxnOrAnd,
        "gprsANDOR": UniqueGPRStringsANDOR,
        "gprDict": dictGPRstring,
        "nRxnTotal": nRxnTotal,
    }

    return [rxnGeneMatrix, summary]


def initiateLogger(name, log_file, level=logging.INFO):
    handler = logging.FileHandler(log_file)
    logger = logging.getLogger(name)
    logger.setLevel(level)
    logger.addHandler(handler)

    return logger

def transformReactionIdsIntoIndexes(gDict, reactionList):
    newGDict = defaultdict(list)
    for key, value in gDict.items():
        reactionIndexes = []
        for reaction in value:
            reactionIndexes.append(reactionList.index(reaction))
        newGDict[key] = reactionIndexes
    return newGDict



        
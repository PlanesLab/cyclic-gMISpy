from cobra import Model, Reaction, Metabolite
from beartype import beartype
import os
from datetime import datetime
import pandas as pd
import re


@beartype
def saveSolutions(
    problemDict: dict, solutionDict: dict, modelStruct: Model, filename: str, path: str
):
    """
    Saves the solutions to a csv file.

    .. warning::
        This function has been design to save the solution of our optimization problems, it may not work for other problems.

    :param problemDict: Dictionary containing the problem information.
    :param solutionDict: Dictionary containing the solution information.
    :param modelStruct: The cobra Model used to build the problem dictionary.
    :param filename: Name of the file to save the solutions.
    :param path: Path to save the file.


    """
    names = []
    orders = []
    reactionSets = []
    for key, data in solutionDict.items():
        vars = []
        for i in data:
            solPos = problemDict["variables"]["variableNames"].index(str(i))
            corresPondingVar = modelStruct.reactions[
                problemDict["variables"]["variableGroups"]["z"].index(solPos)
            ].id
            vars.append(corresPondingVar)
        orderNumber = key.split("order")[1].split("_")[0]
        names.append(key)
        orders.append(orderNumber)
        reactionSets.append(vars)

    df = pd.DataFrame({"solName": names, "order": orders, "Reactions": reactionSets})
    df.to_csv(path + "/" + filename + ".csv", index=False)

@beartype
def saveSolutionDict(dict: dict, filename: str, path: str):
    dataframe = pd.DataFrame.from_dict(dict, orient="columns")
    dataframe.to_csv(path + "/" + filename + ".csv", index=False)

def setSolver(defaultSolver, useIdicators: bool = True):
    """Set the solver to be used for the optimization problem.
    This function is used automatically set the correct interface based on the solver chosen.
    Current options are 'cplex', 'gurobi' and 'scip'.

    :param defaultSolver: The solver to be used for the optimization problem. Must be one of the options above.
    :param indicator: If True, the solver will be set to use indicators. If False, the solver will be set to not use indicators.
        The function also checks if the solver has the capability to use indicators. If the solver does not have the capability
        it returns turns the indicators off and the Big M method is implemented.

    :return:
        - The interface to be used for the optimization problem.
        - A boolean indicating if the solver will use indicators or not.
    """
    if defaultSolver == "cplex":
        try:
            import cplex

            return [cplex, True]
        except ImportError :
            print(
                "cplex could not be loaded, check if it is installed and the path is correct"
            )
    elif defaultSolver == "gurobi":
        try:
            import gurobipy

            return [gurobipy, True]
        except ImportError :
            print(
                "gurobi could not be loaded, check if it is installed and the path is correct"
            )

    elif defaultSolver == "scip":
        try: 
            import pyscipopt

            return [pyscipopt, True]
        except ImportError :
            print(
                "pyscipopt could not be loaded, check if it is installed and the path is correct"
            )
    else:
        print(
            "The default solver is not supported, please select (cplex, gurobi, scip)"
        )

def createLog(folderLogPath: str = "logs/log"):
    """Creates a log file with the current date and time, and returns the filename.
    This functions also creates the folder to save the log file if it does not exist.
    It is used to save the log of the optimization process depending on the verbosity level."""
    filename = folderLogPath
    os.makedirs(os.path.dirname(filename), exist_ok=True)
    print(filename)
    return filename

def loadSolutions(path: str = None):
    """Load the solutions generated from the calculation of gMCSs.

    :param path: Path to the folder with the solutions, if no path is given the most recent calculation in the logs folder is used.

    :return: A dictionary with the solutions.
    """
    if path is None:
        path = 'logs'
        # Get a list of all directories in the given path
        directories = [d for d in os.listdir(path) if os.path.isdir(os.path.join(path, d))]
        
        # Sort the directories by their modification time in descending order
        directories.sort(key=lambda x: os.path.getmtime(os.path.join(path, x)), reverse=True)
        
        if directories:
            # The most recent folder is the first one in the sorted list
            most_recent = directories[0]
            complePath = path + "/" + most_recent + "/solutions/solutions.log"
            return read_sol_log(complePath)
            
        else:
            raise Warning("No solutions found in the given path.")
    else:
        return read_sol_log(path)  

def read_sol_log(logFile):
    # Read the log file
    with open(logFile, "r") as file:
        # First line are headers
        lines = file.readlines()[1:]
        solDict = {}
        for line in lines:
            line = line.strip()
            # Split the line by the separator
            lineList = line.split(":")
            # Get the genes
            geneList = add_quotes_inside_brackets(lineList[0])
            geneList = eval(geneList)
            geneList = frozenset(geneList)
            # Get the reactions in the solution
            reactionList = eval(add_quotes_inside_brackets(lineList[1]))
            solDict[geneList] = reactionList
        
        return solDict
        
def add_quotes_inside_brackets(input_string):
    # Use regular expressions to find text inside square brackets
    pattern = r'\[([^]]*)\]'
    
    # Define a function to add quotes around the found text
    def add_quotes(match):
        # Split the text inside brackets by ',' and add quotes to each element
        elements = match.group(1).split(',')
        quoted_elements = [f'"{element.strip()}"' for element in elements]
        return '[' + ','.join(quoted_elements) + ']'
    
    # Use re.sub() to apply the function to the input string
    result = re.sub(pattern, add_quotes, input_string)
    
    return result

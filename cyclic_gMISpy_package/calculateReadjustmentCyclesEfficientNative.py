from .calculateGMCSwithCycles import getGPRDict, stackLayersWithCycles, checkRegulatoryNetwork, bidict
from .MinimalCutSetClass import tqdm
from .cobraTools import re, initiateLogger
from .ProblemDefinitions import readjustmentProblem
from copy import deepcopy
from datetime import datetime
from .Utilities import createLog
from .calculateGMCSwithCycles import boolTransformer, mpbn
from concurrent.futures import ProcessPoolExecutor, as_completed
import multiprocessing as mp
from functools import partial
import logging as std_logging
from collections import defaultdict


def calculateReadjustmentParallel(gMISList, model, regulatory_dataframe, num_layers=2, solver='gurobi', **kwargs):
    """
    Optimized version that processes gene-gMIS pairs and stops at first readjustment per pair.
    """
    name = (
        "gMCSpy_Readjustment_"
        + model.id
        + "_"
        + datetime.now().strftime("%H-%M-%S--%d-%m-%Y")
    )
    
    pathN = "logs/" + name + "/"
    path = kwargs.get('path', pathN)
    
    createLog(path)
    configName = ("configuration")
    logging = initiateLogger(configName + ".log", path + configName + ".log")
    logging.info(f"solver: {solver}, num_layers: {num_layers}")
    logging.info(f"Model: {model.id}")
    logging.info(f"gMIS to check: {gMISList}")
    handler = logging.handlers[0]
    handler.close()
    logging.removeHandler(handler)
    
    # Results logger setup
    if kwargs.get('path', None) is None:
        results_path = "logs/" + name + "/results/"
        results_name = "results"
        createLog(results_path + results_name + ".log")
        results_logging = initiateLogger(results_name + ".log", results_path + results_name + ".log")
    else:
        results_name = "results.log"
        createLog(path + '/' + results_name)
        results_logging = initiateLogger(results_name, path + '/' + results_name)
        
    kwargs.pop('path', None)
    
    # Prepare the regulatory network
    regulatory_dict = regulatory_dataframe.to_dict('list')
    regulatory_dict, referenceDict = checkRegulatoryNetwork(regulatory_dict)
    gprDict = getGPRDict(model)
    
    # Expand all the rules to the desired number of layers
    rulesDict = {}
    i = 0
    print("Calculating rules...")
    for key, value in tqdm(gprDict.items()):
        GPR = key, value
        rules, genes = calculateRules(GPR, regulatory_dict, num_layers, solver)
        rulesDict[i] = {'genes': genes, 'rules': rules}
        i += 1
    
    # Build the optimized data structure: gene-gMIS pairs with their associated rules
    print("Building gene-gMIS-rule data structure...")
    gene_gMIS_rule_pairs = build_gene_gMIS_rule_structure(gMISList, rulesDict, referenceDict)
    
    print(f'Processing {len(gene_gMIS_rule_pairs)} gene-gMIS-rule pairs...')
    
    # Get number of workers
    numWorkers = kwargs.get('numWorkers', mp.cpu_count())
    print(f'Number of workers: {numWorkers}')
    
    # Process gene-gMIS pairs in parallel with early stopping
    with ProcessPoolExecutor(max_workers=numWorkers) as executor:
        # Create partial function with common parameters
        process_func = partial(process_gene_gMIS_pair_wrapper)
        
        # Submit all tasks
        future_to_pair = {
            executor.submit(process_func, pair_data): pair_data['pair_id']
            for pair_data in gene_gMIS_rule_pairs
        }
        
        # Process results as they complete
        all_results = []
        for future in tqdm(as_completed(future_to_pair), total=len(gene_gMIS_rule_pairs), desc="Processing gene-gMIS pairs"):
            pair_id = future_to_pair[future]
            try:
                result = future.result()
                if result:
                    all_results.append(result)
                    # Log result immediately
                    if result.get('readjustment', False):
                        results_logging.info(result)
            except Exception as exc:
                print(f'Gene-gMIS pair {pair_id} generated an exception: {exc}')
    
    # Close all handlers for logging_results
    for handler in results_logging.handlers[:]:
        handler.close()
        results_logging.removeHandler(handler)
    
    readjustment_count = len([r for r in all_results if r.get('readjustment', False)])
    print(f"Processing complete. Found {readjustment_count} readjustments.")
    return all_results


def build_gene_gMIS_rule_structure(gMISList, rulesDict, referenceDict):
    """
    Build a data structure with gene-gMIS pairs and their associated rules.
    This allows for efficient processing where each pair is processed independently.
    """
    gene_gMIS_rule_pairs = []
    pair_id = 0
    
    # Filter to valid gMIS (more than one gene)
    valid_gMIS = [gMIS for gMIS in gMISList if len(gMIS) > 1]
    
    for gMIS in valid_gMIS:
        # Encode the entire gMIS
        gMISDict = {}
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
        
        # For each gene in the gMIS, create gene-gMIS pairs with relevant rules
        for target_gene in gMIS:
            encoded_target = gMISDict[target_gene]
            
            # Find all rules that contain this target gene
            relevant_rules = []
            for rule_key, rule_value in rulesDict.items():
                if encoded_target in rule_value['genes']:
                    relevant_rules.append({
                        'rule_key': rule_key,
                        'rules': rule_value['rules']
                    })
            
            # Only create pair if there are relevant rules
            if relevant_rules:
                gene_gMIS_rule_pairs.append({
                    'pair_id': pair_id,
                    'target_gene': target_gene,
                    'gMIS': gMIS,
                    'gMISDict': gMISDict,
                    'relevant_rules': relevant_rules
                })
                pair_id += 1
    
    return gene_gMIS_rule_pairs


def process_gene_gMIS_pair_wrapper(pair_data):
    """
    Process a single gene-gMIS pair with early stopping.
    Returns as soon as readjustment is found for this pair.
    """
    try:
        target_gene = pair_data['target_gene']
        gMIS = pair_data['gMIS']
        gMISDict = pair_data['gMISDict']
        relevant_rules = pair_data['relevant_rules']
        pair_id = pair_data['pair_id']
        
        # Check if the gene is knocked out or knocked in
        interaction = 0 if target_gene.endswith('_KO') else 1
        
        # Get the encoded target gene
        encoded_target = gMISDict[target_gene]
        
        # Get the accompanying genes (all genes except the target)
        accompanying_genes = [gene for gene in gMIS if gene != target_gene]
        
        # Early stopping: process rules until we find readjustment
        for rule_info in relevant_rules:
            rule_key = rule_info['rule_key']
            rules = rule_info['rules']
            
            # Skip if target gene not in this rule
            if encoded_target not in rules:
                continue
            
            # Calculate the value of the accompanying genes
            val = calculateSetValues(accompanying_genes, rules)
            if val is None:
                continue  # Only the target gene is relevant
            
            try:
                # Perform the knockout/knockin to the target gene
                network, reachable_from = boolean_Knockout(encoded_target, interaction, rules)
                
                # Prepare constraints for accompanying genes
                constraints = {}
                for gene in accompanying_genes:
                    if gene.endswith('_KI'):
                        gene_without_lastname = gene.split('_KI')[0]
                        tInt = 1
                    else:
                        gene_without_lastname = gene.split('_KO')[0]
                        tInt = 0
                    constraints[gene_without_lastname] = tInt
                
                # Try to find attractors
                act = network.attractors(reachable_from=reachable_from, constraints=constraints, limit=1)
                attractors = list(act)[0]
                
                # If we reach here, no readjustment for this rule, continue to next rule
                
            except Exception:
                # Exception indicates readjustment - return immediately (early stopping)
                return {
                    'readjustment': True,
                    'gene': target_gene,
                    'gMISDict': gMIS,
                    'rule_key': rule_key,
                    'pair_id': pair_id
                }
        
        # If we processed all rules without finding readjustment
        return {
            'readjustment': False,
            'gene': target_gene,
            'gMISDict': gMIS,
            'pair_id': pair_id
        }
        
    except Exception as e:
        print(f"Error processing gene-gMIS pair {pair_data.get('pair_id', 'unknown')}: {e}")
        return None


def calculateRules(GPR, regulatory_dict, num_layers=2, solver='gurobi'):
    """Calculate expanded rules for a given GPR."""
    key, value = GPR
    expandedRules = stackLayersWithCycles(value['genes'], regulatory_dict, num_layers, solver)
    boolform = boolTransformer(expandedRules, value)
    genes = list(boolform.keys())
    return boolform, genes


def calculateSetValues(gene_list, rules):
    """Calculate the sum of intervention values for genes in the list."""
    val = None
    for gene in gene_list:
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


def boolean_Knockout(target: str, interventionVal: int, rules: dict):
    """
    Create a new set of rules with the target gene knocked out.
    """
    newRules = deepcopy(rules)
    newRules.pop(target, None)  # Use pop with default to avoid KeyError
    
    # Creating a boolean network
    network = mpbn.MPBooleanNetwork(newRules)
    # Setting the target intervention value
    reachable_from = {target: interventionVal}
    
    return network, reachable_from


# Additional utility function for analysis
def analyze_gene_gMIS_structure(gene_gMIS_rule_pairs):
    """
    Analyze the gene-gMIS-rule structure for optimization insights.
    """
    print(f"Total gene-gMIS pairs: {len(gene_gMIS_rule_pairs)}")
    
    rules_per_pair = [len(pair['relevant_rules']) for pair in gene_gMIS_rule_pairs]
    print(f"Average rules per pair: {sum(rules_per_pair) / len(rules_per_pair):.2f}")
    print(f"Max rules per pair: {max(rules_per_pair)}")
    print(f"Min rules per pair: {min(rules_per_pair)}")
    
    # Count unique gMIS
    unique_gMIS = set(tuple(sorted(pair['gMIS'])) for pair in gene_gMIS_rule_pairs)
    print(f"Unique gMIS: {len(unique_gMIS)}")
    
    # Count unique genes
    unique_genes = set(pair['target_gene'] for pair in gene_gMIS_rule_pairs)
    print(f"Unique target genes: {len(unique_genes)}")
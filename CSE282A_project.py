import pandas as pd
import numpy as np
import copy
import random

#Converting the data into panda dataframes and visualising
metadata = pd.read_csv("metadata.csv")

filter1 = metadata["group"]=="NORMAL"
filter2 = metadata["group"]=="DEVIANT"

normal_sample_ids = metadata[filter1]["sample-id"].to_numpy()
normal_sample_ids = np.insert(normal_sample_ids, 0, 'taxonomy')

taxonomy = pd.read_csv("taxonomy_400.csv")
taxonomy_normal = taxonomy[normal_sample_ids]
taxonomy_normal.set_index("taxonomy")

# Saving entire list of ASVs or taxons being considered
taxonomy_range = taxonomy_normal['taxonomy']
taxonomy_range = list(taxonomy_range)

# Estimating population averaged relative abundances (mean, standard deviation and defined parameters) 
# of all taxons across the entire set of "NORMAL" samples

reference_data = []

for i in range(len(taxonomy_range)):
    
        ref_taxon = taxonomy_range[i]
        ref_data = taxonomy_normal.loc[(taxonomy_normal['taxonomy'] == ref_taxon)].to_numpy()
        ref_data = np.delete(ref_data, 0)
        ref_data = ref_data.astype(str)
        ref_data = np.char.strip(ref_data, "%")
        ref_data = ref_data.astype(float)
        
        ref_mean = np.mean(ref_data)
        ref_std = np.std(ref_data)
        ref_max = ref_mean + ref_std
        ref_min = ref_mean - ref_std
        
        reference_data.append((ref_taxon, ref_mean, ref_std, ref_max, ref_min))        
        
normal_parameters = pd.DataFrame.from_records(reference_data, columns = ["taxonomy", "mean", "std", "max", "min"])
normal_parameters = normal_parameters.set_index('taxonomy')

deviant_sample_ids = metadata[filter2]["sample-id"].to_numpy()
deviant_sample_ids = np.array(deviant_sample_ids[0])
deviant_sample_ids = np.insert(deviant_sample_ids, 0, 'taxonomy')

taxonomy_deviant = taxonomy[deviant_sample_ids]

# For our calculations, we choose deviant sample ERR1072712 (no particular reason for choice)
taxonomy_deviant["ERR1072712"] = taxonomy_deviant["ERR1072712"].str.strip("%")
taxonomy_deviant.ERR1072712 = taxonomy_deviant.ERR1072712.astype(float)

taxonomy_deviant = taxonomy_deviant.set_index("taxonomy")

# Categorising the taxons in the selected DEVIANT sample into U (under-represented), O (over-represented), 
# and N (normal abundance) categories as per the population averaged values of the respective taxons from 
# the NORMAL samples

class_list = []

for i in range(len(taxonomy_range)):
    
    ref_taxon = taxonomy_range[i]
    ref_parameters = list(normal_parameters.loc[ref_taxon])
    
    test_val = float(taxonomy_deviant.loc[ref_taxon])
    set_var = ""
    
    if test_val < ref_parameters[3]:
        set_var = "U"
    elif test_val > ref_parameters[2]:
        set_var = "O"
    else:
        set_var = "N"
    
    class_list.append(set_var)

taxonomy_deviant = taxonomy_deviant.assign(set_class=class_list)

nim_aminoacids = pd.read_csv("nim-aminoacids_400.csv")
nim_aminoacids = nim_aminoacids.set_index('taxonomy')

filter3 = ['Trp', 'His', 'Pro', 'Leu', 'Arg', 'Ile_Val']
nim_aminoacids = nim_aminoacids[filter3]

# Reading the Nutrient Impact matrices for various categories of nutrients into pandas dataframes

nim_aminoacidsD = pd.read_csv("nim-aminoacidsD_400.csv")
nim_sugars = pd.read_csv("nim-sugars_400.csv")
nim_vitamins = pd.read_csv("nim-vitamins_400.csv")
nim_aminoacidsD = nim_aminoacidsD.set_index('taxonomy')
nim_sugars = nim_sugars.set_index('taxonomy')
nim_vitamins = nim_vitamins.set_index('taxonomy')

nim_total = nim_aminoacids.join(nim_aminoacidsD)
nim_total = nim_total.join(nim_sugars)
nim_total = nim_total.join(nim_vitamins)

# Entire list of nutrients under consideration for dietary intervention
nutrients_range = list(nim_total)

# Based on the taxon classification into U and O sets, we store the corresponding nutrient impact scores from the
# cumulative nutient impact matrix into two dictionaries for further use

dict_unbalanced_U = {}
dict_unbalanced_O = {}

for i in range(len(taxonomy_range)):
    
    ref_taxon = taxonomy_range[i]
    ref_class = taxonomy_deviant.loc[ref_taxon].set_class
            
    if ref_class == 'U':
        dict_unbalanced_U[ref_taxon] = list(nim_total.loc[ref_taxon])
    elif ref_class == 'O':
        dict_unbalanced_O[ref_taxon] = list(nim_total.loc[ref_taxon])

# Defining the reward function that will be called for the score calculation of various nutrient combinations

def reward(a, dict_unbalanced, list_indices, Epsilon):
    
    ## This function calculates the total score for any taxon (a) with respect to the entire set of nutrients 
    ## under consideration
    
    temp_product = 1
    ref_nim_list = dict_unbalanced[a]
    
    for k in range(len(list_indices)):
        temp_product = temp_product * (1 - ref_nim_list[list_indices[k]])

    if (1-temp_product) >= Epsilon:
        return 1
    else:
        return 0

def reward_nutrient(n, dict_unbalanced_O, dict_unbalanced_U, Epsilon_O, Epsilon_U):
    
    ## This function calculates the total score for any nutrinent (n) with respect to all the under-represented and 
    ## over-represented taxons under consideration
    
    sum_U = 0
    sum_O = 0

    for a in O_keys:

        sum_O = sum_O + reward(a,dict_unbalanced_O,[n], Epsilon_O)

    for b in U_keys:

        sum_U = sum_U + reward(b,dict_unbalanced_U,[n], Epsilon_U)
        
    return (sum_U - sum_O)

# Naive Randomized Algorithm - calculates the score for a random selection of nutrients from the entire set and 
# returns the max score out of all random selections

final_score_dict = {}

# Defining parameters

m = [5,10,15,20,25]
l = len(nutrients_range)

# Thresholds for reward calculation

Epsilon_O = 0.9
Epsilon_U = 0.5

O_keys = list(dict_unbalanced_O.keys())
U_keys = list(dict_unbalanced_U.keys())

for k in m:
    
    score_dict = {}
    temp_max_length = 0
    temp_max_score = -10000
    i = 0
    
    while i < 50000:
        # Calculating score for each m across 50,000 iterations
        
        temp_score_dict = {} 
        
        for j in range(1, k+1):
            
            list_indices = random.sample(range(l), j)

            sum_U = 0
            sum_O = 0

            for a in O_keys:

                sum_O = sum_O + reward(a,dict_unbalanced_O,list_indices, Epsilon_O)

            for b in U_keys:
                sum_U = sum_U + reward(b,dict_unbalanced_U,list_indices, Epsilon_U)

            tempscore = sum_U - sum_O

            temp_score_dict[len(list_indices)] = int(tempscore)

        key_max = max(temp_score_dict, key= lambda x: temp_score_dict[x])
        score_max = temp_score_dict[key_max]
        
        if score_max > temp_max_score:
            temp_max_score = score_max
            temp_max_length = key_max
        
        i +=1
    
    final_score_dict[k] = [temp_max_length, temp_max_score]

# Randomized Local Search Algorithm - calculates the score for a random selection of nutrients from the entire set and 
# then, makes localized changes randomly to improve the score of the random selection

# Defining the parameters

l = len(nutrients_range)

# Thresholds for reward calculation

Epsilon_O = 0.9
Epsilon_U = 0.5

O_keys = list(dict_unbalanced_O.keys())
U_keys = list(dict_unbalanced_U.keys())

final_score_dict = {}

for m in [5,10,15,20]:
    
    scores = []
    naive_scores = []
    i = 0
    
    while i < 5000:
        # Calculating randomised selection score for each m across 50,000 iterations

        temp_score_dict = {} 
        list_indices = random.sample(range(l), m)
        sum_U = 0
        sum_O = 0
        
        for a in O_keys:

            sum_O = sum_O + reward(a,dict_unbalanced_O,list_indices, Epsilon_O)

        for b in U_keys:
            sum_U = sum_U + reward(b,dict_unbalanced_U,list_indices, Epsilon_U)
            
        bestscore = sum_U - sum_O
        naive_scores.append(bestscore)
        
        j = 0
        while j<1000:
        # Making localised improvements in each randomised selection for each m across 10,000 iterations

            list_indices_copy = copy.deepcopy(list_indices)
            
            if len(list_indices)==m:
                random_case = random.sample(range(5), 1)[0]
                if random_case==1: # Deleting an element randomly from the list with a 20% probability
                    random_sample = random.sample(range(len(list_indices)), 1)[0]
                    del list_indices[random_sample]
                else: # Replacing an element randomly from the list with a 80% probability
                    random_nutrient = random.choice(range(len(list_indices)))
                    not_used_nutrients = list( set(list(range(l)))-set(list_indices))
                    random_new_nutrient = random.choice(not_used_nutrients)
                    list_indices[random_nutrient] = random_new_nutrient

            else:
                random_case = random.sample(range(5), 1)[0]
                if random_case==1: # Adding an element randomly to the list with a 20% probability
                    not_used_nutrients = list( set(list(range(l)))-set(list_indices))
                    random_new_nutrient = random.sample(not_used_nutrients, 1)[0]
                    list_indices.append(random_new_nutrient)
                else: # Replacing an element randomly from the list with a 80% probability
                    random_nutrient = random.sample(range(len(list_indices)), 1)[0]
                    not_used_nutrients = list( set(list(range(l)))-set(list_indices))
                    random_new_nutrient = random.choice(not_used_nutrients)
                    list_indices[random_nutrient] = random_new_nutrient

            sum_U = 0
            sum_O = 0

            for a in O_keys:

                sum_O = sum_O + reward(a,dict_unbalanced_O,list_indices, Epsilon_O)

            for b in U_keys:

                sum_U = sum_U + reward(b,dict_unbalanced_U,list_indices, Epsilon_U)

            newscore = sum_U - sum_O

            if bestscore<newscore:
                bestscore = newscore
            else:
                list_indices = list_indices_copy
            j+=1

        scores.append(int(bestscore))

        i +=1
        
    final_score_dict[m] = [np.array(naive_scores).max(), np.array(scores).max()]

# Randomised divide and conquer - calculates the score by dividing entire set of nutrients in clusters randomly and
# picking a single nutrient (with the best score) from each cluster to constitute a combination, 
# for which the score is then calculated

#Defining the parameters

l = len(nutrients_range)
final_score_dict = {}

# Thresholds for reward calculation
Epsilon_O = 0.9
Epsilon_U = 0.5

O_keys = list(dict_unbalanced_O.keys())
U_keys = list(dict_unbalanced_U.keys())

for m in [5,10,15,20,25]:
    best_of_all = -10000
    each_cluster_size = l//m+1
    
    i=0
    while i<50000:

        permute_nutrients = np.random.permutation(range(l))
        all_indices = []
        
        for j in range(m):

            permute_partition = permute_nutrients[j*each_cluster_size:(j+1)*each_cluster_size]
            best_tmpscore = -100000
            best_index = -10
            
            for n in permute_partition:
                score_tmp = reward_nutrient(n,dict_unbalanced_O,dict_unbalanced_U, Epsilon_O,Epsilon_U)

                if best_tmpscore<score_tmp:
                    best_tmpscore=score_tmp
                    best_index = n
                    
            all_indices.append(best_index)

        sum_U = 0
        sum_O = 0
        
        for a in O_keys:

            sum_O = sum_O + reward(a,dict_unbalanced_O,all_indices, Epsilon_O)

        for b in U_keys:

            sum_U = sum_U + reward(b,dict_unbalanced_U,all_indices, Epsilon_U)
            
        new_best = sum_U-sum_O
        
        if new_best>best_of_all:  
            best_of_all=new_best
        i+=1
            
    final_score_dict[m] = best_of_all

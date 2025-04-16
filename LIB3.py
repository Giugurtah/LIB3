import ctypes as ct
import numpy as np
import pandas as pd
import math
import json 
import bisect
import matplotlib.pyplot as plt
#from matplotlib.patches import Rectangle, Circle, Polygon
from pathlib import Path

# C types declaration
INT = ct.c_int
PINT = ct.POINTER(INT)
FLO = ct.c_double
PFLO = ct.POINTER(FLO)
PPFLO = ct.POINTER(PFLO)

C_SCRIPT = ct.CDLL(str(Path(__file__).parent.absolute()) + '/LIB3C.so')

gpi_c = getattr(C_SCRIPT,'gpi_c')
gpi_c.restype = FLO

class Node:
    def __init__(self, 
                 gpi=None, 
                 ppi=None, 
                 position=None, 
                 impurity=None,
                 feature=None, 
                 treshold=None, 
                 left=None, 
                 right=None,
                 LIFT_1=None,
                 LIFT_2=None,
                 GCR=None,
                 distribution=None,
                 N=None,
                 labels=None,
                 strat_labels=None,
                 *,value=None):
        
        self.gpi = gpi
        self.ppi = ppi
        self.position = position
        self.feature = feature
        self.treshold = treshold
        self.impurity = impurity
        self.left = left
        self.right = right
        self.value = value
        self.distribution = distribution
        self.N = N
        self.labels = labels
        self.LIFT_1 = LIFT_1
        self.LIFT_2 = LIFT_2
        self.GCR = GCR
        self.strat_labels = strat_labels

    def _is_leaf_node(self):
        return self.value is not None
    
class Tree:
    def __init__(self, 
                 min_ppi = 0.01,
                 min_gpi=0.1, 
                 min_impurity=0.1,
                 compound_feats=False, 
                 impurity_method="gini", 
                 min_samples_split=1, 
                 max_depth=100,
                 feats_viewed=1, 
                 model="twoStage",
                 homogeneity=None,
                 FAST=False):
        
        col_config = {
            "slba": ['id', 'Node Type', 'Value', 'Splitting Variable', 'n', 'Class Distribution', 'Alpha', 'Beta', 'LIFT_K1', 'LIFT_K2', 'GCR', 'Treshold', 'Impurity', 'Gpi', 'Ppi'],
            "lba": ['id', 'Node Type', 'Value', 'Splitting Variable', 'n', 'Class Distribution', 'Alpha', 'Beta', 'LIFT_K1', 'LIFT_K2', 'GCR', 'Treshold', 'Impurity', 'Gpi', 'Ppi'],
            "twoStage": ['id', 'Node Type', 'Value','Splitting Variable', 'n', 'Class Distribution', 'Treshold', 'Impurity', 'Gpi', 'Ppi'],
            "twoing": ['id', 'Node Type', 'Value', 'Splitting Variable', 'n', 'Class Distribution', 'Treshold', 'Impurity', 'Gpi', 'Ppi'],
        }
        col = col_config.get(model, ['id']) 


        self.min_samples_split = min_samples_split
        self.max_depth = max_depth
        self.feats_viewed = feats_viewed
        self.impurity_method = impurity_method
        self.min_gpi = min_gpi
        self.min_ppi = min_ppi
        self.model_name = model
        self.compound_feats = compound_feats
        self.min_impurity = min_impurity
        self.FAST = FAST
        self.homogeneity = homogeneity

        if model == "slba" and homogeneity is not None:
            self.model = getattr(C_SCRIPT, model + "_" + homogeneity + "_c")
        else:
            self.model = getattr(C_SCRIPT, model + "_c")

        self.results = pd.DataFrame(columns=col)
        self.targhet_dist = None
        self.depth = None
        self.l_c = None
        self.r_c = None
        self.root=None

    def fit(self, X, y, x_s = None):
        if(self.impurity_method!="entropy" and self.impurity_method!="gini" and self.impurity_method!="error"):
            print("ERROR. An impurity method that is not gini, error or entropy has been given.\nA gini impurity method has been automatichally given")
            self.impurity_method = "gini"

        self.targhet_dist = [np.unique(y), np.unique(y, return_counts=True)[1]/len(y)]

        if self.model_name == "slba" and x_s is None:
            print("A Simultaneous Latent Budget Analysis was selected but no stratifying variable was given")
            return
        
        if self.model_name == "slba" and x_s is not None:
            self.root = self._grow_tree(X, y, x_s)
            return

        if self.model_name != "slba" and x_s is not None:
            print("A stratifying variable has been given but the Simultaneous Latent Budget Analysis model was not selcted.")
            print("The selected model will be regularly performed without taking in consideration the stratifying variable.")

        #The recursive function is called
        self.root = self._grow_tree(X, y)
        self._pruning(self.root)

    def predict(self, X):
                arr = []
                count = 1
                for i in range(len(X)):
                    print("Riga ispezionata:", count)
                    arr.append(self._traverse_tree(X.iloc[i], self.root))
                    count = count+1
                
                return np.array(arr)

    #Growing functions
    def _grow_tree(self, X, y, x_s=None, depth=0, pos=1):
        print("INIZIO CALCOLO NODO ALLA PROFONDITA':", depth, " ID:", pos) #TODO da cancellare
        # All columns with constant values are dropped
        nunique = X.nunique()
        col_to_drop = nunique[nunique == 1].index
        X = X.drop(col_to_drop, axis=1)

        compound_feature = "" #Compound feature used only if compound=TRUE

        # Generale sizes of X and y
        n_samples, n_feats = X.shape
        n_labels = len(np.unique(y))
        impurity = self._impurity(y)
        distribution =  np.unique(y, return_counts=True)[1]/len(y)
        print("Il sottoinsieme del dataset contiene attualmenete", n_samples, "samples,", n_feats, "feats e", n_labels, "labels") #TODO da cancellare

        # Check the stopping criteria
        if(depth>=self.max_depth or n_labels==1 or n_samples<self.min_samples_split or n_feats==0 or impurity<self.min_impurity):
            leaf_value = y.mode()[0]
            gcr = self._get_gcr(distribution, np.unique(y))
            self._update_plot(depth, pos)
            self._update_results(pos, "Leaf Node", leaf_value, None, len(y), distribution , None, None, None, impurity, None, None, None, None, gcr)

            print("(criteria1)  Nodo foglia raggiunto. Valore di y associato a questo nodo:", leaf_value) #TODO da cancellare
            print("\n")
            return Node(position=pos, value=leaf_value, impurity=impurity, distribution=distribution, N=len(y), labels=np.unique(y), GCR=gcr)

        # The best predictors are found (highest gpi)
        gpi, array_of_Fs = self._gpi(X, y, n_samples)
        gpi_index = np.arange(n_feats)
        gpi, gpi_index = zip(*sorted(zip(gpi, gpi_index), reverse=True))
        if(self.compound_feats and n_feats > 1):
            #This section generates a compound variable using the two variables with the highest gpi
            feature_1 = X.columns[gpi_index[0]]
            feature_2 = X.columns[gpi_index[1]]
            compound_feature = feature_1 + "_" + feature_2
            X[compound_feature] = X[feature_1].astype(str) + "_" + X[feature_2].astype(str)

            F = pd.crosstab(X.loc[:, compound_feature], y, margins=False)
            I = len(F)
            J = len(F.iloc[0])

            FLOARR = FLO * J
            PFLOARR = PFLO * I

            ptr_rc = PFLOARR()
            # Array of pointers initialization
            for i in range(I):
                ptr_rc[i] = FLOARR()
                for j in range(J):
                    ptr_rc[i][j] = F.iloc[i].iloc[j]/n_samples

            array_of_Fs.append(ptr_rc)
            gpi = np.append(gpi, gpi_c(I, J, ptr_rc))
            gpi_index = np.append(gpi_index, n_feats)
            gpi, gpi_index = zip(*sorted(zip(gpi, gpi_index), reverse=True))

        # Amongst the best predictors the one with the highest improvement (highest ppi) is chosen to perform the split
        best_feature, best_treshold, best_ppi, best_idx, alpha, beta = self._find_best_predictor(X, y, x_s, n_labels, array_of_Fs, gpi_index, gpi)

        print("best feature", best_feature)
        print("gpi", gpi[best_idx])
        print("best ppi", best_ppi)

        # Check the stopping criteria
        if(gpi[best_idx]<self.min_gpi or best_ppi<self.min_ppi or best_ppi == 0):
            leaf_value = y.mode()[0]
            gcr = self._get_gcr(distribution, np.unique(y))
            self._update_plot(depth, pos)
            self._update_results(pos, "Leaf Node", leaf_value, None, len(y), distribution, None, None, None, impurity, None, None, None, None, gcr=gcr)

            print("(criteria2) Nodo foglia raggiunto. Valore di y associato a questo nodo:", leaf_value) #TODO da cancellare
            print("\n")
            return Node(position=pos, value=leaf_value, impurity=impurity, distribution=distribution, N=len(y), labels=np.unique(y), GCR=gcr)
        
        if x_s is not None and (self.homogeneity is None or self.homogeneity == "B"): #The homogeneity is None or B then A is composed by t matrices
            local_treshold = []
            for t in range(len(np.unique(x_s))):
                local_treshold.append(np.unique(X[best_feature])[np.array(best_treshold[t])[:] > 0])
            best_treshold = local_treshold
            print("best treshold: ", best_treshold, "\n")

            indexL, indexR = self._splitS(X[best_feature], x_s, best_treshold)
        else: #The homogeneity is A or AB, or the model is not slba then A is a single matrix
            print("best treshold: ", np.unique(X[best_feature])[np.array(best_treshold)[:] > 0], "\n")

            indexL, indexR = self._split(X[best_feature], best_treshold)
        
        if x_s is not None and (self.homogeneity is None or self.homogeneity == "A"): #The homogeneity is None or A then B is composed by t matrices
            lift1 = []
            lift2 = []
            for t in range(len(np.unique(x_s))):
                lift1.append([x[0] for x in beta[t]]/distribution)
                lift2.append([x[1] for x in beta[t]]/distribution)

        if x_s is not None and self.homogeneity != "AB":
            if self.homogeneity is None: #The homogeneity is None then both A and B are composed by t matrices
                self._update_results(pos, "Parent Node", None, best_feature, len(y), distribution, best_treshold, gpi[best_idx], best_ppi, impurity, alpha, beta, lift1, lift2, [0 for _ in range(len(distribution))])
            if self.homogeneity == "A": #The homogeneity is A then B is composed by t matrices
                self._update_results(pos, "Parent Node", None, best_feature, len(y), distribution, np.unique(X[best_feature])[np.array(best_treshold)[:] > 0], gpi[best_idx], best_ppi, impurity, alpha, beta, lift1, lift2, [0 for _ in range(len(distribution))])
            if self.homogeneity == "B": #The homogeneity is B then A is composed by t matrices
                self._update_results(pos, "Parent Node", None, best_feature, len(y), distribution, best_treshold, gpi[best_idx], best_ppi, impurity, alpha, beta, [x[0] for x in beta]/distribution, [x[1] for x in beta]/distribution, [0 for _ in range(len(distribution))])
        else: #The homogeneity is AB or the model is not slba then both A and B are single matrices
            self._update_results(pos, "Parent Node", None, best_feature, len(y), distribution, np.unique(X[best_feature])[np.array(best_treshold)[:] > 0], gpi[best_idx], best_ppi, impurity, alpha, beta, [x[0] for x in beta]/distribution, [x[1] for x in beta]/distribution, [0 for _ in range(len(distribution))])
            
        # Create the child nodes
        if(best_feature == compound_feature):
            X = X.loc[:, X.columns != feature_1]
            X = X.loc[:, X.columns != feature_2]
            left = self._grow_tree(X.loc[indexL, X.columns != best_feature], y[indexL], depth+1, 2*pos)
            right = self._grow_tree(X.loc[indexR, X.columns != best_feature], y[indexR], depth+1, 2*pos+1)
        else:
            X = X.loc[:, X.columns != compound_feature]
            if x_s is not None:
                left = self._grow_tree(X.loc[indexL, X.columns != best_feature], y[indexL], x_s[indexL], depth+1, 2*pos)
                right = self._grow_tree(X.loc[indexR, X.columns != best_feature], y[indexR], x_s[indexR], depth+1, 2*pos+1)
            else:
                left = self._grow_tree(X.loc[indexL, X.columns != best_feature], y[indexL], None, depth+1, 2*pos)
                right = self._grow_tree(X.loc[indexR, X.columns != best_feature], y[indexR], None, depth+1, 2*pos+1)

        if x_s is not None:
            if self.homogeneity is None: #The homogeneity is None then both A and B are composed by t matrices
                return Node(gpi=gpi[best_idx], ppi=best_ppi, position=pos, feature=best_feature, 
                            treshold=best_treshold, 
                            left=left, right=right, impurity=impurity, distribution=distribution, 
                            N=len(y), labels=np.unique(y),
                            LIFT_1=lift1, LIFT_2=lift2, GCR=None, strat_labels=np.unique(x_s))
            if self.homogeneity == "A": #The homogeneity is A then only B is composed by t matrices
                return Node(gpi=gpi[best_idx], ppi=best_ppi, position=pos, feature=best_feature, 
                            treshold=np.unique(X[best_feature])[np.array(best_treshold)[:] > 0], 
                            left=left, right=right, impurity=impurity, distribution=distribution, 
                            N=len(y), labels=np.unique(y),
                            LIFT_1=lift1, LIFT_2=lift2, GCR=None, strat_labels=np.unique(x_s))
            if self.homogeneity == "B": #The homogeneity is B then only A is composed by t matrices
                return Node(gpi=gpi[best_idx], ppi=best_ppi, position=pos, feature=best_feature, 
                            treshold=best_treshold, 
                            left=left, right=right, impurity=impurity, distribution=distribution, 
                            N=len(y), labels=np.unique(y),
                            LIFT_1=[x[0] for x in beta]/distribution, LIFT_2=[x[1] for x in beta]/distribution, GCR=None, strat_labels=np.unique(x_s))
            if self.homogeneity == "AB": #The homogeneity is AB then both A and B are single matrices
                return Node(gpi=gpi[best_idx], ppi=best_ppi, position=pos, feature=best_feature, 
                        treshold=np.unique(X[best_feature])[np.array(best_treshold)[:] > 0], 
                        left=left, right=right, impurity=impurity, distribution=distribution, 
                        N=len(y), labels=np.unique(y),
                        LIFT_1=[x[0] for x in beta]/distribution, LIFT_2=[x[1] for x in beta]/distribution, GCR=None, strat_labels=np.unique(x_s))
        else:
            #The model is not slba
            return Node(gpi=gpi[best_idx], ppi=best_ppi, position=pos, feature=best_feature, 
                        treshold=np.unique(X[best_feature])[np.array(best_treshold)[:] > 0], 
                        left=left, right=right, impurity=impurity, distribution=distribution, 
                        N=len(y), labels=np.unique(y),
                        LIFT_1=[x[0] for x in beta]/distribution, LIFT_2=[x[1] for x in beta]/distribution, GCR=None)

    def _split(self, X_col, best_treshold):
        L = np.unique(X_col)[np.array(best_treshold)[:] > 0]
        indexL = np.isin(X_col, L)
        indexR =  np.isin(X_col, L) == False
        return indexL, indexR

    def _splitS(self, X_col, x_s, best_treshold):
        local_X = X_col.astype(str) + "_" + x_s.astype(str)
        L = []
        for i in range(len(np.unique(x_s))):
            for j in range(len(best_treshold[i])):
                L.append(str(best_treshold[i][j]) + "_" + str(np.unique(x_s)[i]))

        indexL = np.isin(local_X, L)
        indexR =  np.isin(local_X, L) == False
        return indexL, indexR

    def _find_best_predictor(self, X, y, x_s, n_labels, array_of_Fs, gpi_i, gpi):
        best_ppi = 0
        best_treshold = []
        best_alpha = []
        best_beta = []
        best_feature = ""
        best_idx = 0

        for i in range(min(self.feats_viewed, len(gpi))):
            index = gpi_i[i]
            current_feature = X.columns[index]
            n_mod = len(np.unique(X[current_feature]))

            print("Current feature:", current_feature)
            if x_s is not None:
                n_mod_strat = len(np.unique(x_s))
                N = len(X[current_feature])

                #Inizialization of F 
                x_mods = sorted(set(X[current_feature]))
                y_mods = sorted(set(y))

                # A temp dataframe is initialized
                df = pd.DataFrame({
                    'x': pd.Categorical(X[current_feature], categories=x_mods),
                    'y': pd.Categorical(y, categories=y_mods),
                    'x_s': x_s
                })

                PPFLOARRF = PPFLO * (n_mod_strat)
                F = PPFLOARRF()
                for idx, strat in enumerate(df['x_s'].unique()):
                    sub_df = df[df['x_s'] == strat]

                    # The contingency table is evaluated
                    F_local = pd.crosstab(sub_df['x'], sub_df['y'], dropna=False)
                    I = len(F_local)
                    J = len(F_local.iloc[0])

                    FLOARRF = FLO * J
                    PFLOARRF = PFLO * I
                    F[idx] = PFLOARRF()
                    # Array of pointers initialization
                    for i in range(I):
                        F[idx][i] = FLOARRF()
                        for j in range(J):
                            F[idx][i][j] = F_local.iloc[i].iloc[j]/N

                PPFLOARR = PPFLO * (n_mod_strat)
                FLOARR = FLO * 2
                PFLOARRA = PFLO * (n_mod)
                PFLOARRB = PFLO * (n_labels)

                if self.homogeneity is None:
                    # Current treshold will hold the best treshold for the current feature
                    PPFLOTRESH = PFLO * (n_mod_strat + 1)
                    PFLOTRESH = FLO * (n_mod)
                    PFLOTRESH1 = FLO * (1)

                    current_treshold = PPFLOTRESH()
                    current_treshold[0] = PFLOTRESH1()
                    current_treshold[0][0] = 0.0
                    for t in range(n_mod_strat):
                        current_treshold[t+1] = PFLOTRESH()
                        # Array of pointers initialization
                        for k in range(n_mod):
                            current_treshold[t+1][k] = 0.0

                    # Initialization of A
                    alpha = PPFLOARR()
                    for t in range(n_mod_strat):
                        alpha[t] = PFLOARRA()
                        # Array of pointers initialization
                        for k in range(n_mod):
                            alpha[t][k] = FLOARR()
                            for j in range(2):
                                alpha[t][k][j] = 0.0

                    # Initialization of B
                    beta = PPFLOARR()
                    for t in range(n_mod_strat):
                        beta[t] = PFLOARRB()
                        # Array of pointers initialization
                        for k in range(n_labels):
                            beta[t][k] = FLOARR()
                            for j in range(2):
                                beta[t][k][j] = 0.0

                    self.model(n_mod, n_labels, n_mod_strat, array_of_Fs[index], F, current_treshold, alpha, beta)
                    print("\n")

                    if (current_treshold[0][0] > best_ppi):
                        best_idx = index
                        best_treshold = []
                        best_alpha = []
                        best_beta = []

                        best_ppi = current_treshold[0][0]
                        best_feature = current_feature
                        for t in range(n_mod_strat):
                            local_treshold = []
                            local_alpha = np.zeros((n_mod, 2))
                            local_beta = np.zeros((n_labels, 2))

                            for k in range(n_mod):
                                local_treshold.append(current_treshold[t+1][k])
                            best_treshold.append(local_treshold)

                            for k in range(2):
                                for j in range(n_mod):
                                    local_alpha[j][k] = alpha[t][j][k]
                                for j in range(n_labels):
                                    local_beta[j][k] = beta[t][j][k]
                            best_alpha.append(local_alpha)
                            best_beta.append(local_beta)

                elif self.homogeneity == "A":
                    # Current treshold will hold the best treshold for the current feature
                    PPFLOTRESH = PFLO * (2)
                    PFLOTRESH = FLO * (n_mod)
                    PFLOTRESH1 = FLO * (1)

                    current_treshold = PPFLOTRESH()
                    current_treshold[0] = PFLOTRESH1()
                    current_treshold[0][0] = 0.0
                    current_treshold[1] = PFLOTRESH()
                    for k in range(n_mod):
                        current_treshold[1][k] = 0.0

                    # Initialization of A
                    alpha = PFLOARRA()
                    for k in range(n_mod):
                        alpha[k] = FLOARR()
                        for j in range(2):
                            alpha[k][j] = 0.0

                    # Initialization of B
                    beta = PPFLOARR()
                    for t in range(n_mod_strat):
                        beta[t] = PFLOARRB()
                        # Array of pointers initialization
                        for k in range(n_labels):
                            beta[t][k] = FLOARR()
                            for j in range(2):
                                beta[t][k][j] = 0.0

                    self.model(n_mod, n_labels, n_mod_strat, array_of_Fs[index], F, current_treshold, alpha, beta)
                    print("\n")

                    if (current_treshold[0][0] > best_ppi):
                        best_idx = index
                        best_treshold = []
                        best_alpha = np.zeros((n_mod, 2))
                        best_beta = []

                        best_ppi = current_treshold[0][0]
                        best_feature = current_feature

                        for j in range(n_mod):
                            best_alpha[j][0] = alpha[j][0]
                            best_alpha[j][1] = alpha[j][1]
                            best_treshold.append(current_treshold[1][j])

                        for t in range(n_mod_strat):
                            local_beta = np.zeros((n_labels, 2))
                            for k in range(2):
                                for j in range(n_labels):
                                    local_beta[j][k] = beta[t][j][k]
                            best_beta.append(local_beta)

                elif self.homogeneity == "B":
                    # Current treshold will hold the best treshold for the current feature
                    PPFLOTRESH = PFLO * (n_mod_strat + 1)
                    PFLOTRESH = FLO * (n_mod)
                    PFLOTRESH1 = FLO * (1)

                    current_treshold = PPFLOTRESH()
                    current_treshold[0] = PFLOTRESH1()
                    current_treshold[0][0] = 0.0
                    for t in range(n_mod_strat):
                        current_treshold[t+1] = PFLOTRESH()
                        # Array of pointers initialization
                        for k in range(n_mod):
                            current_treshold[t+1][k] = 0.0

                    # Initialization of B
                    beta = PFLOARRB()
                    for k in range(n_labels):
                        beta[k] = FLOARR()
                        for j in range(2):
                            beta[k][j] = 0.0

                    # Initialization of A
                    alpha = PPFLOARR()
                    for t in range(n_mod_strat):
                        alpha[t] = PFLOARRA()
                        # Array of pointers initialization
                        for k in range(n_mod):
                            alpha[t][k] = FLOARR()
                            for j in range(2):
                                alpha[t][k][j] = 0.0

                    self.model(n_mod, n_labels, n_mod_strat, array_of_Fs[index], F, current_treshold, alpha, beta)
                    print("\n")

                    if (current_treshold[0][0] > best_ppi):
                        best_idx = index
                        best_treshold = []
                        best_alpha = []
                        best_beta = np.zeros((n_labels, 2))

                        best_ppi = current_treshold[0][0]
                        best_feature = current_feature
                        for j in range(n_labels):
                            best_beta[j][0] = beta[j][0]
                            best_beta[j][1] = beta[j][1]

                        for t in range(n_mod_strat):
                            local_treshold = []
                            local_alpha = np.zeros((n_mod, 2))
                            local_beta = np.zeros((n_labels, 2))

                            for k in range(n_mod):
                                local_treshold.append(current_treshold[t+1][k])
                            best_treshold.append(local_treshold)

                            for k in range(2):
                                for j in range(n_mod):
                                    local_alpha[j][k] = alpha[t][j][k]
                            best_alpha.append(local_alpha)    

                elif self.homogeneity == "AB":      
                    # Current treshold will hold the best treshold for the current feature
                    PPFLOTRESH = PFLO * (2)
                    PFLOTRESH = FLO * (n_mod)
                    PFLOTRESH1 = FLO * (1)

                    current_treshold = PPFLOTRESH()
                    current_treshold[0] = PFLOTRESH1()
                    current_treshold[0][0] = 0.0
                    current_treshold[1] = PFLOTRESH()
                    for k in range(n_mod):
                        current_treshold[1][k] = 0.0

                    # Initialization of A
                    alpha = PFLOARRA()
                    for k in range(n_mod):
                        alpha[k] = FLOARR()
                        for j in range(2):
                            alpha[k][j] = 0.0

                    # Initialization of B
                    beta = PFLOARRB()
                    for k in range(n_labels):
                        beta[k] = FLOARR()
                        for j in range(2):
                            beta[k][j] = 0.0

                    self.model(n_mod, n_labels, n_mod_strat, array_of_Fs[index], F, current_treshold, alpha, beta)
                    print("\n")

                    if (current_treshold[0][0] > best_ppi):
                        best_idx = index
                        best_treshold = []
                        best_alpha = np.zeros((n_mod, 2))
                        best_beta = np.zeros((n_labels, 2))

                        best_ppi = current_treshold[0][0]
                        best_feature = current_feature

                        for j in range(n_mod):
                            best_alpha[j][0] = alpha[j][0]
                            best_alpha[j][1] = alpha[j][1]
                            best_treshold.append(current_treshold[1][j])

                        for i in range(n_labels):
                            best_beta[i][0] = beta[i][0]
                            best_beta[i][1] = beta[i][1]

            else:
                # Current treshold will hold the best treshold for the current feature
                PFLOTRESH = FLO * (n_mod+1)
                current_treshold = PFLOTRESH()
                # Array of pointers initialization
                for k in range(n_mod + 1):
                    current_treshold[k] = 0.0

                # Alpha and Beta will hold respectivelly the mixing parameters and the latent budgets. 
                # If the model is not lba they remain unused
                FLOARR = FLO * 2
                PFLOARRA = PFLO * (n_mod)
                PFLOARRB = PFLO * (n_labels)

                alpha = PFLOARRA()
                # Array of pointers initialization
                for k in range(n_mod):
                    alpha[k] = FLOARR()
                    for j in range(2):
                        alpha[k][j] = 0.0

                beta = PFLOARRB()
                # Array of pointers initialization
                for k in range(n_labels):
                    beta[k] = FLOARR()
                    for j in range(2):
                        beta[k][j] = 0.0
                
                self.model(n_mod, n_labels, array_of_Fs[index], current_treshold, alpha, beta)

                if (current_treshold[0] > best_ppi):
                    best_idx = i
                    best_treshold = []
                    best_alpha = np.zeros((n_mod, 2))
                    best_beta = np.zeros((n_labels, 2))
                    best_ppi = current_treshold[0]
                    best_feature = current_feature
                    for k in range(1,n_mod+1):
                        best_treshold.append(current_treshold[k]/best_ppi)
                    for k in range(2):
                        for j in range(n_mod):
                            best_alpha[j][k] = alpha[j][k]
                        for j in range(n_labels):
                            best_beta[j][k] = beta[j][k]
            
            if(self.FAST and i<len(gpi)-1 and best_ppi>gpi[i+1]):
                return best_feature, best_treshold, best_ppi, best_idx, best_alpha, best_beta
            
        return best_feature, best_treshold, best_ppi, best_idx, best_alpha, best_beta
    
    def _pruning(self, node):
        if node is None or node.value is not None:
            return node
        
        # Applica il pruning ai sotto-alberi
        node.left = self._pruning(node.left)
        node.right = self._pruning(node.right)

        if node.left and node.right and node.left.value is not None and node.right.value is not None:
            if node.left.value == node.right.value:
                node.gpi = None
                node.ppi = None
                node.feature = None
                node.treshold = None
                node.LIFT_1 = None
                node.LIFT_2 = None
                node.GCR = self._get_gcr(node.distribution, node.labels)
                node.value = node.left.value

                self.results.loc[self.results["id"] == node.position,"Node Type"] = "Leaf Node"
                self.results.loc[self.results["id"] == node.position,"Value"] = node.value
                self.results.loc[self.results["id"] == node.position,"Splitting Variable"] = None
                self.results.loc[self.results["id"] == node.position,"Treshold"] = None
                self.results.loc[self.results["id"] == node.position,"Gpi"] = None
                self.results.loc[self.results["id"] == node.position,"Ppi"] = None

                if self.model_name == "lba":
                    self.results.at[self.results.index[self.results["id"] == node.position][0], "GCR"] = node.GCR
                    self.results.loc[self.results["id"] == node.position,"Alpha"] = None
                    self.results.loc[self.results["id"] == node.position,"Beta"] = None
                    self.results.loc[self.results["id"] == node.position,"LIFT_K1"] = None
                    self.results.loc[self.results["id"] == node.position,"LIFT_K2"] = None

                self.results = self.results[self.results["id"] != node.left.position]
                self.results = self.results[self.results["id"] != node.right.position]

                node.left = None
                node.right = None

        return node

    #Update functions  
    def _get_gcr(self, distribution, labels):
        gcr = [0 for _ in range(len(labels))]
        for i in range(len(labels)):
            for j in range(len(self.targhet_dist[0])):
                if(self.targhet_dist[0][j] == labels[i]):
                    gcr[i] = (distribution[i]/self.targhet_dist[1][j])
        return gcr

    def _traverse_tree(self, x, node):
        if node._is_leaf_node():
            print("Giunto alla foglia. Valore=", node.value)
            return node.value
        
        if x[node.feature] in node.treshold:
            return self._traverse_tree(x, node.left)
        else:
            return self._traverse_tree(x, node.right)

    def _update_plot(self, depth, pos):
        if(self.depth is None or depth>self.depth):
                self.depth = depth
        var = pos
        l_c = 0
        r_c = 0
        incr = 1/pow(2,depth)
        while var>1:
            l_c += (-2*(var%2) + 1)*(incr)
            r_c += (2*(var%2) - 1)*(incr)
            incr = incr*2
            var = var//2
        if(self.l_c is None or l_c>self.l_c):
            self.l_c = l_c
        if(self.r_c is None or r_c>self.r_c):
            self.r_c = r_c
    
    def _update_results(self, id, type, value, feature, n, distribution, treshold, gpi, ppi, impurity, a, b, lift1, lift2, gcr):
        if(self.model_name == "lba" or self.model_name == "slba"):
            self.results.loc[len(self.results)] = [id, type, value, feature, n, distribution, a, b, lift1, lift2, gcr, treshold, impurity, gpi, ppi]
        else:
            self.results.loc[len(self.results)] = [id, type, value, feature, n, distribution, treshold, impurity, gpi, ppi]

    #Metrics functions
    def _gpi(self, X, y, N): 
            gpi = []
            arr = []

            for x in X:
                F = pd.crosstab(X.loc[:, x], y, margins=False)
                I = len(F)
                J = len(F.iloc[0])

                FLOARR = FLO * J
                PFLOARR = PFLO * I
                
                ptr_rc = PFLOARR()
                # Array of pointers initialization
                for i in range(I):
                    ptr_rc[i] = FLOARR()
                    for j in range(J):
                        ptr_rc[i][j] = F.iloc[i].iloc[j]/N
                arr.append(ptr_rc)
                gpi.append(gpi_c(I, J, ptr_rc))
            return gpi, arr

    def _impurity(self, y):
        impurity = 0
        dist = np.unique(y, return_counts=True)[1]/len(y)
        if self.impurity_method == "gini":
            impurity = 1
            for x in dist:
                impurity -= x*x
            return impurity
        elif self.impurity_method == "entropy":
            for x in dist:
                impurity -= x*math.log10(x)
            return impurity
        else:
            impurity = max(dist)
        return impurity

    #Plotting functions
    def plot(self):
        rate = (1)/(self.l_c + self.r_c)
        x_center = rate*self.l_c
        y_center = 1-0.53
        x_dist = rate*0.5
        y_dist = (0.5-0.06)/self.depth

        fig, ax = plt.subplots()
        self._add_node(x_center, y_center, x_dist, y_dist, self.root, ax)
        #ax.set_aspect('equal')
        plt.grid(False)
        plt.axis(False)

        annotation = ax.annotate(
            text='',
            xy=(0, 0),
            xytext=(15, 15), # distance from x, y
            textcoords='offset points',
            bbox={'boxstyle': 'round', 'fc': 'w'},
            arrowprops={'arrowstyle': '->'}
        )

        annotation.set_visible(False)
        #fig.canvas.mpl_connect("pick_event", self._show_annotation)
        plt.show()
        return
        
    def _add_node(self, x_c, y_c, x_dist, y_dist, node, ax):
        if node._is_leaf_node():
            scatter = plt.scatter(x_c, y_c, edgecolor='blue', facecolor='lightblue', marker='s', zorder=2, s=100,  picker=True)
            return
        
        plt.plot([x_c,x_c-x_dist],[y_c,y_c-y_dist],[x_c,x_c+x_dist],[y_c,y_c-y_dist], color="blue")
        scatter = plt.scatter(x_c, y_c, edgecolor='blue', facecolor='lightblue', zorder=2, s=100,  picker=True)

        self._add_node(x_c-x_dist, y_c-y_dist, x_dist/2, y_dist, node.left, ax)
        self._add_node(x_c+x_dist, y_c-y_dist, x_dist/2, y_dist, node.right, ax)
        return

    def _custom_round(self, value, decimals=2):
        if abs(value) >= 0.01:
            return round(value, decimals)
        
        i = 3
        while i>0:
            treshold = pow(10, -1*i)
            if value >= treshold:
                return round(value, i)
            i +=1 

    def _recurse(self, node):
        print("pos: ", node.position)
        dist = ", ".join(f"{label}={self._custom_round(prob, 2)}" for label, prob in zip(node.labels.tolist(), node.distribution.tolist()))
        if node.value is not None:
            if(self.model_name == "lba" or self.model_name == "slba"):
                gcr = ", ".join(f"{label}={self._custom_round(prob, 2)}" for label, prob in zip(node.labels.tolist(), node.GCR))
            else:
                gcr = None
            return {
                    "isLeaf": 1,
                    "distribution" : dist,
                    "distArray": node.distribution.tolist(),
                    "position": node.position,
                    "value": str(node.value),
                    "impurity": node.impurity,
                    "labels": node.N,
                    "labArray": node.labels.tolist(),
                    "gcr" : gcr,
                    }
        
        if(self.model_name == "lba" or (self.model_name == "slba" and (self.homogeneity == "B" or self.homogeneity == "AB"))):
            lift1 = ", ".join(f"{label}={self._custom_round(prob, 2)}" for label, prob in zip(node.labels.tolist(), node.LIFT_1.tolist()))
            lift2 = ", ".join(f"{label}={self._custom_round(prob, 2)}" for label, prob in zip(node.labels.tolist(), node.LIFT_2.tolist()))
        elif(self.model_name == "slba" and (self.homogeneity == "A" or self.homogeneity is None)):
            lift1 = "\n".join(
                    f"{strat_label}: (" + 
                    ", ".join(f"{label}={self._custom_round(prob, 2)}" for label, prob in zip(node.labels, lift_vals)) + 
                    ")" 
                    for strat_label, lift_vals in zip(node.strat_labels, node.LIFT_1)
                )
            lift2 = "\n".join(
                    f"{strat_label}: (" + 
                    ", ".join(f"{label}={self._custom_round(prob, 2)}" for label, prob in zip(node.labels, lift_vals)) + 
                    ")" 
                    for strat_label, lift_vals in zip(node.strat_labels, node.LIFT_2)
                )
        else:
            lift1 = None
            lift2 = None

        if(self.model_name == "slba" and (self.homogeneity == "A" or self.homogeneity == "AB")) or (self.model_name != "slba"):
            treshold = ', '.join(str(x) for x in node.treshold)
        else:
            treshold =  ", ".join(
                        f"{label}:[{', '.join(str(x) for x in values)}]"
                        for label, values in zip(node.strat_labels, node.treshold)
                        )
        return {
            "isLeaf": 0,
            "feature": node.feature,
            "distribution" : dist,
            "distArray": node.distribution.tolist(),
            "lift1" : lift1,
            "lift2" : lift2,
            "treshold" : treshold,
            "position": node.position,
            "gpi" : node.gpi,
            "ppi" : node.ppi,
            "impurity": node.impurity,
            "labels": node.N,
            "labArray": node.labels.tolist(),
            "children": [
                self._recurse(node.left),
                self._recurse(node.right),
            ]
        }
        
    def plot_html(self, output_file = "tree_visualization.html", title="Decision Tree Visualization",
                  color_palet = ["#FF6384", "#36A2EB", "#FFCE56", "#4CAF50", "#9966FF", "#795548", "#D81B60", "#00ACC1", "#8D6E63", "#FF9800"]):
        tree_JSON = json.dumps(self._recurse(self.root), indent=4)
        dataPlot = {
            "l_c": self.l_c,
            "r_c": self.r_c,
            "depth": self.depth,
            "colors": color_palet
        }
        plot_JSON = json.dumps(dataPlot, indent=4)
        html_content = f""" 
        <!DOCTYPE html> 
        <html lang="en"> 
        <head> 
            <meta charset="UTF-8"> 
            <meta name="viewport" content="width=device-width, initial-scale=1.0"> 
            <title>Decision Tree Visualization</title> 
            <style> 
                body {{ font-family: Arial, sans-serif; text-align: center}}
                .tree-container {{ display: flex; justify-content: center; margin-top: 20px}}
                .node {{ border: 2px solid blue; border-radius: 100%; background-color: lightcyan; color: lightcyan ;display: inline-block;  position: absolute; width: 24px; height: 24px}} 
                .square {{ position: absolute; transform: translate(14px, 14px); z-index: -1}}
                .leaf {{ display: inline-block; font-weight: bolder; position: absolute; line-height: 24px}}
                .d_r:after {{ content: ""; width: 100%; height: 100%; position: absolute; top: 0; left: 0; background: linear-gradient(to top right, transparent calc(50% - 1px), blue, transparent calc(50% + 1px))}}
                .d_l:after {{ content: ""; width: 100%; height: 100%; position: absolute; top: 0; left: 0; background: linear-gradient(to top left, transparent calc(50% - 1px), blue, transparent calc(50% + 1px))}}
                .tooltip {{ position: absolute; background-color: rgba(0, 0, 0, 0.8); color: white; padding: 8px; border-radius: 5px; white-space: nowrap; visibility: hidden; opacity: 0; transition: opacity 0.3s; font-size: 14px; width: fit-content}}
                .t_l {{ transform: translate(35px, -5px)}}
                .t_r {{ transform: translateX(-100%) translate(-8px, -5px)}}
                .leaf_value{{ display: inline-block; font-weight: bolder; position: absolute; line-height: 24px; transform: translate(6px, 30px)}}
            </style> 
        </head> 
        <body> 
            <h1>""" + title + f"""</h1> 
            <div class="tree-container" id="tree"></div> 

            <div id="tooltip" class="tooltip"></div>
            
            <script> 
                const treeData = {tree_JSON}; 
                const plotData = {plot_JSON};
                
                function iter(node, top, left, d, h) {{
                    if(node.isLeaf == 1){{
                        let nodeElement = document.createElement("div")
                        nodeElement.classList.add("leaf")
                        nodeElement.style.left = left + "%"
                        nodeElement.style.top = top + "%"

                        nodeElement.onmouseover = function(event) {{ showTooltip(event, node, top, left); }};
                        nodeElement.onmouseout = hideTooltip;
                        
                        let nodeValue = document.createElement("div")
                        nodeValue.classList.add("leaf_value")
                        nodeValue.style.left = left + "%"
                        nodeValue.style.top = top + "%"
                        nodeValue.innerText = node.value

                        let canvas = document.createElement("canvas");
                        canvas.width = 36;
                        canvas.height = 36;
                        canvas.style.position = "absolute";
                        canvas.style.top = "-6px";
                        canvas.style.left = "-6px";
                        nodeElement.appendChild(canvas);

                        tree = document.getElementById("tree")
                        tree.appendChild(nodeElement)
                        tree.appendChild(nodeValue)

                        // Disegna il grafico a torta direttamente
                        if (node.distribution) {{
                            drawPieChart(canvas, node.distArray, node.labArray);
                        }}
                        return
                    }}

                    let nodeElement = document.createElement("div")
                    nodeElement.classList.add("node")
                    nodeElement.style.left = left + "%"
                    nodeElement.style.top = top + "%"
                    nodeElement.onmouseover = function(event) {{ showTooltip(event, node, top, left); }};
                    nodeElement.onmouseout = hideTooltip;

                    let lbranchElement = document.createElement("div")
                    lbranchElement.classList.add("square")
                    lbranchElement.classList.add("d_l")
                    lbranchElement.style.right = (100-left) + "%"
                    lbranchElement.style.top = top + "%"
                    lbranchElement.style.height = h + "%"
                    lbranchElement.style.width = d/2 + "%"

                    let rbranchElement = document.createElement("div")
                    rbranchElement.classList.add("square")
                    rbranchElement.classList.add("d_r")
                    rbranchElement.style.left = left + "%"
                    rbranchElement.style.top = top + "%"
                    rbranchElement.style.height = h + "%"
                    rbranchElement.style.width = d/2 + "%"

                    // Creare un elemento canvas per il grafico a torta
                    let canvas = document.createElement("canvas");
                    canvas.width = 36;
                    canvas.height = 36;
                    canvas.style.position = "absolute";
                    canvas.style.top = "-6px";
                    canvas.style.left = "-6px";
                    nodeElement.appendChild(canvas);

                    tree = document.getElementById("tree")
                    tree.appendChild(nodeElement)
                    tree.appendChild(lbranchElement)
                    tree.appendChild(rbranchElement)

                    if (node.distribution) {{
                        drawPieChart(canvas, node.distArray, node.labArray);
                    }}

                    iter(node.children[0], top+h, left-d/2, d/2, h)
                    iter(node.children[1], top+h, left+d/2, d/2, h)
                    return
                }}

                // Funzione per disegnare un grafico a torta con il canvas
                function drawPieChart(canvas, distArray, labArray) {{
                    let ctx = canvas.getContext("2d");
                    let total = distArray.reduce((sum, val) => sum + val, 0);
                    let startAngle = 0;

                    distArray.forEach((value, index) => {{
                        let sliceAngle = (value / total) * 2 * Math.PI;

                        ctx.beginPath();
                        ctx.moveTo(18, 18);  // Centro del cerchio
                        ctx.arc(18, 18, 16, startAngle, startAngle + sliceAngle);
                        ctx.closePath();
                        ctx.fillStyle = plotData.colors[labArray[index] % plotData.colors.length]; // Usa i colori in loop
                        ctx.fill();

                        startAngle += sliceAngle;
                    }});

                    // Disegna il bordo del cerchio
                    ctx.beginPath();
                    ctx.arc(18, 18, 17, 0, 2 * Math.PI);
                    ctx.strokeStyle = "blue";
                    ctx.lineWidth = 2;
                    ctx.stroke();
                }}

                function showTooltip(event, node, top, left) {{
                    tooltip = document.getElementById('tooltip')

                    if(node.isLeaf == 1){{
                        tooltip.innerText = "Id: " + node.position + "\\nN: " + node.labels + "\\nClass distribution: [" + node.distribution +  "]\\nGCR: [" + node.gcr + "]\\nImpurty: " + node.impurity.toFixed(3) + "\\nPrediction: "+ node.value;
                    }} else {{
                        if(node.lift1 && node.lift2) {{
                            tooltip.innerText = "Id: " + node.position + "\\nN: " + node.labels +  "\\nClass distribution: [" + node.distribution + "]\\nLIFT left: [" + node.lift1 + "]\\nLIFT right: [" + node.lift2 + "]\\nFeature: " + node.feature + "\\nThreshold left: [" + node.treshold + "]\\nGpi: " + node.gpi.toFixed(3) + "\\nPpi: " + node.ppi.toFixed(3)  + "\\nImpurty: " + node.impurity.toFixed(3);
                        }} else {{
                            tooltip.innerText = "Id: " + node.position + "\\nN: " + node.labels +  "\\nClass distribution: [" + node.distribution + "]\\nFeature: " + node.feature + "\\nThreshold left: [" + node.treshold + "]\\nGpi: " + node.gpi.toFixed(3) + "\\nPpi: " + node.ppi.toFixed(3)  + "\\nImpurty: " + node.impurity.toFixed(3);
                        }}
                    }}

                    if(left > 50){{
                        tooltip.style.left = left + "%"
                        tooltip.style.top = top + "%"
                        tooltip.classList.add("t_r")
                    }} else {{
                        tooltip.style.left = left + "%"
                        tooltip.style.top = top + "%"
                        tooltip.classList.add("t_l")
                    }}

                    tooltip.style.visibility = 'visible'
                    tooltip.style.opacity = 1
                }}

                function hideTooltip() {{
                    tooltip = document.getElementById('tooltip')
                    tooltip.classList.remove("t_r")
                    tooltip.classList.remove("t_l")
                    tooltip.style.visibility = 'hidden'
                    tooltip.style.opacity = 0
                }}

                d = 80/(plotData.l_c+plotData.r_c);
                h = 65/(plotData.depth);

                iter(treeData, 15, 10+d*plotData.l_c, d, h)
            </script> 
        </body>
        </html> 
        """

        with open(output_file, "w", encoding="utf-8") as f:
              f.write(html_content)
  
class Categorizer:
    def __init__(self, min_classes=1, max_classes=10, min_variance=0.3, min_clustersize=1, model="kmeans", k=None, k_method="elbow", sequential=False):
        self.min_classes = min_classes
        self.max_classes = max_classes
        self.min_variance = min_variance
        self.min_clustersize = min_clustersize
        self.k = k
        self.k_method = k_method
        self.sequential = sequential

        if model != "kmeans":
            model = "kmeans"

        if k_method!="elbow" and k_method!="silhouette": 
            return
        
        self.base_model = getattr(C_SCRIPT, model + "_c")
        self.model = getattr(C_SCRIPT, model + "_" + k_method + "_c")
    
    def categorize(self, x):
        if self.k_method!="elbow" and self.k_method!="silhouette":
            print("An unknown method for selecting the number K of clusters has been given. Please choose between 'elbow' and 'silhouette")
            return
        
        if self.sequential:
            N = len(x)
            x_ord = sorted(x)
            self.base_model.restype = FLO
            
            n_classes = 1
            highest_i = 1
            highest_impurity = 1

            splits = [-1, len(x)-1]
            while(n_classes < self.min_classes or (highest_impurity > self.min_variance and n_classes < self.max_classes)):
                size = len(x_ord[(splits[highest_i-1]+1):splits[highest_i]])

                XFLO = FLO * size
                x_flo_c = XFLO(*x_ord[(splits[highest_i-1]+1):splits[highest_i]])

                LABINT = INT * size
                x_cat_1 = LABINT()

                split = splits[highest_i-1]+ 1 + int(self.base_model(size, 2, x_flo_c, x_cat_1, FLO(3.0)))
                bisect.insort(splits, split)
                print(splits)

                highest_i = 0
                highest_impurity = 0
                n_classes += 1
                for i in range(1, len(splits)):
                    mean = np.mean(x_ord[splits[i-1]+1:splits[i]])
                    impurity = ((x_ord[splits[i-1]+1:splits[i]] - mean) ** 2).sum()/N   
                    print("impurity:", impurity)
                    if(impurity > highest_impurity):
                        highest_impurity = impurity
                        highest_i = i

            splits[0] = 0
            treshold = [x_ord[splits[_]] for _ in range(len(splits))]
            labels = np.arange(len(splits)-1)
            x_cat = pd.cut(x, bins=treshold, labels=labels, include_lowest=True, right=True, duplicates='drop')
            return x_cat
                

        XFLO = FLO * len(x)
        x_flo_c = XFLO(*x)

        LABINT = INT * len(x)
        x_cat_1 = LABINT()

        if(self.k):
            self.base_model(len(x), self.k, x_flo_c, x_cat_1)
            return list(x_cat_1)

        x_cat_2 = LABINT()
        if(self.k_method == "elbow"):
            k_opt = self.model(len(x), self.max_classes, self.min_classes, self.min_clustersize, x_flo_c, x_cat_1, x_cat_2)

            if (k_opt % 2 == 0):
                return list(x_cat_2)
            else:
                return list(x_cat_1)
        
        if(self.k_method == "silhouette"):
            self.model.restype = FLO
            k_opt = self.model(len(x), self.max_classes, self.min_classes, self.min_clustersize, x_flo_c, x_cat_1, x_cat_2)
            frac = math.modf(k_opt)
            if(frac == 0):
                return list(x_cat_1)
            else:
                return list(x_cat_2)
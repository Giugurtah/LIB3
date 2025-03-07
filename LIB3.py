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
FLO = ct.c_float
PFLO = ct.POINTER(FLO)
PPFLO = ct.POINTER(PFLO)


C_SCRIPT = ct.CDLL(str(Path(__file__).parent.absolute()) + '/LIB3C.so')

gpi_c = getattr(C_SCRIPT,'gpi_c')
gpi_c.restype = FLO
extract = getattr(C_SCRIPT,'extract')
extract.restype = INT

class Node:
    def __init__(self, 
                 gpi=None, 
                 ppi=None, 
                 position=None, 
                 impurity=None,
                 feature=None, 
                 treshold=None, 
                 alpha=None,
                 beta=None,
                 left=None, 
                 right=None,
                 LIFT_1=None,
                 LIFT_2=None,
                 LCR=None,
                 distribution=None,
                 N=None,
                 labels=None,
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
                 FAST=False):
        
        col_config = {
            "lba": ['id', 'Node Type', 'Value', 'Splitting Variable', 'n', 'Class Distribution', 'Alpha', 'Beta', 'LIFT_K1', 'LIFT_K2', 'Treshold', 'Impurity', 'Gpi', 'Ppi'],
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
        self.model = getattr(C_SCRIPT, model + "_c")
        self.compound_feats = compound_feats
        self.min_impurity = min_impurity
        self.FAST = FAST

        self.results = pd.DataFrame(columns=col)
        self.depth = None
        self.l_c = None
        self.r_c = None
        self.root=None

    def fit(self, X, y):
        if(self.impurity_method!="entropy" and self.impurity_method!="gini" and self.impurity_method!="error"):
            print("ERROR. An impurity method that is not gini, error or entropy has been given.\nA gini impurity method has been automatichally given")
            self.impurity_method = "gini"
        
        #The recursive function is called
        self.root = self._grow_tree(X, y)

    def predict(self, X):
                arr = []
                count = 1
                for i in range(len(X)):
                    print("Riga ispezionata:", count)
                    arr.append(self._traverse_tree(X.iloc[i], self.root))
                    count = count+1
                
                return np.array(arr)

    #Growing functions
    def _grow_tree(self, X, y, depth=0, pos=1):
        print("INIZIO CALCOLO NODO ALLA PROFONDITA':", depth, " ID:", pos) #TODO da cancellare
        # All columns with constant values are dropped
        nunique = X.nunique()
        col_to_drop = nunique[nunique == 1].index
        X = X.drop(col_to_drop, axis=1)
        compound_feature = ""

        # Generale sizes of X and y
        n_samples, n_feats = X.shape
        n_labels = len(np.unique(y))
        impurity = self._impurity(y)
        distribution =  np.unique(y, return_counts=True)[1]/len(y)
        print("Il sottoinsieme del dataset contiene attualmenete", n_samples, "samples,", n_feats, "feats e", n_labels, "labels") #TODO da cancellare

        # Check the stopping criteria
        if(depth>=self.max_depth or n_labels==1 or n_samples<self.min_samples_split or n_feats==0 or impurity<self.min_impurity):
            leaf_value = y.mode()[0]
            self._update_plot(depth, pos)
            self._update_results(pos, "Leaf Node", leaf_value, None, len(y), distribution , None, None, None, impurity, None, None, None, None)

            print("(criteria1)  Nodo foglia raggiunto. Valore di y associato a questo nodo:", leaf_value) #TODO da cancellare
            print("\n")
            return Node(position=pos, value=leaf_value, impurity=impurity, distribution=distribution, N=len(y), labels=np.unique(y))

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
        best_feature, best_treshold, best_ppi, best_idx, alpha, beta = self._find_best_predictor(X, n_labels, array_of_Fs, gpi_index, gpi)

        # Check the stopping criteria
        if(gpi[best_idx]<self.min_gpi or best_ppi<self.min_ppi or best_ppi == 0):
            leaf_value = y.mode()[0]
            self._update_plot(depth, pos)
            self._update_results(pos, "Leaf Node", leaf_value, None, len(y), distribution, None, None, None, impurity, None, None, None, None)

            print("(criteria2) Nodo foglia raggiunto. Valore di y associato a questo nodo:", leaf_value) #TODO da cancellare
            print("\n")
            return Node(position=pos, value=leaf_value, impurity=impurity, distribution=distribution, N=len(y), labels=np.unique(y))
        
        print("best feature", best_feature)
        print("gpi", gpi[best_idx])
        print("best ppi", best_ppi)
        print(best_treshold)
        print("best treshold", np.unique(X[best_feature])[np.array(best_treshold)[:] > 0])
        print("\n")

        indexL, indexR = self._split(X[best_feature], best_treshold)
        
        # Create the child nodes
        self._update_results(pos, "Parent Node", None, best_feature, len(y), distribution, best_treshold, gpi[best_idx], best_ppi, impurity, alpha, beta, [x[0] for x in beta]/distribution, [x[1] for x in beta]/distribution)
        if(best_feature == compound_feature):
            X = X.loc[:, X.columns != feature_1]
            X = X.loc[:, X.columns != feature_2]
            left = self._grow_tree(X.loc[indexL, X.columns != best_feature], y[indexL], depth+1, 2*pos)
            right = self._grow_tree(X.loc[indexR, X.columns != best_feature], y[indexR], depth+1, 2*pos+1)
        else:
            X = X.loc[:, X.columns != compound_feature]
            left = self._grow_tree(X.loc[indexL, X.columns != best_feature], y[indexL], depth+1, 2*pos)
            right = self._grow_tree(X.loc[indexR, X.columns != best_feature], y[indexR], depth+1, 2*pos+1)

        return Node(gpi=gpi[best_idx], ppi=best_ppi, position=pos, feature=best_feature, 
                    treshold=np.unique(X[best_feature])[np.array(best_treshold)[:] > 0], 
                    left=left, right=right, impurity=impurity, distribution=distribution, 
                    N=len(y), labels=np.unique(y),
                    LIFT_1=[x[0] for x in beta]/distribution, LIFT_2=[x[1] for x in beta]/distribution)

    def _split(self, X_col, best_treshold):
        L = np.unique(X_col)[np.array(best_treshold)[:] > 0]
        indexL = np.isin(X_col, L)
        indexR =  np.isin(X_col, L) == False
        return indexL, indexR

    def _find_best_predictor(self, X, n_labels, array_of_Fs, gpi_i, gpi):
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

            # This function doesn't return anything, it updates currente_treshold, alpha and beta instead
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
    
    #Update functions   
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
    
    def _update_results(self, id, type, value, feature, n, distribution, treshold, gpi, ppi, impurity, a, b, lift1, lift2):
        if(self.model_name == "lba"):
            self.results.loc[len(self.results)] = [id, type, value, feature, n, distribution, a, b, lift1, lift2, treshold, impurity, gpi, ppi]
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
        dist = ", ".join(f"{label}={self._custom_round(prob, 2)}" for label, prob in zip(node.labels.tolist(), node.distribution.tolist()))
        if node._is_leaf_node():
            return {
                    "isLeaf": 1,
                    "distribution" : dist,
                    "position": node.position,
                    "value": str(node.value),
                    "impurity": node.impurity,
                    "labels": node.N,
                    }
        
        if(self.model_name == "lba"):
            lift1 = ", ".join(f"{label}={self._custom_round(prob, 2)}" for label, prob in zip(node.labels.tolist(), node.LIFT_1.tolist()))
            lift2 = ", ".join(f"{label}={self._custom_round(prob, 2)}" for label, prob in zip(node.labels.tolist(), node.LIFT_2.tolist()))
        else:
            lift1 = None
            lift2 = None

        return {
            "isLeaf": 0,
            "feature": node.feature,
            "distribution" : dist,
            "lift1" : lift1,
            "lift2" : lift2,
            "treshold" : ', '.join(str(x) for x in node.treshold),
            "position": node.position,
            "gpi" : node.gpi,
            "ppi" : node.ppi,
            "impurity": node.impurity,
            "labels": node.N,
            "children": [
                self._recurse(node.left),
                self._recurse(node.right),
            ]
        }
        
    def plot_html(self, output_file = "tree_visualization.html"):
        tree_JSON = json.dumps(self._recurse(self.root), indent=4)
        dataPlot = {
            "l_c": self.l_c,
            "r_c": self.r_c,
            "depth": self.depth
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
                .leaf {{ border: 2px solid blue; background-color: lightcyan; display: inline-block; color: blue; font-weight: bolder; position: absolute; width: 24px; height: 24px; line-height: 24px}}
                .d_r:after {{ content: ""; width: 100%; height: 100%; position: absolute; top: 0; left: 0; background: linear-gradient(to top right, transparent calc(50% - 1px), blue, transparent calc(50% + 1px))}}
                .d_l:after {{ content: ""; width: 100%; height: 100%; position: absolute; top: 0; left: 0; background: linear-gradient(to top left, transparent calc(50% - 1px), blue, transparent calc(50% + 1px))}}
                .tooltip {{ position: absolute; background-color: rgba(0, 0, 0, 0.8); color: white; padding: 8px; border-radius: 5px; white-space: nowrap; visibility: hidden; opacity: 0; transition: opacity 0.3s; font-size: 14px; width: fit-content}}
                .t_l {{ transform: translate(30px, -5px)}}
                .t_r {{ transform: translateX(-100%) translate(-2px, -5px)}}
            </style> 
        </head> 
        <body> 
            <h1>Decision Tree Visualization</h1> 
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
                        nodeElement.innerText = node.value

                        nodeElement.onmouseover = function(event) {{ showTooltip(event, node, top, left); }};
                        nodeElement.onmouseout = hideTooltip;
                        
                        tree = document.getElementById("tree")
                        tree.appendChild(nodeElement)
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

                    tree = document.getElementById("tree")
                    tree.appendChild(nodeElement)
                    tree.appendChild(lbranchElement)
                    tree.appendChild(rbranchElement)

                    iter(node.children[0], top+h, left-d/2, d/2, h)
                    iter(node.children[1], top+h, left+d/2, d/2, h)
                    return
                }}

                function showTooltip(event, node, top, left) {{
                    tooltip = document.getElementById('tooltip')

                    if(node.isLeaf == 1){{
                        tooltip.innerText = "Id: " + node.position + "\\nN: " + node.labels + "\\nClass distribution: [" + node.distribution +  "]\\nImpurty: " + node.impurity.toFixed(3) + "\\nValue: "+ node.value;
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
    def __init__(self, min_classes=1, max_classes=10, min_impurity=0.3, impurity_method="gini"):
        self.min_classes = min_classes
        self.max_classes = max_classes
        self.min_impurity = min_impurity
        self.impurity_method = impurity_method

    def categorize(self, x, y):
        x_sort, y_sort = zip(*sorted(zip(x, y), reverse=False))

        y_encoded, classi_uniche = pd.factorize(y_sort)
        mod_y = np.bincount(y_encoded)  

        MODY = INT * len(mod_y)
        mod_y_c = MODY()

        YENC = INT * len(y_encoded)
        y_enc_c = YENC()
        # Array of pointers initialization
        for k in range(len(y_encoded)):
            y_enc_c[k] = y_encoded[k]

        #Arrays
        splits = [0, len(y_encoded)]
        #Param
        n_classes = 1
        highest_i = 1
        while(n_classes < self.min_classes or (highest_impurity > self.min_impurity and n_classes < self.max_classes)):
            print("ITERAZIONE NUMERO: ", n_classes)
            mod_y = np.bincount(y_encoded[splits[highest_i-1]:splits[highest_i]], minlength=len(mod_y)) 
            for k in range(len(mod_y)):
                mod_y_c[k] = mod_y[k]
            
            print("Start:", splits[highest_i-1], "Stop:", splits[highest_i])
            split = extract(len(y_encoded), len(mod_y), y_enc_c, mod_y_c, splits[highest_i-1], splits[highest_i])   
            
            bisect.insort(splits, split)
            n_classes += 1

            print("Split: ", split, "splits: ", splits)
            highest_impurity = 0
            highest_i = 0
            for i in range(1, len(splits)):
                current_impurity = self._impurity(y_encoded[splits[i-1]:splits[i]]) 
                if current_impurity > highest_impurity:
                    highest_i = i
                    highest_impurity = current_impurity

        #Categorization of x
        labels = np.arange(len(splits)-1)
        treshold = np.zeros(len(splits))
        treshold[0] = x_sort[0]
        treshold[len(treshold)-1] = x_sort[len(x)-1]
        for i in range(1, len(splits) - 1):
            treshold[i] = (x_sort[splits[i]] + x_sort[splits[i-1]])/2

        print(treshold)
        x_cat = pd.cut(x, bins=treshold, labels=labels, right=False)
        return x_cat
    
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

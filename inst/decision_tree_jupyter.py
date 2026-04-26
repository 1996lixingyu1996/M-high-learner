from graphviz import Digraph, Source
import numpy as np
from collections import Counter
from IPython.display import display

The datasets only support numpy array
class ConstrainedDecisionTree:
    ## 
    def __init__(self, max_depth=None, min_samples_split=2, min_samples_leaf=1, feature_names=None):
        self.max_depth = max_depth
        self.min_samples_split = min_samples_split
        self.min_samples_leaf = min_samples_leaf
        self.tree = None
        self.feature_names = feature_names

    def fit(self, X, y):
        if hasattr(X, "iloc"):
            X = X.to_numpy()
        self.tree = self._build_tree(X, y, list(range(X.shape[1])), depth=0)

    def _best_split(self, X, y, available_features):
        best_feature, best_threshold, best_gini = None, None, 1.0
        for feature in available_features:
            thresholds = np.unique(X[:, feature])
            for t in thresholds:
                left_mask = X[:, feature] <= t
                right_mask = ~left_mask
                if left_mask.sum() < self.min_samples_leaf or right_mask.sum() < self.min_samples_leaf:
                    continue
                gini = self._gini_index(y[left_mask], y[right_mask])
                if gini < best_gini:
                    best_gini = gini
                    best_feature = feature
                    best_threshold = t
        return best_feature, best_threshold

    def _gini_index(self, y_left, y_right):
        def gini(y):
            counts = Counter(y)
            impurity = 1.0
            for count in counts.values():
                p = count / len(y)
                impurity -= p ** 2
            return impurity
        total = len(y_left) + len(y_right)
        return (len(y_left) / total) * gini(y_left) + (len(y_right) / total) * gini(y_right)

    def _build_tree(self, X, y, available_features, depth):
        if len(set(y)) == 1 or len(available_features) == 0 or \
           (self.max_depth and depth >= self.max_depth) or \
           len(y) < self.min_samples_split:
            return {"label": Counter(y).most_common(1)[0][0], "samples": len(y)}

        feature, threshold = self._best_split(X, y, available_features)
        if feature is None:
            return {"label": Counter(y).most_common(1)[0][0], "samples": len(y)}

        left_mask = X[:, feature] <= threshold
        right_mask = ~left_mask
        if left_mask.sum() < self.min_samples_leaf or right_mask.sum() < self.min_samples_leaf:
            return {"label": Counter(y).most_common(1)[0][0], "samples": len(y)}

        new_available = [f for f in available_features if f != feature]
        return {
            "feature": feature,
            "threshold": threshold,
            "samples": len(y),
            "left": self._build_tree(X[left_mask], y[left_mask], new_available, depth + 1),
            "right": self._build_tree(X[right_mask], y[right_mask], new_available, depth + 1)
        }

    def predict_one(self, x):
        node = self.tree
        while "feature" in node:
            if x[node["feature"]] <= node["threshold"]:
                node = node["left"]
            else:
                node = node["right"]
        return node["label"]

    def predict(self, X):
        if hasattr(X, "iloc"):
            X = X.to_numpy()
        return np.array([self.predict_one(row) for row in X])

    # ===== Visualization =====
    def _add_nodes(self, dot, node, node_id):
        if "feature" in node:
            feat_name = self.feature_names[node["feature"]] if self.feature_names else f"X[{node['feature']}]"
            label = f"{feat_name} ≤ {node['threshold']}\n(samples={node['samples']})"
            dot.node(str(node_id), label)
            left_id = node_id * 2 + 1
            right_id = node_id * 2 + 2
            self._add_nodes(dot, node["left"], left_id)
            self._add_nodes(dot, node["right"], right_id)
            dot.edge(str(node_id), str(left_id), "Yes")
            dot.edge(str(node_id), str(right_id), "No")
        else:
            label = f"Class={node['label']}\n(samples={node['samples']})"
            dot.node(str(node_id), label, shape="box")

    def visualize_notebook(self):
        dot = Digraph()
        self._add_nodes(dot, self.tree, 0)
        src = Source(dot.source)
        display(src)

##### usage ######
tree = ConstrainedDecisionTree(max_depth=3, feature_names=feature_names,min_samples_split=300, min_samples_leaf=200)
tree.fit(X, y)
tree.visualize_notebook() ## for jupyter

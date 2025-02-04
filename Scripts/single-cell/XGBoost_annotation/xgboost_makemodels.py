'''
    Script to create the models using XGBoost on reference T.Paralytica dataset
'''
import scanpy as sc
import xgboost
from xgboost import XGBClassifier
from xgboost import plot_importance
from sklearn.model_selection import train_test_split
from sklearn.metrics import classification_report, confusion_matrix
import matplotlib.pyplot as plt
import os
import pandas as pd
import numpy as np
from sklearn.feature_selection import mutual_info_classif
from sklearn.preprocessing import LabelEncoder
import argparse
import joblib
import seaborn as sns
from scipy.sparse import issparse

def parse_args():
    parser = argparse.ArgumentParser(description="Process directories for single-cell data.")
    
    # Define the arguments for base_dir and output_base_dir
    parser.add_argument('--base_dir', type=str, required=True, 
                        help="Base directory for input data.")
    parser.add_argument('--output_base_dir', type=str, required=True, 
                        help="Base directory for output data.")
    
    # Parse the arguments
    return parser.parse_args()

args = parse_args()

'''
BASE_DIR='/storage/gge/Quique/Cells2SpineData/Pilot/single-cell/full_seq_aggr/indiv_samples'
OUTPUT_BASE_DIR='/home/quiquevb/Cells2Spine/Cells2Spine/Outputs/single-cell/full_seq_aggr/indiv_samples'
'''

base_dir = args.base_dir
parent_dir = os.path.dirname(base_dir)
model_dir = os.path.join(parent_dir, "XGBoost_models")
os.makedirs(model_dir, exist_ok=True)

output_base_dir = args.output_base_dir
parent_output_dir = os.path.dirname(output_base_dir)
output_dir = os.path.join(parent_output_dir, "XGBoost_outputs")
os.makedirs(output_dir, exist_ok=True)

# Path for reference of single cell data
single_cell_ref_h5 = "/storage/gge/Quique/TabulaeParalytica/single/GSE234774.h5"

def train_models(single_cell_ref_h5, condition):

    adata_ref = sc.read(single_cell_ref_h5)
    if condition == "healthy":
        adata_ref_subs = adata_ref[adata_ref.obs['label'].isin(['uninjured'])].copy()
    elif condition == "injured":
        adata_ref_subs = adata_ref[adata_ref.obs['label'].isin(['7d', '14d'])].copy()
    
    # Preprocess adata: Log1p normalize
    sc.pp.normalize_total(adata_ref_subs, target_sum=1e4)
    sc.pp.log1p(adata_ref_subs)

    # Reduce dimensionality

    # Define X as gene expression matrix and y as celltype annotations
    X = adata_ref_subs.X
    y = adata_ref_subs.obs['cell_type']

    # Compute information gain for each gene
    info_gain_scores = mutual_info_classif(X, y, random_state=42)

    # Select the top-k genes with the highest information gain
    k = 500  # Example: top 500 genes
    top_k_genes = np.argsort(info_gain_scores)[-k:]

    # Subset X to keep only the top-k genes
    X_selected = X[:, top_k_genes]

    # Convert genes to pd dataframe and add gene names
    selected_gene_names = adata_ref_subs.var_names[top_k_genes]
    
    if issparse(X_selected):
        X_selected = X_selected.toarray()

    # Convert to pandas DataFrame with appropriate gene names as columns
    X_selected_df = pd.DataFrame(X_selected, columns=selected_gene_names)  # Convert sparse to dense if necessary

    # Plot the information gain scores
    plt.bar(range(len(info_gain_scores)), sorted(info_gain_scores, reverse=True))
    plt.xlabel('Genes (sorted by IG)')
    plt.ylabel('Information Gain Score')
    plt.title('Information Gain for Genes')
    plt.savefig(os.path.join(output_dir, f"Information_Gain_Genes_{condition}.png"))
    plt.close()

    # Initialize the label encoder
    label_encoder = LabelEncoder()

    # Encode cell types into numeric labels
    y_encoded = label_encoder.fit_transform(y)

    # Check the mapping of labels to integers
    print(dict(zip(label_encoder.classes_, label_encoder.transform(label_encoder.classes_))))

    X_train, X_test, y_train, y_test = train_test_split(X_selected_df, y_encoded, test_size=0.2, random_state=42)

    # Initialize the XGBoost classifier
    xgb_model = XGBClassifier(
        objective='multi:softmax',  # For multi-class classification
        num_class=len(np.unique(y_encoded)),    # Number of cell types (classes)
        eval_metric='mlogloss',    # Multiclass log-loss
        use_label_encoder=False,
        random_state=42
    )
    
    # Step 7: Ensure the feature names are passed explicitly to the model
    # Set the feature names from adata_ref.var_names (genes)

    # Train the model
    xgb_model.fit(X_train, y_train)
    # Set feature names in the booster
    xgb_model.get_booster().feature_names = X_train.columns.tolist()

    path_baseline_model = os.path.join(model_dir, f'baseline_xgb_cell_type_model_{condition}.joblib')
    # Save the baseline model before hyperparameter tuning
    joblib.dump(xgb_model, path_baseline_model)
    print("Baseline model saved")

    # Predict on validation set
    y_pred = xgb_model.predict(X_test)

    # Map predicted numeric labels back to original labels
    y_pred_original = label_encoder.inverse_transform(y_pred)

    # Map actual numeric labels back to original labels for evaluation
    y_test_original = label_encoder.inverse_transform(y_test)

    # Generate classification report
    report = classification_report(y_test_original, y_pred_original, target_names=label_encoder.classes_)

    # Save classification report to a text file
    with open(os.path.join(output_dir, f"classification_report_{condition}.txt"), 'w') as f:
        f.write(report)

    # Generate confusion matrix
    cm = confusion_matrix(y_test_original, y_pred_original)

    # Plot confusion matrix as a heatmap
    plt.figure(figsize=(10, 8))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=label_encoder.classes_, yticklabels=label_encoder.classes_)
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.title(f'Confusion Matrix {condition}')
    plt.tight_layout()

    # Save confusion matrix as a PNG file
    plt.savefig(os.path.join(output_dir, f"confusion_matrix_{condition}.png"))
    plt.close()
    
    from sklearn.model_selection import GridSearchCV

    # Define your parameter grid
    param_grid = {
        'max_depth': [3, 6, 9],
        'learning_rate': [0.01, 0.1, 0.3],
        'n_estimators': [50, 100, 200]
    }

    # GridSearchCV to find the best hyperparameters
    grid_search = GridSearchCV(estimator=xgb_model, param_grid=param_grid, cv=3, scoring='accuracy')
    grid_search.fit(X_train, y_train)

    # Retrieve the best model
    best_model = grid_search.best_estimator_
    
    # Set feature names in the booster of the best model
    best_model.get_booster().feature_names = X_train.columns.tolist()
    
    path_best_model = os.path.join(model_dir, f'best_xgb_cell_type_model_{condition}.joblib')
    # Save the best model
    joblib.dump(best_model, path_best_model)
    print("Best model saved")

    # Use the best model for predictions
    y_pred = best_model.predict(X_test)

    # Map predicted numeric labels back to original labels
    y_pred_original = label_encoder.inverse_transform(y_pred)

    # Map actual numeric labels back to original labels for evaluation
    y_test_original = label_encoder.inverse_transform(y_test)

    # Evaluate model
    print(classification_report(y_test_original, y_pred_original))

    # Save the classification report
    report = classification_report(y_test_original, y_pred_original, target_names=label_encoder.classes_)
    with open(os.path.join(output_dir, f"classification_report_{condition}.txt"), 'w') as f:
        f.write(report)

    # Generate and save confusion matrix
    cm = confusion_matrix(y_test_original, y_pred_original)
    plt.figure(figsize=(10, 8))
    sns.heatmap(cm, annot=True, fmt='d', cmap='Blues', xticklabels=label_encoder.classes_, yticklabels=label_encoder.classes_)
    plt.xlabel('Predicted')
    plt.ylabel('True')
    plt.title(f'Confusion Matrix {condition}')
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"confusion_matrix_{condition}.png"))
    plt.close()

    # Get feature importances and plot them
    feature_importances = best_model.feature_importances_  # Use best_model here
    sorted_idx = np.argsort(feature_importances)[::-1]

    # Get top k features
    top_k_features = sorted_idx[:k]
    top_feature_names = [adata_ref_subs.var_names[i] for i in top_k_features]
    top_feature_importances = feature_importances[top_k_features]

    # Plot the top-k feature importances
    plt.figure(figsize=(10, 6))
    plt.barh(top_feature_names, top_feature_importances, color='skyblue')
    plt.xlabel('Feature Importance')
    plt.title(f'Top {k} Features Based on XGBoost Importance')
    plt.gca().invert_yaxis()
    plt.tight_layout()
    plt.savefig(os.path.join(output_dir, f"Feature_Importance_{condition}.png"))
    plt.close()

# Train model for "healthy"
train_models(single_cell_ref_h5, condition='healthy')

# Train model for "injured"
train_models(single_cell_ref_h5, condition='injured')

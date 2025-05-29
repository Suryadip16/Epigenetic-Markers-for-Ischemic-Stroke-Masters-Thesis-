import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import shap
from sklearn.pipeline import Pipeline
from sklearn.linear_model import LogisticRegression
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.model_selection import StratifiedKFold
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import (
    accuracy_score, precision_score, recall_score, f1_score,
    confusion_matrix, roc_auc_score, precision_recall_curve,
    auc, roc_curve
)
from xgboost import XGBClassifier

# Load data
df = pd.read_csv("/home/ibab/PycharmProjects/Ischemic_stroke/ml_data_array_epi_int_with_DE_epi_all.csv")
X = df.drop(columns=["Group", "Unnamed: 0"]).select_dtypes(include=[np.number])
print(X.applymap(lambda x: isinstance(x, (list, np.ndarray))).any())

y = df["Group"].to_numpy()

models = {
    "Logistic Regression": LogisticRegression(penalty='l2', solver='liblinear', max_iter=1000),
    "Random Forest": RandomForestClassifier(n_estimators=100, random_state=42),
    "SVM": SVC(kernel='linear', probability=True, random_state=42),
    "XGBoost": XGBClassifier(use_label_encoder=False, eval_metric='logloss', random_state=42)
}

cv = StratifiedKFold(n_splits=10, shuffle=True, random_state=42)

for name, clf in models.items():
    print(f"\n=== {name} ===")

    accs, precs, recs, specs, f1s, roc_aucs, pr_aucs = [], [], [], [], [], [], []
    all_y, all_proba = [], []

    for train_idx, test_idx in cv.split(X, y):
        X_train, X_test = X.iloc[train_idx], X.iloc[test_idx]
        y_train, y_test = y[train_idx], y[test_idx]

        # pipe = Pipeline([
        #     ('scaler', StandardScaler()),
        #     ('clf', clf)
        # ])
        pipe = Pipeline([
            ('clf', clf)
        ])

        pipe.fit(X_train, y_train)
        y_pred = pipe.predict(X_test)
        y_proba = pipe.predict_proba(X_test)[:, 1]

        all_y.extend(y_test)
        all_proba.extend(y_proba)

        accs.append(accuracy_score(y_test, y_pred))
        precs.append(precision_score(y_test, y_pred))
        recs.append(recall_score(y_test, y_pred))
        tn, fp, fn, tp = confusion_matrix(y_test, y_pred).ravel()
        specs.append(tn / (tn + fp))
        f1s.append(f1_score(y_test, y_pred))
        roc_aucs.append(roc_auc_score(y_test, y_proba))
        prec_vals, rec_vals, _ = precision_recall_curve(y_test, y_proba)
        pr_aucs.append(auc(rec_vals, prec_vals))

    # Summary statistics
    print(f"Accuracy: {np.mean(accs):.3f} ± {np.std(accs):.3f}")
    print(f"Precision: {np.mean(precs):.3f} ± {np.std(precs):.3f}")
    print(f"Recall (Sensitivity): {np.mean(recs):.3f} ± {np.std(recs):.3f}")
    print(f"Specificity: {np.mean(specs):.3f} ± {np.std(specs):.3f}")
    print(f"F1 Score: {np.mean(f1s):.3f} ± {np.std(f1s):.3f}")
    print(f"ROC AUC: {np.mean(roc_aucs):.3f} ± {np.std(roc_aucs):.3f}")
    print(f"PR AUC: {np.mean(pr_aucs):.3f} ± {np.std(pr_aucs):.3f}")

    # Single PR Curve
    prec_vals, rec_vals, _ = precision_recall_curve(all_y, all_proba)
    pr_auc = auc(rec_vals, prec_vals)
    plt.figure(figsize=(7, 5))
    plt.plot(rec_vals, prec_vals, color='blue', label=f'PR AUC = {pr_auc:.2f}')
    plt.xlabel('Recall')
    plt.ylabel('Precision')
    plt.title(f'Precision-Recall Curve ({name})')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # Single ROC Curve
    fpr, tpr, _ = roc_curve(all_y, all_proba)
    roc_auc = roc_auc_score(all_y, all_proba)
    plt.figure(figsize=(7, 5))
    plt.plot(fpr, tpr, color='green', label=f'ROC AUC = {roc_auc:.2f}')
    plt.plot([0, 1], [0, 1], 'k--')
    plt.xlabel('False Positive Rate')
    plt.ylabel('True Positive Rate (Sensitivity)')
    plt.title(f'ROC Curve ({name})')
    plt.legend()
    plt.grid(True)
    plt.tight_layout()
    plt.show()

    # === Feature Importance ===
    print(f"\nTop 10 Important Features for {name}:")

    if name in ["Logistic Regression", "SVM"]:
        coefs = pipe.named_steps['clf'].coef_.flatten()
        feature_names = X.columns
        feature_importance = pd.Series(coefs, index=feature_names).sort_values(key=np.abs, ascending=False)

        top_n = 30
        top_features = feature_importance.head(top_n)

        # Save to CSV
        top_features.to_csv(f"top_{top_n}_features_{name.replace(' ', '_')}.csv", header=["Coefficient"])

        # Bar colors based on sign of coefficient
        colors = top_features.apply(lambda x: 'red' if x > 0 else 'green')

        plt.figure(figsize=(8, 6))
        top_features.plot(kind='barh', color=colors)
        plt.xlabel("Coefficient Value")
        plt.title(f"Top {top_n} Important Features ({name})")
        plt.gca().invert_yaxis()
        plt.grid(True)
        plt.tight_layout()
        plt.show()

    elif name in ["Random Forest", "XGBoost"]:
        pipe.fit(X, y)
        try:
            # SHAP for Tree-Based Models
            explainer = shap.TreeExplainer(pipe.named_steps['clf'])
            shap_values = explainer.shap_values(X)

            # Beeswarm plot for feature importance
            shap.summary_plot(shap_values, X, plot_type="dot", show=False, max_display=30)

            # Mean absolute SHAP values
            shap_importance = pd.DataFrame({
                'Feature': X.columns,
                'Mean_Abs_SHAP': np.abs(shap_values).mean(axis=0)
            }).sort_values(by='Mean_Abs_SHAP', ascending=False)

            # Save top N
            shap_importance.head(30).to_csv(f"top_30_shap_features_{name.replace(' ', '_')}.csv", index=False)

            plt.title(f"SHAP Feature Importance - {name}")
            plt.tight_layout()
            plt.show()
        except Exception as e:
            print(f"SHAP failed for {name}: {e}")
            importance = pd.Series(pipe.named_steps['clf'].feature_importances_, index=X.columns).sort_values(
                ascending=False)
            print(importance.head(10))
            importance.head(10).plot(kind='barh')
            plt.title(f"Top 10 Important Features - {name}")
            plt.gca().invert_yaxis()
            plt.tight_layout()
            plt.show()

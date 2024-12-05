import numpy as np

def save_feature_matrix(feature_matrix, output_file):
    """
    Save the feature matrix to a CSV file.
    """
    np.savetxt(output_file, feature_matrix, delimiter=",")
    print(f"Feature matrix saved to {output_file}")

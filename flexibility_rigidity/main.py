from data_loader import load_datasets
from distance_calculator import calculate_features
from output import save_feature_matrix


def main():
    # Load training data
    train_list, test_list = load_datasets()

    # Calculate feature matrix
    feature_matrix = calculate_features(train_list)

    # Save the feature matrix to CSV
    save_feature_matrix(feature_matrix, "Lor_2kerneltrain_2.5-1-12_2.5-5.5-21.csv")


if __name__ == "__main__":
    main()

import os
import pandas as pd
import numpy as np

def simulate_coral_expression(num_genes=100, num_samples=6):
    # Build an absolute path to data/mock_expression_data.csv
    base_dir = os.path.dirname(os.path.dirname(__file__))      # project root
    data_dir = os.path.join(base_dir, "data")
    os.makedirs(data_dir, exist_ok=True)
    out_path = os.path.join(data_dir, "mock_expression_data.csv")

    # Simulate counts
    genes = [f"Gene_{i:03d}" for i in range(1, num_genes+1)]
    ctrl = np.random.poisson(lam=50, size=(num_genes, num_samples//2))
    stress = ctrl.copy()
    stress[:10] += np.random.poisson(lam=40, size=(10, num_samples//2))
    stress[10:20] -= np.random.poisson(lam=20, size=(10, num_samples//2))
    stress = np.clip(stress, 1, None)

    df = pd.DataFrame(
        np.hstack([ctrl, stress]),
        index=genes,
        columns=[f"Ctrl_{i+1}" for i in range(num_samples//2)]
                + [f"Stress_{i+1}" for i in range(num_samples//2)]
    )

    # Always overwrite, never leave an empty file
    df.to_csv(out_path)
    return out_path

if __name__ == "__main__":
    print("Generating mock data at:", simulate_coral_expression())

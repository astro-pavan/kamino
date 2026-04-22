import torch
import torch.nn as nn

class ChemistryEmulator(nn.Module):
    def __init__(self):
        super().__init__()
        # Shared trunk to learn the general thermodynamic state
        self.trunk = nn.Sequential(
            nn.Linear(11, 256),
            nn.SiLU(),
            nn.Linear(256, 256),
            nn.SiLU(),
            nn.Linear(256, 128),
            nn.SiLU()
        )
        
        # Head 1: Elemental equilibrium concentrations (7 outputs)
        self.b_eq_head = nn.Linear(128, 7)
        
        # Head 2: pH (1 output)
        self.pH_head = nn.Linear(128, 1)

    def forward(self, x):
        features = self.trunk(x)
        b_eq = self.b_eq_head(features)
        pH = self.pH_head(features)
        return b_eq, pH

def custom_loss(pred_b_eq, pred_pH, true_b_eq, true_pH):
    # Weight the pH heavily, as speciation depends entirely on it
    mse = nn.MSELoss()
    return 0.2 * mse(pred_b_eq, true_b_eq) + 0.8 * mse(pred_pH, true_pH)
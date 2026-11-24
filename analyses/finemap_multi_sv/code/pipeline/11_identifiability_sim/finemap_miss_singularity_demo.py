import numpy as np
np.random.seed(0)

n, p = 500, 8
X = np.random.randn(n, p)
X[:, 1] += 0.8 * X[:, 0]; X[:, 2] += 0.6 * X[:, 0]; X[:, 5] += 0.7 * X[:, 4]
Xs = (X - X.mean(0)) / X.std(0)
R = np.corrcoef(Xs, rowvar=False)
print("R is full rank:", np.linalg.matrix_rank(R), "of", p,
      "| min eig(R) = %.4g" % np.linalg.eigvalsh(R).min())

w = np.random.randn(p)

denom = np.sqrt(w @ R @ w)
c = (R @ w) / denom
var_zsv = (w @ R @ w) / (denom**2)
print("Var(z_SV) in reference = %.6f (expected 1)" % var_zsv)

def schur(A_block_var, cov_row, R):
    return A_block_var - cov_row @ np.linalg.solve(R, cov_row)

A1 = np.block([[R, c[:, None]], [c[None, :], np.array([[var_zsv]])]])
print("\n[INFO=1]  Schur complement = %.3e" % schur(var_zsv, c, R))
print("[INFO=1]  min eigenvalue(A_aug) = %.3e   (==0 -> SINGULAR)" % np.linalg.eigvalsh(A1).min())
print("[INFO=1]  rank(A_aug) = %d of %d" % (np.linalg.matrix_rank(A1), p+1))

cv_R2 = 0.30
s = np.sqrt(cv_R2)
c2 = s * c
var2 = s*s * var_zsv
print("\n[INFO=cv_R2=%.2f]  Schur complement = %.3e   (still 0 -> weighting does NOT help)"
      % (cv_R2, schur(var2, c2, R)))
A2 = np.block([[R, c2[:, None]], [c2[None, :], np.array([[var2]])]])
print("[INFO=cv_R2]  min eigenvalue(A_aug) = %.3e" % np.linalg.eigvalsh(A2).min())

tau2 = 0.5
A3_var = var_zsv + tau2
print("\n[genuine residual tau^2=%.2f]  Schur complement = %.3e   (>0 -> COHERENT, non-degenerate)"
      % (tau2, schur(A3_var, c, R)))
A3 = np.block([[R, c[:, None]], [c[None, :], np.array([[A3_var]])]])
print("[genuine residual]  min eigenvalue(A_aug) = %.3e" % np.linalg.eigvalsh(A3).min())
print("\nCONCLUSION: deterministic ImpG z_SV => Schur=0 => singular at ANY INFO weight;")
print("only genuine SNP-orthogonal residual variance (tau^2>0) makes the model coherent.")

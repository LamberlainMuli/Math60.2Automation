# %%
from module4 import *
from scipy.stats import f_oneway,t
import numpy as np
import pandas as pd
import math
def print_hypothesis_anova():
    print_latex("$H_0 : $ All means are equal")
    print_latex("$H_A : $ means are not all equal")


def get_lsd(levels, mse, df_residuals, alpha=0.05,df=None):
    t_critical = t.ppf(1 - alpha/2, df_residuals)
    print(f"t-critical value for alpha={alpha}, df={df_residuals} is {t_critical:.4f}")
    print(f"Mean Square Error = {mse:.4f}")
    if df is None:
        df=pd.DataFrame(levels)
    groups = list(df.columns)
    lsd_results = []
    
    for i in range(len(groups)):
        for j in range(i+1, len(groups)):
            groupA, groupB = groups[i], groups[j]
            nA = len(df[groupA].dropna())
            nB = len(df[groupB].dropna())
            xA = df[groupA].dropna().mean()
            xB = df[groupB].dropna().mean()
            lsd = t_critical *((mse * (1/nA + 1/nB)))**0.5
            print(nA,nB,xA,xB,mse,lsd)
            
            diff_means = abs(xA - xB)
            print()
            print(f"Currently comparing {groupA} vs {groupB}...")
            print()
            print("Using the formula for confidence interval:")            
            print_latex(f"$\\bar{{x}}_A - \\bar{{x}}_B \\pm t_{{\\alpha/2,df}} \\sqrt{{\\frac{{MSE}}{{n_A}} + \\frac{{MSE}}{{n_B}}}}$")
            print("Substituting the values we have:")
            print_latex(f"$\\bar{{x}}_A - \\bar{{x}}_B \\pm {t_critical:.4f} \\sqrt{{\\frac{{{mse:.4f}}}{{{nA}}} + \\frac{{{mse:.4f}}}{{{nB}}}}}$")
            print(f"LSD = {t_critical*math.sqrt(((mse)*((1/nA)+(1/nB))))}")
            print("The confidence intervals are:")
            print_latex(f"Lower: {diff_means:.4f} - {lsd:.4f} = {diff_means - lsd:.4f}")
            print_latex(f"Upper: {diff_means:.4f} + {lsd:.4f} = {diff_means + lsd:.4f}")
            
            conclusion = "significant difference" if lsd < diff_means else "no significant difference"
            
            lsd_results.append([groupA, groupB, lsd, diff_means, conclusion])
    
    lsd_df = pd.DataFrame(lsd_results, columns=['Group A', 'Group B', 'LSD', '|xA-xB|', 'Conclusion'])
    
    print("\nHere are the results of the LSD test for all pairs of groups:")
    print(lsd_df)

def get_anova(levels,alpha=0.05):
    if type(levels)==pd.DataFrame:
        levels = {col: levels[col].dropna().values for col in levels.columns}
    F,p=f_oneway(*levels.values())
    print(f"F={F:.4f}, p={p:.4f}")
    
    if p < alpha:
        print(f"Since, {p:.4f} < {alpha}, we reject the null hypothesis that all means are equal.")
        reject_null = True
    else:
        print(f"Since, {p:.4f} > {alpha}, we fail to reject the null hypothesis that all means are equal.")
        reject_null = False
        
    for level in levels.keys():
        print(f"Participants in the {level} group had a mean of {np.mean(levels[level]):.2f} and a standard deviation of {np.std(levels[level],ddof=1):.2f}.")
    
    # Calculate the Mean Square Error
    residuals = []
    grand_mean = np.mean(np.concatenate(list(levels.values())))
    
    for level in levels.keys():
        residuals.extend(levels[level] - np.mean(levels[level]))
    
    ss_residuals = sum([r**2 for r in residuals])
    df_residuals = len(np.concatenate(list(levels.values()))) - len(levels)
    mse = ss_residuals / df_residuals
    
    if reject_null:
        print("\nSince we have rejected the null hypothesis in the ANOVA test, let's perform the LSD test.")
        get_lsd(levels, mse, df_residuals, alpha)
    else:
        print("\nSince we have failed to reject the null hypothesis in the ANOVA test, there's no need to perform the LSD test.")

    return F, p, mse, df_residuals, reject_null
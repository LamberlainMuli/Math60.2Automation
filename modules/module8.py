# %%
import pandas as pd
import matplotlib.pyplot as plt
import statsmodels.api as sm
from IPython.display import display

def simple_linear_regression(x, y):
    # Add a constant to the independent variable
    X = sm.add_constant(x)

    # Fit the model
    model = sm.OLS(y, X)
    results = model.fit()
    # Get Peearson R
    r = results.rsquared**0.5
    print()
    print(f"The pearson's R is {r}")
    print()
    # Print the summary statistics
    print(results.summary2(float_format="%.4f"))

    # Calculate the residuals
    residuals = results.resid

    # Create a new DataFrame for the results
    results_df = pd.DataFrame({
        'x': x,
        'y': y,
        'a + bx': results.fittedvalues,
        'Residual Error': residuals
    })
    
    print("\nTable with residuals and predicted values:")
    display(results_df)
    return results_df


def print_result(results,x,y,xlabel,ylabel,title):
    
    # Create a scatter plot of x and y and plot the regression line
    plt.figure(figsize=(10, 6))
    plt.scatter(x, y, label='Data')
    plt.plot(x, results['a + bx'], 'r', label='Best fit line')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title('title')
    plt.legend()
    plt.show()

    # Create a residual plot
    plt.figure(figsize=(10, 6))
    plt.scatter(x, results['Residual Error'], label='Residuals')
    plt.axhline(y=0, color='r', linestyle='--', label='y=0')
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    plt.legend()
    plt.show()




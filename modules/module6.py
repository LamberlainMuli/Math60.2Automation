# %%
from scipy.stats import norm
from module4 import *
from module5 import *

def test_hypothesis(z, alternative, alpha,null=""):
    if alternative == 'left':
        # Left-tailed test
        print_latex(f"$P(Z\\le {{z}}) = norm.cdf({z}) = {norm.cdf(z)}$")
        p_value = norm.cdf(z)
    elif alternative == 'right':
        # Right-tailed test
        print_latex(f"$P(Z\\ge {{z}}) = 1 - P(Z\\le {{z}}) = 1-norm.cdf({z}) = {1 - norm.cdf(z)}$")
        p_value = 1 - norm.cdf(z)
    elif alternative == 'two-sided':
        # Two-tailed test
        print_latex(f"$P(Z\\le -{{z}})+ P(Z\\ge {{z}}) = 2*(1-norm.cdf(|{z}|)) = {2*(1 - norm.cdf(abs(z)))}$")
        p_value = 2*(1 - norm.cdf(abs(z)))
    else:
        return "Invalid alternative hypothesis. It should be 'less', 'greater' or 'two-sided'"
    print(f"the test statistic is {z:.4f} and the p-value is {p_value:.4f}.")
    if p_value < alpha:
        decision = "reject the null hypothesis"
        print(f"Since, {p_value} < {alpha}, we reject the null hypothesis that {null}.")
    elif p_value == alpha:
        decision = "reject the null hypothesis"
        print(f"Since, {p_value} = {alpha}, we reject the null hypothesis that {null}.")
    else:
        decision = "fail to reject the null hypothesis"
        print(f"Since, {p_value} > {alpha}, we fail to reject the null hypothesis that {null}.")

    return p_value, decision

def get_test_statistic(x, u, se):
    print("Given the formula for test statistic:")
    print_latex("$z = \\frac{\\bar{x} - \\mu}{SE}$")
    print("We can calculate the test statistic as:")
    print_latex(f"$z = \\frac{{{x:.2f} - {u:.2f}}}{{{se}}} = {(x-u)/(se)}$")
    return (x - u) / se

def print_hypothesis(about="",u=0,alternative="\\neq"):
    if about=="mean":
        print_latex(f"$H_0: \\mu = {u}$")
        print_latex(f"$H_A: \\mu {alternative} {u}$")
    elif about=="two means":
        print_latex(f"$H_0: \\mu_1 - \\mu_2 = {u}$")
        print_latex(f"$H_A: \\mu_1 - \\mu_2 {alternative} {u}$")
    elif about=="proportion":
        print_latex(f"$H_0: p = {u}$")
        print_latex(f"$H_A: p {alternative} {u}$")
    elif about=="two proportions":
        print_latex(f"$H_0: p_1 - p_2 = {u}$")
        print_latex(f"$H_A: p_1 - p_2 {alternative} {u}$")
    else:
        print("Invalid hypothesis")

def get_standard_error_two_proportions(p_1=None, p_2=None, p=None, n_1=None, n_2=None,hypothesis=True):
    print("Given the formula for standard error of difference in proportions:")
    print_latex("$SE = \\sqrt{\\frac{p_1(1-p_1)}{n_1} + \\frac{p_2(1-p_2)}{n_2}}$")
    if hypothesis:
        print("It can be shown that the standard error of difference in proportions under the null hypothesis is:")
        print_latex("$SE = \\sqrt{p(1-p) \\left( \\frac{1}{n_1} + \\frac{1}{n_2} \\right)}$")
        print("where p is the pooled proportion:")
        print_latex("$p = \\frac{p_1n_1 + p_2n_2}{n_1 + n_2}$")
        if p==None:
            p = (p_1 * n_1 + p_2 * n_2) / (n_1 + n_2)
            print_latex(f"$p = \\frac{{{p_1} \\times {n_1} + {p_2} \\times {n_2}}}{{{n_1} + {n_2}}} = {p:.4f}$")
        else:
            print_latex(f"$p = {p:.4f}$")
        se = (p * (1 - p) * (1 / n_1 + 1 / n_2)) ** 0.5
        print_latex(f"$SE = \\sqrt{{{p} \\times (1 - {p}) \\left( \\frac{{1}}{{{n_1}}} + \\frac{{1}}{{{n_2}}} \\right)}} = {se:.4f}$")
        return se
    else:
        p = None
        se = (p_1 * (1 - p_1) / n_1 + p_2 * (1 - p_2) / n_2) ** 0.5
        print("Solving for SE we get:")
        print_latex(f"$SE = \\sqrt{{\\frac{{{p_1} \\times (1 - {p_1})}}{{{n_1}}} + \\frac{{{p_2} \\times (1 - {p_2})}}{{{n_2}}}}} = {se:.4f}$")
        return se

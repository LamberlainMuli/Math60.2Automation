# %% [markdown]
# Module 5.1: Estimators and Point Estimation

# %%
import numpy as np
from module4 import *
def estimator_mean(sample):
    return np.mean(sample)
def estimator_proportion(sample):
    return np.mean(sample)
def estimator_sd(sample):
    return np.std(sample,ddof=1)
def estimate_mean(mean,n,sd):
    se=get_standard_error_mean(sd=sd,n=n)
    print_latex(f"Point estimate of the population mean is: $\\hat \\mu \\approx N({mean},{se})$")
    return mean,se
def estimate_proportion(n,p):
    se=get_standard_error_proportion(n=n,p=p)
    print_latex(f"Point estimate of the population proportion is: $\\hat p \\approx N({p},{se})$")
    return p,se

def get_alpha(alpha=None,confidence=None):
    print("Calculating alpha...")
    if alpha is None and confidence is None:
        print("Please provide a value for alpha or confidence")
    elif alpha is None:
        print_latex(f"From the formula $\\alpha = 1-$ confidence $ \\alpha $ is $ {1-confidence}$")
        return 1-confidence,confidence
    elif confidence is None:
        print_latex(f"From the formula confidence $= 1-\\alpha$ confidence is ${1-alpha}$")
        return alpha,1-alpha
    return None,None


def get_z_alpha(za=None,alpha=None,confidence=None,half=False):
    print("Calculating za...")
    alpha,confidence = get_alpha(alpha=alpha,confidence=confidence)
    if alpha is None and za is None:
        print("Please provide a value for alpha or za")
    elif za is None:
        if half:
            za_answer=norm.ppf(1-alpha/2)
            print_latex(f"From the formula $Z_{{\\alpha/2}} = F_X^-1(alpha/2)$, we can simply substitute $\\alpha$")
            print_latex(f"$Z_{{\\alpha/2}} = F_X^-1({alpha}/2) \\approx {za_answer}$")
        else:
            za_answer=norm.ppf(1-alpha)
            print_latex(f"From the formula $Z_{{\\alpha}} = F_X^-1(alpha)$, we can simply substitute $\\alpha$")
            print_latex(f"$Z_{{\\alpha}} = F_X^-1({alpha}) \\approx {za_answer}$")
        return za_answer,alpha
    elif alpha is None:
        alpha_answer=norm.cdf(za)
        print_latex(f"From the formula $Z_{{\\alpha}} = F_X^-1(alpha)$, we can isolate $\\alpha$ by using the cdf")
        print_latex(f"$\\alpha = F_X(Z_{{\\alpha}})$")
        print_latex(f"$\\alpha = F_X({za}) \\approx {alpha_answer}$")
        return za,alpha_answer
    else:
        return za,alpha
    return None,None

def get_margin_of_error_mean(moe_mean=None,se=None,sd=None,n=None,za=None,alpha=None,confidence=None):
    print("Calculating margin of error...")
    se = get_standard_error_mean(se=se,sd=sd,n=n)
    za,alpha = get_z_alpha(za=za,alpha=alpha,confidence=confidence,half=True)
    if moe_mean is not None and [za,se,moe_mean].count(None)>1:
        return moe_mean,za,se
    elif [za,moe_mean,se].count(None)!=1:
        print("Please provide exactly two of the following: margin of error, se, za")
        return None,None,None
    elif moe_mean is None:
        print_latex("From the formula margin of error $ = Z_{\\alpha/2} \\times SE$, we simply substitute to get:")
        print_latex(f" margin of error $= {za} \\times {se} \\approx {za*se}$")
        return za*se,za,se
    elif se is None:
        print_latex("From the formula $MOE_{mean} = Z_{\\alpha/2} \\times SE$, we can isolate se by dividing both sides by $Z_{\\alpha/2}$")
        print_latex("SE = $\\frac{MOE_{mean}}{Z_{\\alpha/2}}$")
        print_latex(f"SE = $\\frac{{{moe_mean}}}{{{za}}} \\approx {moe_mean/za}$")
        return moe_mean,za,moe_mean/za
    elif za is None:
        print_latex("From the formula $MOE_{mean} = Z_{\\alpha/2} \\times SE$, we can isolate $Z_{\\alpha/2}$ by dividing both sides by $SE$")
        print_latex("$Z_{\\alpha/2} = \\frac{MOE_{mean}}{SE}$")
        print_latex(f"$Z_{{\\alpha/2}} = \\frac{{{moe_mean}}}{{{se}}} \\approx {moe_mean/se}$")
        return moe_mean,moe_mean/se,se
    else:
        return moe_mean,za,se



def get_confidence_interval_mean(lower=None,upper=None,mean=None,moe_mean=None,se=None,sd=None,n=None,za=None,alpha=None,confidence=None):
    print("Calculating confidence interval...")
    moe_mean,za,se=get_margin_of_error_mean(moe_mean=moe_mean,se=se,sd=sd,n=n,za=za,alpha=alpha,confidence=confidence)
    if lower is None and upper is None:
        # If no interval provided, require both mean and MOE
        if mean is None or moe_mean is None:
            print("Please provide both the mean and the margin of error, or a complete interval.")
            return None,None,None,None
        else:
            print_latex(f"To get the lower and upper bounds, we simply add and subtract the margin of error from the mean:")
            print_latex(f"Lower bound: ${mean}-{moe_mean} \\approx {mean-moe_mean}$")
            print_latex(f"Upper bound: ${mean}+{moe_mean} \\approx {mean+moe_mean}$")
            print(f"We are {confidence*100 if confidence else 'unkown'}% confident that the population mean is within ({mean-moe_mean},{mean+moe_mean})")
            return mean,moe_mean,mean-moe_mean,mean+moe_mean
    else:
        # If interval is provided, handle different scenarios
        if mean is None and moe_mean is not None:
            if lower is not None:
                mean = lower + moe_mean
                upper = lower + 2*moe_mean
                print_latex(f"To get the mean, we add the margin of error to the lower bound: ${lower}+{moe_mean} \\approx {mean}$")
            elif upper is not None:
                mean = upper - moe_mean
                lower = upper - 2*moe_mean
                print_latex(f"To get the mean, we subtract the margin of error from the upper bound: ${upper}-{moe_mean} \\approx {mean}$")
            return mean,moe_mean,lower,upper
        elif mean is None and moe_mean is None:
            if lower is not None and upper is not None:
                mean = (lower + upper) / 2
                moe_mean=upper-mean
                print_latex(f"To get the mean, we take the midpoint of the interval: $\\frac{{{upper}+{lower}}}{{2}} \\approx {mean}$")
                return mean,moe_mean,lower,upper
            else:
                print("Please provide: the mean or the margin of error or the remanining bound.")
                return None,None,None,None
        elif moe_mean is None and mean is not None:
            if upper is not None:
                moe_mean = upper - mean
                lower = mean - moe_mean
                print_latex(f"To get the margin of error, we subtract the mean from the upper bound: ${upper}-{mean} \\approx {moe_mean}$")
            elif lower is not None:
                moe_mean = mean - lower
                upper = mean + moe_mean
                print_latex(f"To get the margin of error, we subtract the lower bound from the mean: ${mean}-{lower} \\approx {moe_mean}$")
            return mean,moe_mean,lower,upper
        else:
            print("Please provide: the mean or the margin of error or the remanining bound.")
            return mean,moe_mean,lower,upper
          


# %% [markdown]
# $SE = \sqrt{\frac{\sigma_1^2}{n_1} + \frac{\sigma_2^2}{n_2}}$
# 
# $SE^2 = \frac{\sigma_1^2}{n_1} + \frac{\sigma_2^2}{n_2}$
# 
# $\sigma_1 = SE^2 \cdot n_1 - \frac{\sigma_2^2 \cdot n_1}{n_2}$
# 
# $n_1 = \frac{\sigma_1^2}{SE^2 - \frac{\sigma_2^2}{n_2}}$
# 
# $\sigma_2 = SE^2 \cdot n_2 - \frac{\sigma_1^2 \cdot n_2}{n_1}$
# 
# $n_2 = \frac{\sigma_2^2}{SE^2 - \frac{\sigma_1^2}{n_1}}$
# 

# %% [markdown]
# $MOE = z_{\alpha/2}SE$
# 
# $z_{\alpha/2} = MOE/SE$
# 
# $SE = MOE/z_{\alpha/2}$
# 
# diff=$\bar{x_1}-\bar{x_2}$
# 
# diff = upper-MOE
# diff = lower+MOE
# 
# MOE = upper-diff
# MOE = lower+diff
# 
# upper = diff + MOE
# lower = diff - MOE
# 
# 
# 
# 
# 

# %%

def get_standard_error_two_means(se=None,sd_1=None,sd_2=None,n_1=None,n_2=None):
    if se is not None and [se,sd_1,sd_2,n_1,n_2].count(None)>1:
        return se,sd_1,sd_2,n_1,n_2
    if [se,sd_1,sd_2,n_1,n_2].count(None)!=1:
        print("Please provide exactly four of the following: se, sd_1, sd_2, n_1, n_2")
        return se,sd_1,sd_2,n_1,n_2
    elif se==None:
        se_answer=math.sqrt(sd_1**2/n_1+sd_2**2/n_2)
        print_latex(f"We can substitute to the formula $SE = \\sqrt{{\\frac{{SD_1^2}}{{n_1}}+\\frac{{SD_2^2}}{{n_2}}}}$")
        print_latex(f"$SE = \\sqrt{{\\frac{{{sd_1}^2}}{{{n_1}}}+\\frac{{{sd_2}^2}}{{{n_2}}}}} \\approx {se_answer}$")
        return se_answer,sd_1,sd_2,n_1,n_2
    elif sd_1==None:
        sd_1_answer=np.sqrt(se**2 * n_1 - (sd_2**2 * n_1) / n_2)
        print_latex(f"We can substitute to the formula $SD_1 = \\sqrt{{SE^2 * n_1 - \\frac{{SD_2^2 * n_1}}{{n_2}}}}$")
        print_latex(f"$SD_1 = \\sqrt{{{se}^2 * {n_1} - \\frac{{{sd_2}^2 * {n_1}}}{{{n_2}}}}} \\approx {sd_1_answer}$")
        return se,sd_1_answer,sd_2,n_1,n_2

    elif n_1==None:
        n_1_answer=(sd_1**2) / (se**2 - (sd_2**2 / n_2))
        print_latex(f"We can substitute to the formula $n_1 = \\frac{{SD_1}}{{SE^2 - \\frac{{SD_2^2}}{{n_2}}}}$")
        print_latex(f"$n_1 = \\frac{{{sd_1}^2}}{{{se}^2 - \\frac{{{sd_2}^2}}{{{n_2}}}}} \\approx {n_1_answer}$")
        return se,sd_1,sd_2,n_1_answer,n_2

    elif sd_2==None:
        sd_2_answer=np.sqrt(se**2 * n_2 - (sd_1**2 * n_2) / n_1)
        print_latex(f"We can substitute to the formula $SD_2 = \\sqrt{{SE^2 * n_2 - \\frac{{SD_1^2 * n_2}}{{n_1}}}}$")
        print_latex(f"$SD_2 = \\sqrt{{{se}^2 * {n_2} - \\frac{{{sd_1}^2 * {n_2}}}{{{n_1}}}}} \\approx {sd_2_answer}$")
        return se,sd_1,sd_2_answer,n_1,n_2
    
    elif n_2==None:
        n_2_answer=(sd_2**2) / (se**2 - (sd_1**2 / n_1))
        print_latex(f"We can substitute to the formula $n_2 = \\frac{{SD_2}}{{SE^2 - \\frac{{SD_1^2}}{{n_1}}}}$")
        print_latex(f"$n_2 = \\frac{{{sd_2}^2}}{{{se}^2 - \\frac{{{sd_1}^2}}{{{n_1}}}}} \\approx {n_2_answer}$")
        return se,sd_1,sd_2,n_1,n_2_answer




def get_diff(diff=None, x_1=None, x_2=None):
    if [diff, x_1, x_2].count(None) != 1:
        print("Please provide exactly two of the following: diff, x_1, x_2")
        return diff, x_1, x_2
    elif diff is None:
        diff = x_1 - x_2
        print_latex(f"$Diff = x_1 - x_2 = {x_1} - {x_2} = {diff}$")
    elif x_1 is None:
        x_1 = diff + x_2
        print_latex(f"$x_1 = Diff + x_2 = {diff} + {x_2} = {x_1}$")
    elif x_2 is None:
        x_2 = x_1 - diff
        print_latex(f"$x_2 = x_1 - Diff = {x_1} - {diff} = {x_2}$")
    return diff, x_1, x_2



def get_margin_of_error_two_means(moe_two_means=None,se=None,sd_1=None,sd_2=None,n_1=None,n_2=None, za=None,alpha=None,confidence=None):
    print("Calculating margin of error...")
    se,sd_1,sd_2,n_1,n_2 = get_standard_error_two_means(se=se, sd_1=sd_1, sd_2=sd_2, n_1=n_1, n_2=n_2)
    za,alpha = get_z_alpha(za=za,alpha=alpha,confidence=confidence,half=True)
    if moe_two_means is not None and [za, se, moe_two_means].count(None) > 1:
        return moe_two_means,za,se
    if [za, moe_two_means, se].count(None) != 1:
        print("Please provide exactly two of the following: margin of error, se, za")
        return moe_two_means,za,se
    elif moe_two_means is None:
        print_latex("From the formula margin of error $ = Z_{\\alpha/2} \\times SE$, we simply substitute to get:")
        print_latex(f" margin of error $= {za} \\times {se} \\approx {za*se}$")
        return za*se,za,se
    elif se is None:
        print_latex("From the formula $MOE = Z_{\\alpha/2} \\times SE$, we can isolate se by dividing both sides by $Z_{\\alpha/2}$")
        print_latex("SE = $\\frac{MOE}{Z_{\\alpha/2}}$")
        print_latex(f"SE = $\\frac{{{moe_two_means}}}{{{za}}} \\approx {moe_two_means/za}$")
        return moe_two_means,za,moe_two_means/za
    elif za is None:
        print_latex("From the formula $MOE = Z_{\\alpha/2} \\times SE$, we can isolate $Z_{\\alpha/2}$ by dividing both sides by $SE$")
        print_latex("$Z_{\\alpha/2} = \\frac{MOE}{SE}$")
        print_latex(f"$Z_{{\\alpha/2}} = \\frac{{{moe_two_means}}}{{{se}}} \\approx {moe_two_means/se}$")
        return moe_two_means,moe_two_means/se,se
    else:
        return moe_two_means,moe_two_means/se,se



def get_confidence_interval_difference_means(
    diff=None,
        mean_1=None,
        mean_2=None,
    lower=None,
    upper=None,
    moe_two_means=None,
        se=None,
            sd_1=None,
            sd_2=None,
            n_1=None,
            n_2=None,
        za=None,
            confidence=None,
            alpha=None
    ):
    moe_two_means,za,se=get_margin_of_error_two_means(
        moe_two_means=moe_two_means,
        se=se,sd_1=sd_1,sd_2=sd_2,n_1=n_1,n_2=n_2,
        za=za,confidence=confidence,alpha=alpha)
    
    diff, x_1, x_2=get_diff(diff=diff, x_1=mean_1, x_2=mean_2)
    print(diff,moe_two_means,lower,upper)
    if lower is None and upper is None:
        if diff is None or moe_two_means is None:
            print("Please provide both the difference of the means and the margin of error, or a complete interval.")
            return diff, moe_two_means, lower, upper
        else:
            lower = diff - moe_two_means
            upper = diff + moe_two_means
            print_latex(f"$Lower = Diff - MOE = {diff} - {moe_two_means} = {lower}$")
            print_latex(f"$Upper = Diff + MOE = {diff} + {moe_two_means} = {upper}$")
            print(f"We are {confidence*100 if confidence else 'unknown'}% confident that the difference in population means is within ({lower},{upper})")
            return diff, moe_two_means, lower, upper
    else:
        if diff is None and moe_two_means is not None:
            if lower is not None:
                diff = lower + moe_two_means
                print_latex(f"$Diff = Lower + MOE = {lower} + {moe_two_means} = {diff}$")
                if upper is None:
                    upper = lower + 2*moe_two_means
                    print_latex(f"upper = Lower + 2*MOE = {lower} + 2*{moe_two_means} = {upper}$")
            elif upper is not None:
                diff = upper - moe_two_means
                print_latex(f"$Diff = Upper - MOE = {upper} - {moe_two_means} = {diff}$")
                if lower is None:
                    lower = upper - 2*moe_two_means
                    print_latex(f"lower = Upper - 2*MOE = {upper} - 2*{moe_two_means} = {lower}$")
            return diff, moe_two_means, lower, upper
        elif diff is None and moe_two_means is None:
            if lower is not None and upper is not None:
                diff = (lower + upper) / 2
                moe_two_means = upper - diff
                return diff, moe_two_means, lower, upper
            else:
                print("Please provide: the difference or the margin of error or the remaining bound.")
                return diff, moe_two_means, lower, upper
        elif moe_two_means is None and diff is not None:
            if upper is not None:
                moe_two_means = upper - diff
                print_latex(f"$MOE = Upper - Diff = {upper} - {diff} = {moe_two_means}$")
                if lower is None:
                    lower = diff - moe_two_means
                    print_latex(f"lower = Diff - MOE = {diff} - {moe_two_means} = {lower}$")
            elif lower is not None:
                moe_two_means = diff - lower
                print_latex(f"$MOE = Diff - Lower = {diff} - {lower} = {moe_two_means}$")
                if upper is None:
                    upper = diff + moe_two_means
                    print_latex(f"upper = Diff + MOE = {diff} + {moe_two_means} = {upper}$")
            return diff, moe_two_means, lower, upper
        else:
            print("Please provide: the difference or the margin of error or the remaining bound.")
            return diff, moe_two_means, lower, upper
        return diff, moe_two_means, lower, upper
    



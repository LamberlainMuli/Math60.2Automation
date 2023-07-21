import math
from scipy.stats import norm
from IPython.display import display, Latex

def print_latex(string):
    display(Latex(string))


def solution_standardization_population(lower=None,upper=None,probability=None):
    if lower is None and upper is not None:
        print_latex(f"$P(X \\le {{{upper}}}) = {probability}$")
        print_latex(f"$P(\\frac{{X - \\mu}}{{\\sigma}} \\le \\frac{{{{{upper}}}-\\mu}}{{\\sigma}})  = {probability}$")
        print_latex(f"Let $Z=\\frac{{X - \\mu}}{{\\sigma}}$, then $Z\\approx N(0,1)$. Substituting Z, we get")
        print_latex(f"$P(Z \\le \\frac{{{{{upper}}}-\\mu}}{{\\sigma}})  = {probability}$")
        print("Taking the inverse of the cumulative distribution function, we get:")
        print_latex(f"$\\frac{{{{{upper}}}-\\mu}}{{\\sigma}} = $ norm.ppf(${probability}$)")
    elif lower is not None and upper is None:    
        print_latex(f"$P(X \\ge {{{lower}}}) = {probability}$")
        print_latex(f"$P(X \\le {{{lower}}}) = 1-{probability}$")
        print_latex(f"$P(\\frac{{X - \\mu}}{{\\sigma}} \\le \\frac{{{{{lower}}}-\\mu}}{{\\sigma}})  = 1 - {probability}$")
        print_latex(f"Let $Z=\\frac{{X - \\mu}}{{\\sigma}}$, then $Z\\approx N(0,1)$. Substituting Z, we get")
        print_latex(f"$P(Z \\le \\frac{{{{{lower}}}-\\mu}}{{\\sigma}})  = 1 - {probability}$")
        print("Taking the inverse of the cumulative distribution function, we get:")
        print_latex(f"$\\frac{{{{{lower}}}-\\mu}}{{\\sigma}} = $ norm.ppf($1 - {probability}$)")

def probability_population(mean=None,sd=None,lower=None,upper=None,probability=None):
    if mean is None:
        if probability is not None and sd is not None and (lower is not None or upper is not None):
            if lower is not None and upper is not None:
                print(f"Cannot calculate sample size with both lower and upper bounds")
            elif lower is not None:
                mean_answer=lower + norm.ppf(probability, loc=0, scale=1) * sd 
                solution_standardization_population(lower=lower,probability=probability)
                print("Solving for the mean, we get:")
                print_latex(f"$\\mu= {lower} - $ norm.ppf( $1-{probability} $) $* {sd} \\approx {mean_answer}$")
                print(f"The mean is {mean_answer}")
                return mean_answer
            elif upper is not None:
                mean_answer=upper - norm.ppf(probability) * sd
                solution_standardization_population(upper=upper,probability=probability)
                print("Solving for the mean, we get:")
                print_latex(f"$\\mu= {upper} - $ norm.ppf( $ {probability} $) $* {sd} \\approx {mean_answer}$")
                print(f"The mean is {mean_answer}")
                return mean_answer
        else:
            print("To calculate the mean, please provide the probability, standard deviation, and upper or lower bound.")
    elif sd is None:
        if probability is not None and mean is not None and (lower is not None or upper is not None):
            if lower is not None and upper is not None:
                print(f"Cannot calculate sample size with both lower and upper bounds")
            elif lower is not None:
                sd_answer=abs(abs(lower - mean)/norm.ppf(1-probability, loc=0, scale=1) )
                solution_standardization_population(lower=lower,probability=probability)
                print("Solving for the standard deviation, we get:")
                print_latex(f"$ \\sigma= \\frac{{ {lower} - {mean} }} {{ norm.ppf(1- {probability}) }} \\approx {sd_answer} $")
                print(f"The standard deviation is {sd_answer}")
                return sd_answer
            elif upper is not None:
                sd_answer=abs(abs(upper - mean)/norm.ppf(probability, loc=0, scale=1) )
                solution_standardization_population(upper=upper,probability=probability)
                print("Solving for the standard deviation, we get:")
                print_latex(f"$ \\sigma= \\frac{{{upper} - {mean}}}{{ norm.ppf({probability}) }}  \\approx {sd_answer} $")
                print(f"The standard deviation is {sd_answer}")
                return sd_answer
        else:
            print("To calculate the standard deviation, please provide the probability, mean, and upper or lower bound.")
    elif probability is None:
        if mean is not None and sd is not None and (lower is not None or upper is not None):
            if lower is None and upper is None:
                print(f"To calculate the probability, provide lower bound or upper bound.")
            elif lower is not None and upper is None:
                probability_answer=1-norm.cdf(lower,loc=mean,scale=sd)
                print_latex(f"$P(X \\ge {{{lower}}}) = 1-P(X \le {{{lower}}})$")
                print_latex(f"$P(X \\ge {{{lower}}}) = 1- $ norm.cdf(${lower}$,${mean}$,${sd}$) ")
                print_latex(f"$P(X \\ge {{{lower}}}) = {probability_answer}$")
                print(f"The probability that the X is greater than {lower} is {probability_answer}")
                return probability_answer
            elif lower is None and upper is not None:
                probability_answer=norm.cdf(upper,loc=mean,scale=sd)
                print_latex(f"$P(X \\le {{{upper}}}) = $ norm.cdf(${upper}$,${mean}$,${sd}$) ")
                print_latex(f"$P(X \\le {{{upper}}}) = {probability_answer}$")
                print(f"The probability that the X is less than {upper} is {probability_answer}")
                return probability_answer
            elif lower is not None and upper is not None:
                probability_answer=norm.cdf(upper,loc=mean,scale=sd)-norm.cdf(lower,loc=mean,scale=sd)
                print_latex(f"$P({{{lower}}} \\le X \\le {{{upper}}}) = P(X \\le {{{upper}}}) - P(X \\le {{{lower}}})$")
                print_latex(f"$P({{{lower}}} \\le X \\le {{{upper}}}) = $ norm.cdf(${upper}$,${mean}$,${sd}$) - norm.cdf(${lower}$,${mean}$,${sd}$) ")
                print_latex(f"$P({{{lower}}} \\le X \\le {{{upper}}}) = {probability_answer}$")
                print(f"The probability that the X is between {lower} and {upper} is {probability_answer}")
                return probability_answer
        else:
            print("To calculate the probability, please provide the mean, standard deviation, sample size, and lower bound or upper bound.")
    elif lower is None and upper is not None:
        if probability is not None and mean is not None and sd is not None:
            lower_answer=mean - norm.ppf(probability, loc=0, scale=1) * sd 
            solution_standardization_population(lower="lower",probability=probability)
            print("Solving for the lower bound we get:")
            print_latex(f"$lower={{{mean}}}-norm.ppf({probability})*{sd} \\approx {lower_answer}$")
            print(f"The lower bound is {lower_answer}")
            return lower_answer
        else:
            print("To calculate the lower bound, please provide the probability, mean, standard deviation, and sample size.")
    elif lower is not None and upper is None:
        if probability is not None and mean is not None and sd is not None:
            upper_answer=mean + norm.ppf(probability, loc=0, scale=1) * sd
            solution_standardization_population(upper="upper",probability=probability)
            print("Solving for the upper bound we get:")
            print_latex(f"$upper={{{mean}}}+norm.ppf({probability})*{sd} \\approx {upper_answer}$")
            print(f"The upper bound is {upper_answer}")
            return upper_answer
        else:
            print("To calculate the upper bound, please provide the probability, mean, standard deviation, and sample size.")
    else:
        print("Please provide any 3 quantities from (probability, mean, standard deviation, sample size) and either upper or lower bound.")  

def get_standard_error_mean(se=None,sd=None,n=None):
    if se is not None and [se,sd,n].count(None)>1:
        return se
    elif [se,sd,n].count(None)!=1:
        print(f"Please provide two of the three values: se, sd, n")
    elif se is None and sd is not None and n is not None:
        se_answer=sd/math.sqrt(n)
        print_latex(f"Given the formula for the standard error of the mean: $SE=\\frac{{\\sigma}}{{\\sqrt{{n}}}}$. We substitute the values as follows:")
        print_latex(f"$SE=\\frac{{{sd}}}{{\\sqrt{{{n}}}}} \\approx {se_answer}$ ")
        return se_answer
    elif sd is None and n is not None and se is not None:
        sd_answer=se*math.sqrt(n)
        print_latex(f"Given the formula for the standard error of the mean: $SE=\\frac{{\\sigma}}{{\\sqrt{{n}}}}$. We can rearrange as follows to isolate the standard deviation:")
        print_latex(f"$\\sigma=SE\\sqrt n$")
        print("Substituting the values we get:")
        print_latex(f"$\\sigma={se}\\sqrt{{{n}}} \\approx {sd_answer}$")
        return sd_answer
    elif n is None and se is not None and sd is not None:
        n_answer=math.pow(sd/se,2)
        print_latex(f"Given the formula for the standard error of the mean: $SE=\\frac{{\\sigma}}{{\\sqrt{{n}}}}$. We can rearrange as follows to isolate the sample size:")
        print_latex(f"$n=(\\frac{{\\sigma}}{{SE}})^2$")
        print("Substituting the values we get:")
        print_latex(f"$n=(\\frac{{{sd}}}{{{se}}})^2 \\approx {n_answer}$")
        return n_answer
        

def solution_standardization_mean(lower=None,upper=None,probability=None):
    if lower is None and upper is not None:
        print_latex(f"$P(\\bar X \\le {{{upper}}}) = {probability}$")
        print_latex(f"$P(\\frac{{\\bar X - \\mu}}{{\\sigma}} \\le \\frac{{{{{upper}}}-\\mu}}{{\\sigma}})  = {probability}$")
        print_latex(f"Let $Z=\\frac{{\\bar X - \\mu}}{{\\sigma}}$, then $Z\\approx N(0,1)$. Substituting Z, we get")
        print_latex(f"$P(Z \\le \\frac{{{{{upper}}}-\\mu}}{{\\sigma}})  = {probability}$")
        print("Taking the inverse of the cumulative distribution function, we get:")
        print_latex(f"$\\frac{{{{{upper}}}-\\mu}}{{\\sigma}} = $ norm.ppf(${probability}$)")
    elif lower is not None and upper is None:    
        print_latex(f"$P(\\bar X \\ge {{{lower}}}) = {probability}$")
        print_latex(f"$P(\\bar X \\le {{{lower}}}) = 1-{probability}$")
        print_latex(f"$P(\\frac{{\\bar X - \\mu}}{{\\sigma}} \\le \\frac{{{{{lower}}}-\\mu}}{{\\sigma}})  = 1 - {probability}$")
        print_latex(f"Let $Z=\\frac{{\\bar X - \\mu}}{{\\sigma}}$, then $Z\\approx N(0,1)$. Substituting Z, we get")
        print_latex(f"$P(Z \\le \\frac{{{{{lower}}}-\\mu}}{{\\sigma}})  = 1 - {probability}$")
        print("Taking the inverse of the cumulative distribution function, we get:")
        print_latex(f"$\\frac{{{{{lower}}}-\\mu}}{{\\sigma}} = $ norm.ppf($1 - {probability}$)")
                    
def probability_mean(probability=None,mean=None,sd=None,n=None,lower=None,upper=None):
    if mean is None:
        if probability is not None and sd is not None and n is not None and (lower is not None or upper is not None):
            if lower is None:
                mean_answer=upper-norm.ppf(probability,loc=0,scale=1)*sd/math.sqrt(n) 
                solution_standardization_mean(upper=upper,probability=probability)
                print("Solving for the mean, we get:")
                print_latex(f"$\\mu={{{upper}}}$ - norm.ppf(${probability}$)$\\times{{{sd}}} \\approx {mean_answer}$")
                print(f"The mean is {mean_answer}")
                return mean_answer
            elif upper is None:
                mean_answer=lower-norm.ppf(1-probability,loc=0,scale=1)*sd/math.sqrt(n)
                solution_standardization_mean(lower=lower,probability=probability)
                print("Solving for the mean, we get:")
                print_latex(f"$\\mu={{{lower}}}$ - norm.ppf(1 - ${probability}$)$\\times{{{sd}}} \\approx {mean_answer}$")
                print(f"The mean is {mean_answer}")
                return mean_answer
            elif lower is not None and upper is not None:
                print(f"The mean cannot be calculated with both lower and upper bounds")
        else:
            print(f"To calculate the mean, please provide the probability, standard deviation, sample size, and either upper or lower bound")
    elif sd is None:
        if probability is not None and mean is not None and n is not None and (lower is not None or upper is not None):
            if lower is None and upper is not None:
                sd_answer=abs((upper - mean)*math.sqrt(n) / (norm.ppf(probability))) 
                solution_standardization_mean(upper=upper,probability=probability)
                print("Solving for the standard deviation (se), we get:")
                print_latex(f"$\\frac{{sd}}{{\\sqrt{{n}}}}=\\frac{{{{{upper}}}-\\mu}}{{norm.ppf({probability})}} \\approx {abs((upper - mean)/ (norm.ppf(probability))) }$")
                print("Solving for the standard deviation of the sample, we get:")
                print_latex(f"$sd=|\\frac{{\\sqrt{{{n}}}\cdot({{{upper}}}-{{{mean}}})}}{{norm.ppf({probability})}}| \\approx {sd_answer} $")
                print(f"The standard deviation is {sd_answer}")
                return sd_answer
            elif upper is None:
                sd_answer=abs((mean - lower)*math.sqrt(n)/ (norm.ppf(1-probability)))
                solution_standardization_mean(lower=lower,probability=probability)
                print("Solving for the standard deviation (se), we get:")
                print_latex(f"$\\frac{{sd}}{{\\sqrt{{n}}}}=\\frac{{{{{lower}}}-\\mu}}{{norm.ppf(1-{probability})}}  \\approx {abs((lower - mean)/ (norm.ppf(probability))) }$")
                print("Solving for the standard deviation of the sample, we get:")
                print_latex(f"$sd=|\\frac{{\\sqrt{{{n}}}\cdot({{{lower}}}-{{{mean}}})}}{{norm.ppf(1-{probability})}}| \\approx {sd_answer} $")
                print(f"The standard deviation is {sd_answer}")
                return sd_answer
            elif lower is not None and upper is not None:
                print(f"Cannot calculate standard deviation with both lower and upper bounds")
        else:
            print(f"To calculate the standard deviation, please provide the probability, mean, sample size, and either upper or lower bound")
    elif n is None:
        if probability is not None and mean is not None and sd is not None and (lower is not None or upper is not None):
            if lower is None:
                n_answer=((sd*norm.ppf(probability))/(upper - mean))**2 
                solution_standardization_mean(upper=upper,probability=probability)
                print("Solving for the standard deviation (se), we get:")
                print_latex(f"$\\frac{{sd}}{{\\sqrt{{n}}}}=\\frac{{{{{upper}}}-\\mu}}{{norm.ppf({probability})}}  \\approx {abs((upper - mean)/ (norm.ppf(probability))) }$")
                print("Solving for the sample size, we get:")
                print_latex(f"n = $[\\frac{{{{{sd}}}*norm.ppf({probability})}}{{({{{upper}}}-{{{mean}}})}}]^2 \\approx {n_answer}$")
                print(f"The sample size is {n_answer}")
                return n_answer
            elif upper is None:
                n_answer=((sd*norm.ppf(1-probability))/(mean - lower))**2 
                solution_standardization_mean(lower=lower,probability=probability)
                print("Solving for the standard deviation (se), we get:")
                print_latex(f"$\\frac{{sd}}{{\\sqrt{{n}}}}=\\frac{{{{{lower}}}-{{{mean}}}}}{{norm.ppf(1-{probability})}} \\approx {abs((lower - mean)/ (norm.ppf(probability))) }$")
                print("Solving for the sample size, we get:")
                print_latex(f"n = $[\\frac{{{{{sd}}}*norm.ppf(1-{probability})}}{{({{{lower}}}-{{{mean}}})}}]^2 \\approx {n_answer}$")              
                print(f"The sample size is {n_answer}")
                return n_answer
            elif lower is not None and upper is not None:
                print(f"Cannot calculate sample size with both lower and upper bounds")
        else:
            print(f"To calculate the sample size, please provide the probability, mean, standard deviation, and either upper or lower bound")    
    elif probability is None:
        if mean is not None and sd is not None and n is not None:
            if lower is None and upper is None:
                print(f"To calculate the probability, lower bound or upper bound.")
            elif lower is not None and upper is None:
                probability_answer=1-norm.cdf(lower,loc=mean,scale=sd/math.sqrt(n)) 
                print_latex(f"$P(\\bar X \\ge {{{lower}}}) = 1-P(\\bar X \le {{{lower}}})$")
                print_latex(f"$P(\\bar X \\ge {{{lower}}}) = 1- $ norm.cdf(${lower}$,${mean}$,$\\frac{{{sd}}}{{\\sqrt{{{n}}}}} \\approx {sd/math.sqrt(n)}$) ")
                print_latex(f"$P(\\bar X \\ge {{{lower}}}) = {probability_answer}$")
                print(f"The probability that the sample mean is greater than {lower} is {probability_answer}")
                return probability_answer
            elif lower is None and upper is not None:
                probability_answer=norm.cdf(upper,loc=mean,scale=sd/math.sqrt(n)) 
                print_latex(f"$P(\\bar X \\le {{{upper}}}) = $ norm.cdf(${upper}$,${mean}$,$\\frac{{{sd}}}{{\\sqrt{{{n}}}}} \\approx {sd/math.sqrt(n)}$) ")
                print_latex(f"$P(\\bar X \\le {{{upper}}}) = {probability_answer}$")
                print(f"The probability that the sample mean is less than {upper} is {probability_answer}")
                return probability_answer
            elif lower is not None and upper is not None:
                probability_answer=norm.cdf(upper,loc=mean,scale=sd/math.sqrt(n))-norm.cdf(lower,loc=mean,scale=sd/math.sqrt(n)) 
                print_latex(f"$P({{{lower}}} \\le \\bar X \\le {{{upper}}}) = P(\\bar X \\le {{{upper}}}) - P(\\bar X \\le {{{lower}}})$")
                print_latex(f"$P({{{lower}}} \\le \\bar X \\le {{{upper}}}) = $ norm.cdf(${upper}$,${mean}$,$\\frac{{{sd}}}{{\\sqrt{{{n}}}}} \\approx {sd/math.sqrt(n)}$) - norm.cdf(${lower}$,${mean}$,$\\frac{{{sd}}}{{\\sqrt{{{n}}}}} \\approx {sd/math.sqrt(n)}$) ")
                print_latex(f"$P({{{lower}}} \\le \\bar X \\le {{{upper}}}) = {probability_answer}$")
                print(f"The probability that the sample mean is between {lower} and {upper} is {probability_answer}")
                return probability_answer
        else:
            print("To calculate the probability, please provide the mean, standard deviation, sample size, and lower bound or upper bound.")
            
    elif lower is None and upper is not None:
        if probability is not None and mean is not None and sd is not None and n is not None:
            lower_answer=mean + norm.ppf(1 - probability, loc=0, scale=1) * sd/math.sqrt(n) 
            solution_standardization_mean(lower="lower",probability=probability)
            print("Solving for the lower bound we get:")
            print_latex(f"$lower={{{mean}}}-norm.ppf({probability})*\\frac{{{sd}}}{{\\sqrt{{{n}}}}} \\approx {lower_answer}$")
            print(f"The lower bound is {lower_answer}")
            return lower_answer
        else:
            print("To calculate the lower bound, please provide the probability, mean, standard deviation, and sample size.")
    elif lower is not None and upper is None:
        if probability is not None and mean is not None and sd is not None and n is not None:
            upper_answer=mean + norm.ppf(probability, loc=0, scale=1) * sd/math.sqrt(n) 
            solution_standardization_mean(upper="upper",probability=probability)
            print("Solving for the upper bound we get:")
            print_latex(f"$upper={{{mean}}}+norm.ppf({probability})*\\frac{{{sd}}}{{\\sqrt{{{n}}}}} \\approx {upper_answer}$")
            print(f"The upper bound is {upper_answer}")
            return upper_answer
        else:
            print("To calculate the upper bound, please provide the probability, mean, standard deviation, and sample size.")
    else:
        print("Please provide any 3 quantities from (probability, mean, standard deviation, sample size) and either upper or lower bound.")

def get_standard_error_proportion(se=None, p=None, n=None):
    if se is None and p is not None and n is not None:
        se_answer=math.sqrt((p*(1-p))/(n))
        print_latex(f"Given the formula for the standard error of the p: $SE=\\sqrt{{\\frac{{ p(1-p) }}{{ n }} }}$. We substitute the values as follows:")
        print_latex(f"$SE=\\sqrt{{\\frac{{ {p}(1-{p}) }}{{ {n} }} }} \\approx {se_answer}$ ")
        return se_answer
    elif p is None and n is not None and se is not None:
        p_answer_1=(1+math.sqrt(1-4*(math.sqrt(n)*(se**2))))/(2)
        p_answer_2=(1-math.sqrt(1-4*(math.sqrt(n)*(se**2))))/(2)
        print_latex(f"Given the formula for the standard error of the p: $SE=\\sqrt{{\\frac{{ p(1-p) }}{{ n }} }}$. We can rearrange as follows to isolate the p using the quadratic formula:")
        print_latex(f"$p = \\frac{{1 \\pm \\sqrt{{1-4\\left(SE^2 \\cdot \\sqrt{{n}}\\right)}}}}{{2}}$")
        print("Substituting the values we get:")
        print_latex(f"$p = \\frac{{1 \\pm \\sqrt{{1-4\\left({se}^2 \\cdot \\sqrt{n}\\right)}}}}{{2}} \\approx {p_answer_1} $ or ${p_answer_2} $")
        return p_answer_1, p_answer_2
    elif n is None and se is not None and p is not None:
        n_answer=((p-p**2)/(se**2))**2
        print_latex(f"Given the formula for the standard error of the p: $SE=\\sqrt{{\\frac{{ p(1-p) }}{{ n }} }}$. We can rearrange as follows to isolate the sample size:")
        print_latex(f"$ n = \\left(\\frac{{p-p^2}}{{SE^2}}\\right)^2$")
        print("Substituting the values we get:")
        print_latex(f"$ n = \\left(\\frac{{{p}-{p}^2}}{{{se}^2}}\\right)^2 \\approx {n_answer}$")
        return n_answer
    else:
        print(f"Please provide two of the three values: se, p, n")

def solution_standardization_proportion(lower=None, upper=None, probability=None):
    if lower is None and upper is not None:
        print_latex(f"$P(\hat p \\leq {upper}) = {probability}$")
        print_latex(f"$P( \\frac{{{{\hat p}} - p}}{{\sqrt{{p*(1-p)/n}}}} \\leq \\frac{{{{upper}}-p}}{{\sqrt{{p*(1-p)/n}}}}) = {probability}$")
        print_latex(f"Let $Z = \\frac{{{{\hat p}} - p}}{{\sqrt{{p*(1-p)/n}}}}$, then $Z \\sim N(0, 1)$")
        print_latex(f"$P(Z \\leq \\frac{{{{upper}}-p}}{{\sqrt{{p*(1-p)/n}}}}) = {probability}$")
        print_latex(f"Taking inverse of the cumulative distribution function, we get:")
        print_latex(f"$\\frac{{{{upper}} - p}}{{\sqrt{{p*(1-p)/n}}}} = $ norm.ppf(${probability}$)")
    elif lower is not None and upper is None:    
        print_latex(f"$P(\hat p \\geq {lower}) = {probability}$")
        print_latex(f"$P(\hat p \\leq {lower}) = 1 - {probability}$")
        print_latex(f"$P( \\frac{{{{\hat p}} - p}}{{\sqrt{{p*(1-p)/n}}}} \\leq \\frac{{{{lower}} - p}}{{\sqrt{{p*(1-p)/n}}}}) = 1 - {probability}$")
        print_latex(f"Let $Z = \\frac{{{{\hat p}} - p}}{{\sqrt{{p*(1-p)/n}}}}$, then $Z \\sim N(0, 1)$")
        print_latex(f"$P(Z \\leq \\frac{{{{lower}} - p}}{{\sqrt{{p*(1-p)/n}}}}) = 1 - {probability}$")
        print_latex(f"Taking inverse of the cumulative distribution function, we get:")
        print_latex(f"$\\frac{{{{lower}} - p}}{{\sqrt{{p*(1-p)/n}}}} = $ norm.ppf($1 - {probability}$)")

def probability_proportion(probability=None,proportion=None,n=None,lower=None,upper=None):
    if proportion is None:
        if probability is not None and n is not None and (lower is not None or upper is not None):
            if lower is None:
                solution_standardization_proportion(upper=upper,probability=probability)
                print("Solving for p, we get:")
                d=norm.ppf(probability)
                d_latex=f"norm.ppf({probability})"
                a=n+d**2
                a_latex=f"{n}+{d_latex}^2"
                b=-d**2-2*upper*n
                b_latex=f"-{d_latex}^2-2*({upper}*{n})"
                c=upper**2*n
                c_latex=f"{upper}^2*{n}"
                p1=(-b+math.sqrt(b**2-4*a*c))/(2*a)
                p2=(-b-math.sqrt(b**2-4*a*c))/(2*a)
                print_latex(f"$p = \\frac{{-({b_latex}) + \\sqrt{{({b_latex})^2-4({a_latex})({c_latex})}} }}{{2\\cdot ({a_latex})}} \\approx {p1}$")
                print(f"The proportion is {p1}")
                return p1
            elif upper is None:
                solution_standardization_proportion(lower=lower,probability=probability)
                print("Solving for p, we get:")
                d=norm.ppf(1-probability)
                d_latex=f"norm.ppf(1-{probability})"
                a=n+d**2
                a_latex=f"{n}+{d_latex}^2"
                b=-d**2-2*lower*n
                b_latex=f"-{d_latex}^2-2*({lower}*{n})"
                c=lower**2*n
                c_latex=f"{lower}^2*{n}"
                p1=(-b+math.sqrt(b**2-4*a*c))/(2*a)
                p2=(-b-math.sqrt(b**2-4*a*c))/(2*a)
                print_latex(f"$p = \\frac{{-({b_latex}) - \\sqrt{{({b_latex})^2-4({a_latex})({c_latex})}} }}{{2\\cdot ({a_latex})}} \\approx {p2}$")
                print(f"The proportion is {p2}")
                return p2
            elif lower is not None and upper is not None:
                print(f"The proportion cannot be calculated with both lower and upper bounds")
        else:
            print(f"To calculate the proportion, please provide the probability, sample size, and either upper or lower bound")
    elif n is None:
        if probability is not None and proportion is not None and (lower is not None or upper is not None):
            if lower is None:
                solution_standardization_proportion(upper=upper,probability=probability)
                print("Solving for n, we get:")
                d=norm.ppf(probability)
                d_latex=f"norm.ppf({probability})"
                n_answer=((proportion*(1-proportion))/((upper - proportion)/norm.ppf(probability))**2) 
                print_latex(f"$n = \\frac{{({d_latex})^2(p)(1-p)}}{{({upper}-p)^2}} \\approx {(d**2*proportion*(1-proportion))/((upper-proportion)**2)}$")
                print(f"The sample size is {((proportion*(1-proportion))/((upper - proportion)/norm.ppf(probability))**2) :.0f}")
                return n_answer
            elif upper is None:
                solution_standardization_proportion(lower=lower,probability=probability)
                print("Solving for n, we get:")
                d=norm.ppf(1-probability)
                d_latex=f"norm.ppf(1-{probability})"
                n_answer=((proportion*(1-proportion))/((lower - proportion)/norm.ppf(probability))**2)
                print_latex(f"$n = \\frac{{({d_latex})^2(p)(1-p)}}{{({lower}-p)^2}} \\approx {(d**2*proportion*(1-proportion))/((lower-proportion)**2)}$")
                print(f"The sample size is {((proportion*(1-proportion))/((lower - proportion)/norm.ppf(probability))**2) :.0f}")
                return n_answer
            elif lower is not None and upper is not None:
                print(f"Cannot calculate sample size with both lower and upper bounds")
        else:
            print(f"To calculate the sample size, please provide the probability, proportion, and either upper or lower bound")
    elif probability is None:
        if proportion is not None and n is not None:
            if lower is None and upper is None:
                print(f"To calculate the probability, provide lower bound or upper bound.")
            elif lower is not None and upper is None:
                probability_answer=1-norm.cdf(lower,loc=proportion,scale=math.sqrt((proportion*(1-proportion))/n)) 
                print_latex(f"$P(\\hat p \\ge {lower}) = 1 - P(\\hat p \\le {lower})$")
                print_latex(f"$P(\\hat p \\ge {lower}) = 1 - $ norm.cdf($ {lower} $ ,loc=${proportion}$,scale=$\\sqrt{{ \\frac{{ {proportion} \\cdot (1-{proportion})}} {{{n}}} }} \\approx {math.sqrt((proportion*(1-proportion))/n)}$)")
                print_latex(f"$P(\\hat p \\ge {lower}) = {probability_answer}$")
                print(f"The probability that the sample proportion is greater than {lower} is {probability_answer}")
                return probability_answer
            elif lower is None and upper is not None:
                probability_answer=norm.cdf(upper,loc=proportion,scale=math.sqrt((proportion*(1-proportion))/n)) 
                print_latex(f"$P(\\hat p \\le {upper}) = $ norm.cdf($ {upper} $ ,loc=${proportion}$,scale=$\\sqrt{{ \\frac{{ {proportion} \\cdot (1-{proportion})}} {{{n}}} }} \\approx {math.sqrt((proportion*(1-proportion))/n)}$)")
                print_latex(f"$P(\\hat p \\le {upper}) = {probability_answer}$")
                print(f"The probability that the sample proportion is less than {upper} is {probability_answer}")
                return probability_answer
            elif lower is not None and upper is not None:
                probability_answer=norm.cdf(upper,loc=proportion,scale=math.sqrt((proportion*(1-proportion))/n))-norm.cdf(lower,loc=proportion,scale=math.sqrt((proportion*(1-proportion))/n)) 
                print_latex(f"$P({lower} \le \\hat p \le {upper}) = P(\\hat p \le {upper} ) - P(\\hat p \le {lower}) $")
                print_latex(f"$P({lower} \le \\hat p \le {upper}) = $ norm.cdf ($ {upper} $ ,loc=${proportion}$,scale=$\\sqrt{{ \\frac{{ {proportion} \\cdot (1-{proportion})}} {{{n}}} }} \\approx {math.sqrt((proportion*(1-proportion))/n)}$) - norm.cdf($ {lower} $ ,loc=${proportion}$,scale=$\\sqrt{{ \\frac{{ {proportion} \\cdot (1-{proportion})}} {{{n}}} }} \\approx {math.sqrt((proportion*(1-proportion))/n)}$)")
                print_latex(f"$P(\\hat p \\le {upper}) = {probability_answer}$")
                print(f"The probability that the sample proportion is between {lower} and {upper} is {probability_answer}")
                return probability_answer
        else:
            print("To calculate the probability, please provide the proportion, sample size, and lower bound or upper bound.")
    elif lower is None and upper is not None:
        if probability is not None and proportion is not None and n is not None:
            lower_answer=proportion + norm.ppf(1 - probability, loc=0, scale=1) * math.sqrt((proportion*(1-proportion))/n) 
            solution_standardization_proportion(lower="lower",probability=probability)
            print("Solving for the lower bound, we get:")
            print_latex(f"$lower= {proportion} + $ norm.ppf( $ 1 - {probability} $) $* \\sqrt{{ \\frac{{ {proportion} \\cdot (1-{proportion})}} {{{n}}} }} \\approx {lower_answer}$")
            print(f"The lower bound is {lower_answer}")
            return lower_answer
        else:
            print("To calculate the lower bound, please provide the probability, proportion, and sample size.")
    elif lower is not None and upper is None:
        if probability is not None and proportion is not None and n is not None:
            upper_answer=proportion + norm.ppf(probability, loc=0, scale=1) * math.sqrt((proportion*(1-proportion))/n) 
            solution_standardization_proportion(upper="upper",probability=probability)
            print("Solving for the upper bound, we get:")
            print_latex(f"$upper= {proportion} + $ norm.ppf( $ {probability} $) $* \\sqrt{{ \\frac{{ {proportion} \\cdot (1-{proportion})}} {{{n}}} }} \\approx {upper_answer}$")
            print(f"The upper bound is {upper_answer}")
        else:
            print("To calculate the upper bound, please provide the probability, proportion, and sample size.")
    else:
        print("Please provide any 3 quantities from (probability, proportion, sample size) and either upper or lower bound.")


def test_clt_mean(sample_size=30):
    normal=input("Is the sampled population normal? (y/n)")
    if normal=='y':
        print_latex("Since the sampled population is normally distributed, the sampling distribution of $\\bar{X}$ will also be normal, no matter what sample size you choose")
    else:
        symmetric = input("Is the sampled population approximately symmetric? (y/n)")
        if symmetric=='y':
            print_latex("Since the sampled population is approximately symmetric, the sampling distribution of $\\bar{X}$ becomes approximately normal for relatively small values of $n$")
        else:
            print_latex("the sample size $n$ must be larger, with $n$ at least $30$ before the sampling distribution of $\\bar{X}$ becomes approximately normal.")
            if sample_size>=30:
                print_latex("In this case since $n \ge 30$, the sampling distribution of $\\bar{X}$ becomes approximately normal.")
            else:
                print_latex("In this case since $n < 30$, the sampling distribution of $\\bar{X}$ does not become approximately normal.")

def test_clt_proportion(n,p):
    normal=input("Is the sampled population normal? (y/n)")
    if normal=='y':
        print_latex(f"Since the sampled population is normally distributed, the sampling distribution of $\\hat P \\approx N({p},{get_standard_error_proportion(n,p)})$, no matter what sample size you choose")
    else:
        symmetric = input("Is the sampled population approximately symmetric? (y/n)")
        if symmetric=='y':
            print_latex(f"Since the sampled population is approximately symmetric, the sampling distribution of $\\hat P \\approx N({p},{get_standard_error_proportion(n,p)})$, for relatively small values of $n$")
        else:
            print_latex("The sample size $n$ must be larger, with $np > 5$ and $n(1 - p) > 5$ before the sampling distribution of $\\hat{P}$ becomes approximately normal.")
            if n*p>5 and n*(1 - p) > 5:
                print_latex(f"In this case since $np = {n*p :.2f} > 5$ and $n(1 - p) = {n*(1-p) :.2f} > 5$, the sampling distribution of $\\hat P \\approx N({p},{get_standard_error_proportion(p,n)})$")
            else:
                print_latex(f"In this case since $np = {n*p :.2f} \\le 5$ or $n(1-p) = {n*(1-p) :.2f} \\le 5$, the sampling distribution of  $\\hat P \\not\\approx N({p},{get_standard_error_proportion(p,n)})$")

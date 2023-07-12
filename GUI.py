import tkinter as tk
from tkinter import messagebox


import math
def get_standard_error_mean(sd,n):
    return sd/math.sqrt(n)

from scipy.stats import norm
def get_probability_greater_than_mean(greater_than, mean, sd, n):
    return 1-norm.cdf(greater_than,loc=mean,scale=get_standard_error_mean(sd,n))
def get_probability_less_than_mean(less_than, mean, sd, n):
    return norm.cdf(less_than,loc=mean,scale=get_standard_error_mean(sd,n))
def get_probability_between_mean(lower_bound, upper_bound, mean, sd, n):
    return norm.cdf(upper_bound,loc=mean,scale=get_standard_error_mean(sd,n))-norm.cdf(lower_bound,loc=mean,scale=get_standard_error_mean(sd,n))

print(get_probability_less_than_mean(less_than=7,mean=8,sd=4,n=30))
print(get_probability_greater_than_mean(greater_than=7,mean=8,sd=4,n=30))
print(get_probability_between_mean(lower_bound=7,upper_bound=9,mean=8,sd=4,n=30))

from IPython.display import display, Latex

def print_latex(string):
    display(Latex(string))

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

def get_sampling_distribution_mean_variable(mean,sd,n):
    print_latex(f"Sampling distribution is: $N({mean},{get_standard_error_mean(sd,n)})$")


def get_standard_error_proportion(n,p):
    return math.sqrt(p*(1-p)/n)

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

from scipy.stats import norm

def get_probability_greater_than_proportion(greater_than,p,n):
    return 1-norm.cdf(greater_than,loc=p,scale=get_standard_error_proportion(n,p))

def get_probability_less_than_proportion(less_than,p,n):
    return norm.cdf(less_than,loc=p,scale=get_standard_error_proportion(n,p))

def get_probability_between_proportion(lower_bound,upper_bound,p,n):
    return get_probability_less_than_proportion(upper_bound,p,n)-get_probability_less_than_proportion(lower_bound,p,n)

print(get_probability_greater_than_proportion(greater_than=0.4,p=0.3,n=64))
print(get_probability_less_than_proportion(less_than=25/100,p=0.3,n=100))
print(get_probability_less_than_proportion(less_than=100/120,p=0.7,n=120))


def get_sample_standard_deviation(sample):
    return math.sqrt(sum([(xi-sum(sample)/len(sample))**2 for xi in sample])/(len(sample)-1))

def get_point_estimate_mean(sample=[],n=1,total=0):
    if sample==[]:
        return total/n
    else:
        return sum(sample)/len(sample)
    
def estimate_mean(mean,n,sd):
    print_latex(f"Point estimate of the population mean is: $\\hat \\mu \\approx N({mean},{get_standard_error_mean(sd,n)})$")
    return mean,get_standard_error_mean(sd,n)

def compare_point_estimates_mean(mean_1,se_1,mean_2,se_2):
    if se_1<=se_2:
        print_latex(f"Since the standard error of the first point estimate is smaller than the second point estimate, the first point estimate is more precise")
    else:
        print_latex(f"Since the standard error of the first point estimate is larger than the second point estimate, the second point estimate is more precise")

def get_point_estimate_proportion(n,p):
    return p

def estimate_proportion(n,p):
    print_latex(f"Point estimate of the population proportion is: $\\hat p \\approx N({p},{get_standard_error_proportion(n,p)})$")
    return p,get_standard_error_proportion(n,p)

def get_alpha_from_confidence(confidence):
    return 1-confidence
def get_confidence_from_alpha(alpha):
    return 1-alpha
def get_z_alpha(confidence=None,alpha=None):
    if alpha==None:
        alpha=get_alpha_from_confidence(confidence)
    return norm.ppf(1-alpha)
def get_margin_of_error_mean(sd=None,n=None,confidence=None,alpha=None,z_alpha=None,se=None):
    if z_alpha==None:
        if alpha==None:
            alpha=get_alpha_from_confidence(confidence)
        z_alpha=get_z_alpha(alpha=alpha/2)
    if se==None:
        se=get_standard_error_mean(sd=sd,n=n)
    return z_alpha*se

def get_confidence_interval_mean(confidence,n,mean,sd,alpha=None):
    if alpha==None:
        alpha=get_alpha_from_confidence(confidence)
    z_alpha=get_z_alpha(alpha=alpha/2)
    se=get_standard_error_mean(sd=sd,n=n)
    margin_of_error=get_margin_of_error_mean(z_alpha=z_alpha,se=se)
    return f"We are {confidence*100}% confident that the population mean is within ({mean-margin_of_error},{mean+margin_of_error})"

def get_standard_error_two_means(sd_1,sd_2,n_1,n_2):
    return math.sqrt(sd_1**2/n_1+sd_2**2/n_2)
def get_confidence_interval_difference_means(mean_1,mean_2,sd_1,sd_2,n_1,n_2,confidence=None,alpha=None,z_alpha=None,se=None):
    if z_alpha==None:
        if alpha==None:
            alpha=get_alpha_from_confidence(confidence)
        z_alpha=get_z_alpha(alpha=alpha/2)
    
    if se==None:
        se=get_standard_error_two_means(sd_1=sd_1,sd_2=sd_2,n_1=n_1,n_2=n_2)
    
    lower_bound=mean_1-mean_2-z_alpha*se
    upper_bound=mean_1-mean_2+z_alpha*se
    print_latex(f"We can be ${100*(1-alpha)}$% confident that the true difference between the population means within $({lower_bound},{upper_bound})$")
    if lower_bound<=0<=upper_bound:
        return f"Since the interval ({lower_bound},{upper_bound}) contains 0, we can say with ${100*(1-alpha)}$% confidence that there is no significant difference between the two population means"
    elif upper_bound<0:
        return f"Since ${lower_bound} \\le {upper_bound} < 0$ We can be ${100*(1-alpha)}$% confident that $\\mu_1$ is significantly less than $\\mu_2$"
    elif lower_bound>0:
        return f"Since $0<{lower_bound} \\le {upper_bound}$ We can be ${100*(1-alpha)}$% confident that $\\mu_1$ is significantly greater than $\\mu_2%"
    
def get_margin_of_error_proportion(p=None,n=None,confidence=None,alpha=None,z_alpha=None,se=None):
    if z_alpha==None:
        if alpha==None:
            alpha=get_alpha_from_confidence(confidence)
        z_alpha=get_z_alpha(alpha=alpha/2)
    if se==None:
        se=get_standard_error_proportion(p=p,n=n)
    return z_alpha*se

def get_confidence_interval_proportion(confidence,n,p,alpha=None):
    if alpha==None:
        alpha=get_alpha_from_confidence(confidence)
    z_alpha=get_z_alpha(alpha=alpha/2)
    se=get_standard_error_proportion(p=p,n=n)
    margin_of_error=get_margin_of_error_proportion(z_alpha=z_alpha,se=se)
    return f"We are ${confidence*100}$% confident that the population mean is within (${p-margin_of_error}$,${p+margin_of_error}$)"


def get_standard_error_two_proportions(p_1,p_2,n_1,n_2):
    return math.sqrt(p_1*(1-p_1)/n_1+p_2*(1-p_2)/n_2)

def get_confidence_interval_difference_proportions(p_1,p_2,n_1,n_2,confidence=None,alpha=None,z_alpha=None,se=None):
    if z_alpha==None:
        if alpha==None:
            alpha=get_alpha_from_confidence(confidence)
        z_alpha=get_z_alpha(alpha=alpha/2)
    
    if se==None:
        se=get_standard_error_two_proportions(p_1=p_1,p_2=p_2,n_1=n_1,n_2=n_2)
    
    lower_bound=p_1-p_2-z_alpha*se
    upper_bound=p_1-p_2+z_alpha*se

    print_latex(f"We can be ${100*(1-alpha)}$% confident that the true difference between the population means within $({lower_bound},{upper_bound})$")
    if lower_bound<=0<=upper_bound:
        return f"Since the interval ({lower_bound},{upper_bound}) contains 0, we can say with ${100*(1-alpha)}$% confidence that there is no significant difference between the two population means"
    elif upper_bound<0:
        return f"Since ${lower_bound} \\le {upper_bound} < 0$ We can be ${100*(1-alpha)}$% confident that $\\mu_1$ is significantly less than $\\mu_2$"
    elif lower_bound>0:
        return f"Since $0<{lower_bound} \\le {upper_bound}$ We can be ${100*(1-alpha)}$% confident that $\\mu_1$ is significantly greater than $\\mu_2%"
  

def calculate():
    try:
        mean = float(entry_mean.get())
        sd = float(entry_sd.get())
        n = int(entry_n.get())
        less_than = float(entry_less_than.get())
        greater_than = float(entry_greater_than.get())
        lower_bound = float(entry_lower_bound.get())
        upper_bound = float(entry_upper_bound.get())
    except ValueError:
        messagebox.showerror("Invalid input", "Please enter valid numbers.")
        return

    se = get_standard_error_mean(sd, n)
    p_less_than = get_probability_less_than_mean(less_than, mean, sd, n)
    p_greater_than = get_probability_greater_than_mean(greater_than, mean, sd, n)
    p_between = get_probability_between_mean(lower_bound, upper_bound, mean, sd, n)

    result_text = f"Standard error: {se}\n" \
                  f"P(X < {less_than}): {p_less_than}\n" \
                  f"P(X > {greater_than}): {p_greater_than}\n" \
                  f"P({lower_bound} < X < {upper_bound}): {p_between}"

    messagebox.showinfo("Results", result_text)

window = tk.Tk()
window.title("Sampling Distribution of the Mean")

label_mean = tk.Label(window, text="Mean:")
label_mean.grid(row=0, column=0)
entry_mean = tk.Entry(window)
entry_mean.grid(row=0, column=1)

label_sd = tk.Label(window, text="Standard deviation:")
label_sd.grid(row=1, column=0)
entry_sd = tk.Entry(window)
entry_sd.grid(row=1, column=1)

label_n = tk.Label(window, text="Sample size:")
label_n.grid(row=2, column=0)
entry_n = tk.Entry(window)
entry_n.grid(row=2, column=1)

label_less_than = tk.Label(window, text="Less than:")
label_less_than.grid(row=3, column=0)
entry_less_than = tk.Entry(window)
entry_less_than.grid(row=3, column=1)

label_greater_than = tk.Label(window, text="Greater than:")
label_greater_than.grid(row=4, column=0)
entry_greater_than = tk.Entry(window)
entry_greater_than.grid(row=4, column=1)

label_lower_bound = tk.Label(window, text="Lower bound:")
label_lower_bound.grid(row=5, column=0)
entry_lower_bound = tk.Entry(window)
entry_lower_bound.grid(row=5, column=1)

label_upper_bound = tk.Label(window, text="Upper bound:")
label_upper_bound.grid(row=6, column=0)
entry_upper_bound = tk.Entry(window)
entry_upper_bound.grid(row=6, column=1)

button_calculate = tk.Button(window, text="Calculate", command=calculate)
button_calculate.grid(row=7, column=0, columnspan=2)

window.mainloop()

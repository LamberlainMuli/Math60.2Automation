# %% [markdown]
# Module 4.1: Sampling Plans and Experimental Designs

# %%
import nltk
from nltk.corpus import stopwords
from nltk.tokenize import word_tokenize

# Define the keywords for each type of sampling
judgement_keywords = ['chooses', 'selects', 'opinion', 'experience', 'judgement', 'believes', 'thinks', 'most', 'best', 'worst', 'top', 'bottom', 'likely', 'unlikely', 'risk', 'at risk', 'expert', 'subjective', 'decision', 'quality', 'specific', 'particular', 'criteria', 'standard', 'belief', 'based on']
convenience_keywords = ['easy', 'accessible', 'close', 'nearby', 'available', 'convenience', 'quick', 'fast', 'immediate', 'now', 'current', 'first', 'last', 'asks', 'his', 'her', 'their', 'volunteers', 'willing', 'accepts', 'agrees', 'convenience', 'without difficulty', 'readily available', 'without special effort', 'particular day', 'visits', 'conducts', 'website', 'store', 'call', 'helpline']
quota_keywords = ['proportion', 'first', 'quota', 'percentage', 'ratio', 'equal', 'same', 'balance', 'representative', 'demographic', 'gender', 'age', 'race', 'income', 'education', 'occupation', 'region', 'quota', 'divides', 'splits', 'categorizes', 'groups', 'segments', 'population', 'sample', 'ensuring', 'specific proportion', 'includes']

# Define the stop words
stop_words = set(stopwords.words('english'))

def classify_non_probabilistic_sampling(text):
    # Tokenize the text
    text = text.lower()
    words = word_tokenize(text)

    # Remove stop words
    words = [word for word in words if word not in stop_words]

    # Count the number of keywords for each type of sampling
    judgement_count = sum(word in words for word in judgement_keywords)
    convenience_count = sum(word in words for word in convenience_keywords)
    quota_count = sum(word in words for word in quota_keywords)

    # Classify the text based on which type has the most keywords
    # print(judgement_count, convenience_count, quota_count)
    max_count = max(judgement_count, convenience_count, quota_count)
    if max_count == judgement_count:
        return 'Judgement Sampling'
    elif max_count == convenience_count:
        return 'Convenience Sampling'
    elif max_count == quota_count:
        return 'Quota Sampling'
    else:
        return 'Unknown'
    
# Test the function
print(classify_non_probabilistic_sampling('A doctor chooses six of his most “at risk” patients to participate in a clinical trial.'))
print(classify_non_probabilistic_sampling('A student waits until Sunday night to complete a survey for his sociology class.He asks 10 of his orgmates to participate.'))
print(classify_non_probabilistic_sampling('A population consists of 40% women and 60% men. The researcher decides to choose a sample consisting of 20 women and 30 men.'))
print(classify_non_probabilistic_sampling('A professor chooses the 20 students in his class whom he thinks will most likely return his questionnaire'))

test_string='''A researcher selects a sample of 50 people from a population of 500, ensuring that the sample includes 25 men and 25 women.
quota sampling
A journalist interviews the first 10 people she encounters at a political rally.
convenience sampling
A teacher selects the 5 students who are most likely to complete an extra credit assignment.
judgement sampling
A company conducts a survey on customer satisfaction by interviewing customers who visit their store on a particular day.
convenience sampling
A pollster decides to survey a sample of 100 voters, ensuring that 60 are Democrats and 40 are Republicans.
quota sampling
A scientist selects 10 plants from his greenhouse that he believes are most likely to respond to a new fertilizer.
judgement sampling
A student conducts a survey on study habits by asking his friends who are readily available.
convenience sampling
A researcher studying obesity selects participants who she believes are most likely to provide useful data.
judgement sampling
A company decides to survey a sample of its customers, ensuring that the sample includes both users of its old product and its new product.
quota sampling
A journalist interviews people at a protest who are readily available for comments.
convenience sampling
A teacher selects a group of students who she believes are most likely to benefit from a new teaching method.
judgement sampling
A researcher studying sleep patterns decides to survey people who are readily available at a late-night diner.
convenience sampling
A pollster decides to survey a sample of voters, ensuring that the sample includes both urban and rural residents.
quota sampling
A scientist selects a group of lab rats that he believes are most likely to respond to a new drug.
judgement sampling
A student conducts a survey on food preferences by asking his dorm mates who are readily available.
convenience sampling
A researcher studying climate change selects participants who she believes are most likely to provide insightful data.
judgement sampling
A company decides to survey a sample of its customers, ensuring that the sample includes both online and in-store shoppers.
quota sampling
A journalist interviews people at a music festival who are readily available for comments.
convenience sampling
A teacher selects a group of students who she believes are most likely to excel in a math competition.
judgement sampling
A researcher studying caffeine consumption decides to survey people who are readily available at a coffee shop.
convenience sampling
A researcher selects a sample of 100 employees from a company, ensuring that the sample includes 50 managers and 50 non-managers.
quota sampling
A journalist interviews the first 10 people she encounters at a book fair.
convenience sampling
A fitness coach selects the 5 clients who are most likely to adhere to a new workout regimen.
judgement sampling
A company conducts a survey on product preferences by interviewing customers who visit their website on a particular day.
convenience sampling
A pollster decides to survey a sample of 100 citizens, ensuring that 70 are employed and 30 are unemployed.
quota sampling
A botanist selects 10 trees from a forest that he believes are most likely to show signs of a particular disease.
judgement sampling
A student conducts a survey on video game habits by asking his online friends who are readily available.
convenience sampling
A researcher studying mental health selects participants who she believes are most likely to provide comprehensive data.
judgement sampling
A company decides to survey a sample of its users, ensuring that the sample includes both Android and iOS users.
quota sampling
A journalist interviews people at a movie premiere who are readily available for comments.
convenience sampling
A coach selects a group of athletes who she believes are most likely to improve from a new training method.
judgement sampling
A researcher studying late-night eating habits decides to survey people who are readily available at a 24-hour fast food restaurant.
convenience sampling
A pollster decides to survey a sample of residents, ensuring that the sample includes both homeowners and renters.
quota sampling
A scientist selects a group of lab mice that he believes are most likely to respond to a new diet.
judgement sampling
A student conducts a survey on music preferences by asking his classmates who are readily available.
convenience sampling
A researcher studying digital literacy selects participants who she believes are most likely to provide accurate data.
judgement sampling
A company decides to survey a sample of its customers, ensuring that the sample includes both long-term and new customers.
quota sampling
A journalist interviews people at a tech conference who are readily available for comments.
convenience sampling
A teacher selects a group of students who she believes are most likely to perform well in a science fair.
judgement sampling
A researcher studying social media usage decides to survey people who are readily available at a public Wi-Fi hotspot.
convenience sampling
A researcher selects a sample of 100 students from a university, ensuring that the sample includes 50 undergraduates and 50 postgraduates.
quota sampling
A journalist interviews the first 10 people she encounters at a farmers market.
convenience sampling
A dietitian selects the 5 clients who are most likely to adhere to a new diet plan.
judgement sampling
A company conducts a survey on website usability by interviewing users who visit their website on a particular day.
convenience sampling
A pollster decides to survey a sample of 100 people, ensuring that 50 are smokers and 50 are non-smokers.
quota sampling
A zoologist selects 10 animals from a zoo that he believes are most likely to show signs of a particular behavior.
judgement sampling
A student conducts a survey on reading habits by asking his book club members who are readily available.
convenience sampling
A researcher studying physical fitness selects participants who she believes are most likely to provide detailed data.
judgement sampling
A company decides to survey a sample of its customers, ensuring that the sample includes both in-person and online shoppers.
quota sampling
A journalist interviews people at a music festival who are readily available for comments.
convenience sampling
A music teacher selects a group of students who she believes are most likely to benefit from a new teaching method.
judgement sampling
A researcher studying breakfast habits decides to survey people who are readily available at a breakfast diner.
convenience sampling
A pollster decides to survey a sample of residents, ensuring that the sample includes both pet owners and non-pet owners.
quota sampling
A scientist selects a group of lab rabbits that he believes are most likely to respond to a new treatment.
judgement sampling
A student conducts a survey on fashion preferences by asking his friends who are readily available.
convenience sampling
A researcher studying online shopping habits selects participants who she believes are most likely to provide accurate data.
judgement sampling
A company decides to survey a sample of its customers, ensuring that the sample includes both frequent and infrequent shoppers.
quota sampling
A journalist interviews people at a car show who are readily available for comments.
convenience sampling
A sports coach selects a group of athletes who she believes are most likely to excel in a tournament.
judgement sampling
A researcher studying coffee consumption decides to survey people who are readily available at a coffee shop.
convenience sampling
A researcher selects a sample of 100 households from a city, ensuring that the sample includes 50 with children and 50 without children.
quota sampling
A journalist interviews the first 10 people she encounters at a city park.
convenience sampling
A psychologist selects the 5 clients who are most likely to benefit from a new therapy technique.
judgement sampling
A company conducts a survey on customer service by interviewing customers who call their helpline on a particular day.
convenience sampling
A pollster decides to survey a sample of 100 people, ensuring that 70 are meat-eaters and 30 are vegetarians.
quota sampling
A biologist selects 10 birds from a sanctuary that he believes are most likely to exhibit a particular behavior.
judgement sampling
A student conducts a survey on movie preferences by asking his film class peers who are readily available.
convenience sampling
A researcher studying dietary habits selects participants who she believes are most likely to provide comprehensive data.
judgement sampling
A company decides to surveya sample of its users, ensuring that the sample includes both free and premium subscribers.
quota sampling
A journalist interviews people at a food festival who are readily available for comments.
convenience sampling
A dance instructor selects a group of students who she believes are most likely to excel in a dance competition.
judgement sampling
A researcher studying late-night TV viewing habits decides to survey people who are readily available at a late-night diner.
convenience sampling
A pollster decides to survey a sample of residents, ensuring that the sample includes both car owners and non-car owners.
quota sampling
A scientist selects a group of lab fish that he believes are most likely to respond to a new environment.
judgement sampling
A student conducts a survey on travel preferences by asking his classmates who are readily available.
convenience sampling
A researcher studying work-from-home habits selects participants who she believes are most likely to provide accurate data.
judgement sampling
A company decides to survey a sample of its customers, ensuring that the sample includes both domestic and international customers.
quota sampling
A journalist interviews people at a tech expo who are readily available for comments.
convenience sampling
A chess coach selects a group of players who she believes are most likely to perform well in a chess tournament.
judgement sampling
A researcher studying smartphone usage decides to survey people who are readily available at a public transportation station.
convenience sampling
A researcher decides to interview the first 10 people who walk into a clinic, ensuring that 5 are men and 5 are women.
quota sampling
A journalist selects the first 10 people she encounters at a concert, based on her belief that they will provide the most interesting stories.
judgement sampling
A teacher asks the students who stay after class to participate in a study, ensuring that the group includes both high-performing and low-performing students.
quota sampling
A company surveys the customers who make a purchase on a particular day, based on the assumption that these customers will provide the most relevant feedback.
judgement sampling
A student decides to survey his classmates who are present in the library, ensuring that the sample includes both science and arts students.
quota sampling
A scientist chooses to study the first 10 plants he encounters on a field trip, based on his belief that these plants are most likely to exhibit a certain trait.
judgement sampling
A researcher decides to survey the people who are currently in a coffee shop, ensuring that the sample includes both coffee drinkers and non-coffee drinkers.
quota sampling
A journalist interviews the people who are readily available at a city park, based on her belief that they will provide the most diverse viewpoints.
judgement sampling
A company decides to survey its employees who are present on a particular day, ensuring that the sample includes both managers and non-managers.
quota sampling
A student asks his friends who are currently online to participate in a survey, based on his belief that they will provide the most honest responses.
judgement sampling
A researcher decides to survey the people who are currently at a fitness center, ensuring that the sample includes both regular members and first-time visitors.
quota sampling
A journalist interviews the people who are readily available at a political rally, based on her belief that they will provide the most insightful comments.
judgement sampling
A company decides to survey its customers who make a purchase on a particular day, ensuring that the sample includes both first-time and repeat customers.
quota sampling
A student asks his classmates who are currently in the study room to participate in a survey, based on his belief that they will provide the most accurate data.
judgement sampling
A researcher decides to survey the people who are currently at a museum, ensuring that the sample includes both locals and tourists.
quota sampling
A journalist interviews the people who are readily available at a sports event, based on her belief that they will provide the most engaging stories.
judgement sampling
A company decides to survey its users who are online on a particular day, ensuring that the sample includes both free and premium users.
quota sampling
A student asks his friends who are currently at a party to participate in a survey, based on his belief that they will provide the most candid responses.
judgement sampling
A researcher decides to survey the people who are currently at a shopping mall, ensuring that the sample includes both shoppers and non-shoppers.
quota sampling
A journalist interviews the people who are readily available at a book fair, based on her belief that they will provide the most interesting perspectives.
judgement sampling'''

test_answers=(test_string.split('\n'))
tests=[]
answers=[]
# print(test_answers)
# print(len(test_answers))
for i in range(0,len(test_answers),2):
    tests.append(test_answers[i])
    answers.append(test_answers[i+1])

print("Test has started: \n")
num=len(tests)
correct=0
for test,answer in zip(tests,answers):
    # print("Statement:\t", test)
    # print("Model Answer: \t\t",classify_non_probabilistic_sampling(test).lower())
    # print("Correct Answer: \t",answer)
    if classify_non_probabilistic_sampling(test).strip().lower()!=answer.lower().strip():
        print("Statement:\t", test)
        print("Model Answer: \t\t",classify_non_probabilistic_sampling(test).lower())
        print("Correct Answer: \t",answer)
        print('Wrong!')

    else:
        correct+=1
        # print('Correct!')
print(f"The model has {correct} correct answers out of {num} tests. The accuracy is {correct/num*100}%.")

# %%
import spacy

# Load the spaCy model
nlp = spacy.load('en_core_web_sm')

def classify_probabilistic_sampling(statement):
    # Parse the statement using spaCy
    doc = nlp(statement)

    # Extract the named entities
    named_entities = [(ent.text, ent.label_) for ent in doc.ents]
    print(named_entities)
    # Initialize variables to store the population, data, and sample
    population = None
    data = None
    sample = None

    # Loop over the named entities to identify the population and data
    for entity, label in named_entities:
        if label == 'GPE':  # Geo-political entity, e.g., a city or country
            population = entity
        elif label == 'ORG':  # Organization, e.g., a company or institution
            data = entity

    # Loop over the tokens in the document to identify the sample
    for token in doc:
        if token.dep_ == 'dobj':  # Direct object, e.g., the object that the action is being done to
            sample = doc[token.i :].text  # Get the text from the current token to the end of the sentence

    # Use the relationships between the population, data, and sample to classify the sampling method
    if population and data and sample:
        if 'within each' in sample:
            return 'Stratified Sampling'
        elif 'each' in sample:
            return 'Simple Random Sampling'
        elif 'all' in sample:
            return 'Cluster Sampling'
        elif 'every' in sample:
            return 'Systematic Sampling'

    return 'Unknown'



print(classify_probabilistic_sampling('Divide Quezon City into barangays and take a simple random sample of people within each barangay.'))
print(classify_probabilistic_sampling('Label the different barangays of Quezon City from 1 to N, and use a random number generator to take a sample of 10 barangays.'))
print(classify_probabilistic_sampling('Divide Quezon City into city blocks, choose a simple random sample of 10 city blocks, and interview all who live there.'))
print(classify_probabilistic_sampling('Choose an entry at random from the city register, and select every 50th resident thereafter.'))


# %% [markdown]
# Module 4.2: Sampling Distribution of the Mean and the CLT

# %%
import math
def get_standard_error_mean(sd,n):
    return sd/math.sqrt(n)

# %%
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



# %% [markdown]
# When is the sample size large enough to use the CLT?

# %%
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

test_clt_mean(sample_size=30)

# %%
def get_sampling_distribution_mean_variable(mean,sd,n):
    print_latex(f"Sampling distribution is: $N({mean},{get_standard_error_mean(sd,n)})$")

get_sampling_distribution_mean_variable(mean=5.9,sd=0.7,n=16)

# %% [markdown]
# Module 4.3: Sampling Distribution of the Sample Proportion
# 
# $X \approx B(n, p)$, mean of $X$ is $np$ and the standard deviation of $X$ is $\sqrt{np(1-p)}$

# %%

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

# test_clt_proportion(p=0.9,n=1)

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

# %% [markdown]
# Module 5.1: Estimators and Point Estimation

# %%
def get_sample_standard_deviation(sample):
    return math.sqrt(sum([(xi-sum(sample)/len(sample))**2 for xi in sample])/(len(sample)-1))

def get_point_estimate_mean(sample=[],n=1,total=0):
    if sample==[]:
        return total/n
    else:
        return sum(sample)/len(sample)

print(get_point_estimate_mean(n=50,total=49000))
print(get_standard_error_mean(n=50,sd=105))

def estimate_mean(mean,n,sd):
    print_latex(f"Point estimate of the population mean is: $\\hat \\mu \\approx N({mean},{get_standard_error_mean(sd,n)})$")
    return mean,get_standard_error_mean(sd,n)
mean_1,se_1=estimate_mean(mean=980,sd=105,n=50)
mean_2,se_2=estimate_mean(mean=975.8,sd=108,n=100)

def compare_point_estimates_mean(mean_1,se_1,mean_2,se_2):
    if se_1<=se_2:
        print_latex(f"Since the standard error of the first point estimate is smaller than the second point estimate, the first point estimate is more precise")
    else:
        print_latex(f"Since the standard error of the first point estimate is larger than the second point estimate, the second point estimate is more precise")

compare_point_estimates_mean(mean_1,se_1,mean_2,se_2)

# %%
def get_point_estimate_proportion(n,p):
    return p

def estimate_proportion(n,p):
    print_latex(f"Point estimate of the population proportion is: $\\hat p \\approx N({p},{get_standard_error_proportion(n,p)})$")
    return p,get_standard_error_proportion(n,p)

p_1,se_1=estimate_proportion(n=50,p=0.2)


# %% [markdown]
# Module 5.2: Confidence Intervals for Population Means

# %%
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

print(get_confidence_interval_mean(confidence=0.95,n=50,mean=8.0,sd=0.3))
print(get_confidence_interval_mean(confidence=0.99,n=50,mean=67.3,sd=0.64))


print(get_margin_of_error_mean(sd=15,n=50,confidence=0.95))


# %%
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
    

print_latex(get_confidence_interval_difference_means(mean_1=62428,sd_1=12500,mean_2=57762,sd_2=13330,confidence=0.95,n_1=50,n_2=50))
print_latex(get_confidence_interval_difference_means(mean_1=220.2,sd_1=14.3,n_1=40,mean_2=236.8,sd_2=18.8,n_2=40,confidence=0.99))


# %% [markdown]
# Module 5.3: Confidence Intervals for Population Proportions

# %%
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

print_latex(get_confidence_interval_proportion(confidence=0.95,n=150,p=134/150,))
print_latex(get_confidence_interval_proportion(confidence=0.98,n=2625,p=735/2625,))


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
    

print_latex(get_confidence_interval_difference_proportions(p_1=13/100,p_2=6/100,n_1=100,n_2=100,confidence=0.98))
print_latex(get_confidence_interval_difference_proportions(p_1=0.93,n_1=200,p_2=0.96,n_2=450,confidence=0.99))




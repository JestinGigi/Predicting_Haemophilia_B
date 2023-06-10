# %% [markdown]
# Importing Libraries

# %%
import pandas as pd # for data manipulation and analysis
import numpy as np # linear algebra

# For visualization 
import matplotlib.pyplot as plt 
import seaborn as sns 
sns.set_style("whitegrid")

# Ignoring warnings
import warnings
warnings.filterwarnings('ignore')

# Generating random values
import random as rd 

# For regex search
import re

#Importing libraries ML models
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB as GNB
from sklearn.ensemble import RandomForestClassifier
from sklearn.ensemble import GradientBoostingClassifier
from sklearn.metrics import accuracy_score

#Importing dataset splitting library
from sklearn.model_selection import train_test_split

#Import Joblib
import joblib

# %% [markdown]
# Exploratory Data analysis

# %%
#Read Dataset file 
df = pd.read_csv('Dataset/Point_variant_factorIX.csv')

# %%
df.sample(10)

# %%
df.shape

# %%
df.info()

# %%
df.describe(include='all')

# %%
#Plot non-null values of each features using Horizontal Barplot
features, count = list(df.columns), list(df.notnull().sum())
plt.figure(dpi=150)
plt.barh(y=features, width=count)
plt.title("Non Null values of features")
plt.show()

# %% [markdown]
# Data Preprocessing

# %%
# Change all values to lowercase
df = df.applymap(lambda s : s.lower() if type(s) == str else s)

# %%
# Remove all the irrevalent columns
new_df = df.drop(columns=[i[0] for i in zip(df.columns, range(len(df))) if i[1] not in [3, 7, 8, 9, 10, 11, 12, 15]])

# Drop rows with null values
new_df = new_df.dropna()

# %%
# Drop ambiguious severity values
new_df = new_df[new_df.Severity.isin(['-']) == False]

# %%
# Drop "-" values from FIX:C% and FIX:Ag% 
new_df = new_df[(new_df['FIX:Ag %'] != '-') & (new_df['FIX:C %'] != '-')]

# %%
# Drop duplicate values
new_df = new_df.drop_duplicates()

# %%
#Correct values in FIX:C%
# max(FIX:C%) = 100%
# min(FIX:C%) = 0%
for index, row in new_df.iterrows(): #iternating over rows
    less_match = re.search(r'(?<=<)\d*\.*\d*', row['FIX:C %']) # searching for <100/<.34/<2.34 type values
    to_match = re.search(r'(\d*\.*\d*)\s*(to|-)\s*(\d*\.*\d*)', row['FIX:C %']) # searching for 10 to 5/10to14/10-15 type values 
    if less_match:
        new_df.at[index, 'FIX:C %'] = float(less_match[0])*rd.random() # converting to floating value in the range(0,less_match[0])
    if to_match:
        new_df.at[index, 'FIX:C %'] = rd.uniform(float(to_match[1]), float(to_match[3])) # converting to floating value in the range(less_match[0],less_match[1])


#Correct values in FIX:Ag%
# max(FIX:Ag%) = 200%
# min(FIX:Ag%) = 0%
for index, row in new_df.iterrows(): #iternating over rows
    less_match = re.search(r'(?<=<)\d*\.*\d*', row['FIX:Ag %']) # searching for <100/<.34/<2.34 type values
    to_match = re.search(r'(\d*\.*\d*)\s*(to|-)\s*(\d*\.*\d*)', row['FIX:Ag %']) # searching for 10 to 5/10to14/10-15 type values 
    grt_match = re.search(r'(?<=>)\d*\.*\d*', row['FIX:Ag %']) # Searcing for >100 values
    if less_match:
        new_df.at[index, 'FIX:Ag %'] = float(less_match[0])*rd.random() 
    if to_match:
        new_df.at[index, 'FIX:Ag %'] = rd.uniform(float(to_match[1]), float(to_match[3])) 
    if grt_match:
        new_df.at[index, 'FIX:Ag %'] = rd.uniform(float(grt_match[0]), 200) # converting to floating value in range(grt_match[0], 200)

# %%
# Convert values in FIX:C % and FIX:Ag % to float
new_df['FIX:C %'] = new_df['FIX:C %'].astype(float)
new_df['FIX:Ag %'] = new_df['FIX:Ag %'].astype(float)

# %%
#Write unique values and their count of each feature to a file
with open('unique_val_count.txt', 'w') as f:
    for col in new_df.columns:
        f.write(col+'\n')
        for i in zip(list(new_df[col].value_counts().index), list(new_df[col].value_counts().values)):
            f.write(str(i))
            f.write("\n")
        f.write("\n")
f.close()

# %%
#Separate "Variant" into 3 columns

#Initializing lists
act_nucleo = [] 
mut_nucleo = []
nucleo_pos = []

for index, row in new_df.iterrows(): # iterating over rows
    match = re.search(r'(\d+)(\w)>{1}(\w)', row['Variant']) # regex to extract actual, mutated nucleotide and mutation position
    if match:
        nucleo_pos.append(match[1]) # group 1: mutation position of nucleotide
        act_nucleo.append(match[2]) # group 2: actual nucleotide
        mut_nucleo.append(match[3]) # group 3: mutated nucleotide

# dropping existing Variant column
new_df = new_df.drop('Variant', axis=1) 

#Adding new Variant columns with its values
new_df['act_nucleo'] = act_nucleo 
new_df['nucleo_pos'] = nucleo_pos
new_df['mut_nucleo'] = mut_nucleo


# %%
#Separate "Protein change" into 3 columns
#Initializing lists
act_amino = [] 
mut_amino = []
amino_pos = []

for index, row in new_df.iterrows(): # iterating over rows
    match = re.search(r'([a-z]+)(\d+)([a-z]+|\*)', row['Protein Change']) # regex to extract actual, mutated nucleotide and mutation position
    if match:
        act_amino.append(match[1]) # group 1: mutation position of amino acid
        amino_pos.append(match[2]) # group 2: actual amino acid
        # group 3: mutated amino acide
        if match[3] == "*":
            mut_amino.append("stp_cdn")    
        else:
            mut_amino.append(match[3])  

# dropping existing Variant column
new_df = new_df.drop('Protein Change', axis=1) 

#Adding new Variant columns with its values
new_df['act_amino'] = act_amino 
new_df['amino_pos'] = amino_pos
new_df['mut_amino'] = mut_amino

# %%
# Move severity column to end
cols = list(new_df.columns)
cols.insert(len(cols), cols.pop(cols.index('Severity')))
new_df = new_df.loc[:, cols]

# %%
#Converting to int
cols = ['nucleo_pos', 'amino_pos']
for i in cols:
    new_df[i] = new_df[i].astype(int)

# %%
new_df.info()

# %%
sns.countplot(data=new_df, x='Effect', hue='Severity')

# %%
plt.figure(figsize=(15,7))
sns.countplot(data=new_df, x='Domain', hue='Severity')

# %%
sns.countplot(data=new_df, x='Location in gene', hue='Severity')

# %%
new_df.info()

# %%
sns.countplot(data=new_df, x='Severity')

# %%
sns.countplot(data=new_df, x='mut_nucleo', hue='Severity')

# %%
plt.figure(figsize=(15,6))
sns.countplot(data=new_df, x='mut_amino', hue='Severity')

# %% [markdown]
# Data Visualization

# %%
#No.of cases vs FIX: C% 
plt.figure(figsize=(10, 10))
sns.displot(new_df['FIX:C %'], color='b')
plt.title('No. of Cases Vs FIX:C %')
plt.ylabel('Number of cases')
plt.show()

# %%
#No.of cases vs FIX: Ag% 
plt.figure(figsize=(10, 10))
sns.displot(new_df['FIX:Ag %'], color='b')
plt.title('No. of Cases Vs FIX:Ag %')
plt.ylabel('Number of cases')
plt.show()

# %%
new_df.info()

# %%
plt.figure(figsize=(10,6))
sns.boxplot(x=new_df['FIX:C %'])
plt.title('Distribution of FIX:C %')
plt.show()

# %%
Q1 = new_df['FIX:C %'].quantile(0.25)
Q3 = new_df['FIX:C %'].quantile(0.75)
IQR = Q3 - Q1
print(IQR)

# %%
new_df[(new_df['FIX:C %'] < Q1 - 1.5*IQR) | (new_df['FIX:C %'] > Q3 + 1.5*IQR)]

# %%
plt.figure(figsize=(10,6))
sns.boxplot(x=new_df['FIX:Ag %'])
plt.title('Distribution of FIX: Ag %')
plt.show()

# %%
Q1 = new_df['FIX:Ag %'].quantile(.25)
Q3 = new_df['FIX:Ag %'].quantile(.75)
IQR = Q3 - Q1

# %%
new_df[(new_df['FIX:Ag %'] < Q1 - 1.5*IQR) | (new_df['FIX:Ag %'] > Q3 + 1.5*IQR)]

# %%
new_df.info()

# %%
n_df = new_df.iloc[:, :-1]
t_df = new_df.iloc[:, -1]

# %%
from sklearn.preprocessing import OrdinalEncoder
from sklearn.preprocessing import LabelEncoder
from sklearn.compose import ColumnTransformer
from sklearn.preprocessing import StandardScaler

# Label Encode severity
lb = LabelEncoder()
t_df = lb.fit_transform(t_df)

# Feature transformation
categorial_cols = [i for i in n_df.columns if n_df.dtypes[i] in ['O']]
numerical_cols = [i for i in n_df.columns if n_df.dtypes[i] in ['float64', 'int64']]

feature_cols = categorial_cols + numerical_cols

transformer = [
    ('encode', OrdinalEncoder(), categorial_cols),
    ('scale', StandardScaler(), numerical_cols)
]

colTrans = ColumnTransformer(transformer, remainder='passthrough')
n_df_raw = colTrans.fit_transform(n_df)

n_df = pd.DataFrame(
    n_df_raw, columns=colTrans.get_feature_names_out())


# %%
# Adding severity
n_df['Severity'] = t_df

# %%
sns.heatmap(n_df.corr())

# %% [markdown]
# Model Training

# %%
# Split data into target and features
x = n_df.iloc[:, :-1]
y = n_df.iloc[:, -1]

X_train, X_test, y_train, y_test = train_test_split(x, y, test_size=.25, random_state=42)

# correcting imbalance dataset
from imblearn.over_sampling import SMOTE
smote = SMOTE(random_state=42)
X_train, y_train = smote.fit_resample(X_train, y_train)

# %%
#object creation
nb_classifier = GNB()
log_reg = LogisticRegression(random_state=42)
rf = RandomForestClassifier(n_estimators=100)
gb = GradientBoostingClassifier(random_state=42)

# %%
#Training
log_reg.fit(X_train, y_train)
nb_classifier.fit(X_train, y_train)
rf.fit(X_train, y_train)
gb.fit(X_train, y_train)

# %%
# Calculate the accuracy
print(f"Logistic regression: {accuracy_score(y_test, log_reg.predict(X_test))}")
print(f"Naive Bayes: {accuracy_score(y_test, nb_classifier.predict(X_test))}")
print(f"Random Forest: {accuracy_score(y_test, rf.predict(X_test))}")
print(f"Gradient boosting: {accuracy_score(y_test, gb.predict(X_test))}")

# %%
n_df.info()

# %%
joblib.dump(gb, 'gradient_model.joblib')
joblib.dump(colTrans, 'colTrans.joblib')
joblib.dump(lb, 'labelEncoder.joblib')

# %%
user_input = [['silent', 'linker', 'exon 3', 'c', 'g', 'leu', 'arg', 1.5, 90, 328, 355]]
user_df = pd.DataFrame(user_input, columns=feature_cols)
user_df = colTrans.transform(user_df)
print(lb.inverse_transform(gb.predict(user_df)))



# Haemophilia B Severity Prediction with Machine Learning Models

This project focuses on predicting the severity of Haemophilia B using various machine learning models, including Logistic Regression, Gaussian Naive Bayes (GNB), Gradient Boosting, and Random Forest Classifier. Additionally, a simple web application is created using Flask to showcase the predictions. The dataset used for this project is the F9 gene mutation dataset available at the public repository of EAHAD.

## Dataset
The F9 gene mutation dataset from the EAHAD public repository is utilized in this project. This dataset contains relevant features related to Haemophilia B and the severity levels associated with them.

## Exploratory Analysis and Data Preprocessing
Before training the machine learning models, a thorough exploratory analysis is performed on the dataset. This analysis involves understanding the data, identifying missing values, handling outliers, and visualizing the distributions and relationships between different variables. Proper data preprocessing techniques are applied to ensure the dataset is ready for model training.

## Machine Learning Models
The following machine learning models are employed for predicting the severity of Haemophilia B:

1. Logistic Regression: This model applies logistic function to estimate the probability of the severity levels based on the input features.
2. Gaussian Naive Bayes (GNB): GNB assumes that the features are independently distributed and calculates the posterior probability of the severity levels using Bayes' theorem.
3. Gradient Boosting: This ensemble model combines multiple weak learners (decision trees) to create a strong predictive model by minimizing the errors of the previous models.
4. Random Forest Classifier: Random Forest creates an ensemble of decision trees and provides predictions based on the majority vote of the individual trees.

## Web Application using Flask
To demonstrate the predictions made by the machine learning models, a simple web application is developed using Flask, a Python web framework. The application allows users to input the relevant features, and the trained models provide predictions regarding the severity of Haemophilia B based on those inputs.

The web application provides an interactive and user-friendly interface to facilitate the prediction process. Users can easily access and utilize the model's functionality.

## Usage
1. Move to repository
2. Install all dependencies using pip
```bash
pip install -r requirements.txt
```
3. Run `Haemophilia_prediction_B.ipynb` 
4. Run the flask app
```
flask --app app run
```
5. Visit the web application on `http://127.0.0.1:5000`

## Conclusion
This project utilizes machine learning models, including Logistic Regression, GNB, Gradient Boosting, and Random Forest Classifier, to predict the severity of Haemophilia B. By employing proper exploratory analysis and data preprocessing techniques, the F9 gene mutation dataset from the EAHAD public repository is prepared for model training. Additionally, a simple web application using Flask is created to showcase the predictions, enhancing accessibility and usability.
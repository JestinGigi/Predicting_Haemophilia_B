from flask import Flask, render_template, request, jsonify
import sklearn
import pandas as pd
import joblib

app = Flask(__name__)
gb = joblib.load('gradient_model.joblib')
colTrans = joblib.load('colTrans.joblib')
lb = joblib.load('labelEncoder.joblib')


@app.route("/")
def index():
    return render_template("index.html")

@app.route("/predict", methods=["POST"])
def predict():
    effect = request.json["effect"]
    domain = request.json["domain"]
    gene_location = str(request.json["gene_location"])
    act_nucleo = str(request.json["act_nucleo"])
    mut_nucleo = str(request.json["mut_nucleo"])
    act_amino = str(request.json["act_amino"])
    mut_amino = str(request.json["mut_amino"])
    fix_c = str(request.json["fix_c"])
    fix_ag = str(request.json["fix_ag"])
    nucleo_pos = str(request.json["nucleo_pos"])
    amino_pos = str(request.json["amino_pos"])
    
    #Creating a 2d array of user input
    user_input = [[effect, domain, gene_location, act_nucleo, mut_nucleo, act_amino, mut_amino, fix_c, fix_ag, nucleo_pos, amino_pos]]
    user_df = pd.DataFrame(user_input, columns=['Effect', 'Domain', 'Location in gene', 'act_nucleo', 'mut_nucleo', 'act_amino', 'mut_amino', 'FIX:C %', 'FIX:Ag %', 'nucleo_pos', 'amino_pos'])
    user_df = colTrans.transform(user_df)
    prediction = lb.inverse_transform(gb.predict(user_df))
    preds = prediction[0]
    return {'severity': preds}

@app.route('/suggest/<pred>', methods=["GET"])
def suggest(pred):
    if pred == 'mild':
        return render_template('mild.html')
    elif pred == 'moderate':
        return render_template('moderate.html')
    else:
        return render_template('severe.html')

if __name__ == '__main__':
    app.run()
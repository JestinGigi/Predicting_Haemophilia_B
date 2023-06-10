document.getElementById("prediction-form").addEventListener("submit", function(event) {
    event.preventDefault(); // Prevent form submission
  
    var effect = document.getElementById("effect").value;
    var domain = document.getElementById("domain").value;
    var gene_location = document.getElementById("gene_location").value;
    var act_nucleo = document.getElementById("act_nucleo").value;
    var mut_nucleo = document.getElementById("mut_nucleo").value;
    var act_amino = document.getElementById("act_amino").value;
    var mut_amino = document.getElementById("mut_amino").value;
    var fix_c = document.getElementById("fix_c").value;
    var fix_ag = document.getElementById("fix_ag").value;
    var nucleo_pos = document.getElementById("nucleo_pos").value;
    var amino_pos = document.getElementById("amino_pos").value;

  
    // Send the data to the server for prediction
    fetch("/predict", {
      method: "POST",
      headers: {
        "Content-Type": "application/json"
      },
      body: JSON.stringify({
        effect: effect,
        domain: domain,
        gene_location: gene_location,
        act_nucleo: act_nucleo,
        mut_nucleo: mut_nucleo,
        act_amino: act_amino,
        mut_amino: mut_amino,
        fix_c: fix_c,
        fix_ag: fix_ag,
        nucleo_pos: nucleo_pos,
        amino_pos: amino_pos
      })
    })
    .then(response => response.json())
    .then(data => {
      // Display the prediction result
      document.getElementById("result").textContent = data.severity;
      document.getElementById('link').href = "/suggest/" + data.severity;
      document.getElementById('output').textContent = "Predicted Severity: " + data.severity;
    })
    .then((json) => console.log(json));
  });
  
import sys, requests

def mutalyzer_variantConversion( hgvs ):

    # Set up the parameters we want to pass to the API.
    # This is the variant HGVS nomenclature and genome build hg19 [other values: hg18, mm10]
    parameters = {"variant": hgvs, "build": "hg19"}

    # Make a get request with the parameters.
    response = requests.get("https://mutalyzer.nl/json/numberConversion?", params=parameters)
    data = response.json()

    for p in data: print(p)

    return;

mutalyzer_variantConversion(str(sys.argv[1]))
